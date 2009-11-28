/*
 *             Automatically Tuned Linear Algebra Software v3.8.3
 *                    (C) Copyright 2000 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "atlas_type.h"

#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )
#define Mmin(x, y) ( (x) > (y) ? (y) : (x) )
#define Mmax(x, y) ( (x) > (y) ? (x) : (y) )
#define ATL_MVIsMM(iflag_) ( ((iflag_) | 8) == (iflag_) )
#define ATL_MVIsAxpy(iflag_) ( ((iflag_) | 16) == (iflag_) )
#define ATL_MVIsDot(iflag_) ( !ATL_MVIsMM(iflag_) && !ATL_MVIsAxpy(iflag_) )
#define ATL_MVNoBlock(iflag_) ( ((iflag_) | 32) == (iflag_) )
int GetL1CacheSize()
{
   FILE *L1f;
   int L1Size;

   L1f = fopen("res/L1CacheSize", "r");
   if (!L1f)
   {
      assert(system("make res/L1CacheSize\n") == 0);
      L1f = fopen("res/L1CacheSize", "r");
      assert(L1f != NULL);
   }
   fscanf(L1f, "%d", &L1Size);
   fclose(L1f);
   fprintf(stderr, "\n      Read in L1 Cache size as = %dKB.\n",L1Size);
   return(L1Size);
}

int GetL1Elts(char pre, double l1mul)
{
   int l1;
   l1 = 1024 * GetL1CacheSize();
   if (pre == 's') l1 /= ATL_ssize;
   else if (pre == 'd') l1 /= ATL_dsize;
   else if (pre == 'c') l1 /= ATL_csize;
   else if (pre == 'z') l1 /= ATL_zsize;
   l1 = l1mul * l1;
   return(l1);
}
void emit_mvhead(char pre, int l1mul)
{
   char fnam[64];
   int l1;
   FILE *fp;

   l1 = GetL1Elts(pre, 0.01*l1mul);
   sprintf(fnam, "atlas_%cmv.h", pre);
   fp = fopen(fnam, "w");
   fprintf(fp, "#ifndef ATLAS_%cMV_H\n", toupper(pre));
   fprintf(fp, "#define ATLAS_%cMV_H\n\n", toupper(pre));

   fprintf(fp, "#define ATL_L1mvelts %d\n", l1);
   fprintf(fp, "#include \"atlas_%cmvN.h\"\n", pre);
   fprintf(fp, "#include \"atlas_%cmvT.h\"\n", pre);
   fprintf(fp, "#include \"atlas_%cmvS.h\"\n", pre);

   fprintf(fp, "\n#endif\n");
   fclose(fp);
}

void emit_transhead(char pre, char TA, int flag, int mu, int nu, int l1mul)
{
   char fnam[128];
   FILE *fp;
   int l1, mb;

   l1 = GetL1Elts(pre, 0.01*l1mul);
   sprintf(fnam, "atlas_%cmv%c.h", pre, TA);
   fp = fopen(fnam, "w");
   assert(fp);
   fprintf(fp, "#ifndef ATLAS_MV%c_H\n", TA);
   fprintf(fp, "#define ATLAS_MV%c_H\n\n", TA);

   fprintf(fp, "#include \"atlas_misc.h\"\n");
   if (ATL_MVIsMM(flag))
   {
      fprintf(fp, "#include \"atlas_%cNCmm.h\"\n", pre);
      fprintf(fp, "#include \"atlas_lvl3.h\"\n");
   }
   fprintf(fp, "\n");

   if (!ATL_MVIsMM(flag))
   {
      fprintf(fp, "#define ATL_mv%cMU %d\n", TA, mu);
      fprintf(fp, "#define ATL_mv%cNU %d\n", TA, nu);
   }
   else
   {
      fprintf(fp, "#define ATL_mv%cMU ATL_mmMU\n", TA);
      fprintf(fp, "#define ATL_mv%cNU NB\n", TA);
      fprintf(fp, "#define ATL_mv%cCallsGemm\n", TA);
   }
   if (ATL_MVNoBlock(flag)) fprintf(fp, "#define ATL_mv%cNOBLOCK\n", TA);
   fprintf(fp, "#ifndef ATL_L1mvelts\n");
   fprintf(fp, "   #define ATL_L1mvelts ((3*ATL_L1elts)>>2)\n");
   fprintf(fp, "#endif\n");

   if (TA == 'S') /* symm/hemm stuff */
   {
      fprintf(fp, "#ifndef ATL_mvNNU\n");
      fprintf(fp, "   #include \"atlas_%cmvN.h\"\n", pre);
      fprintf(fp, "#endif\n");
      fprintf(fp, "#ifndef ATL_mvTNU\n");
      fprintf(fp, "   #include \"atlas_%cmvT.h\"\n", pre);
      fprintf(fp, "#endif\n");
      if (ATL_MVIsMM(flag))
      {
         fprintf(fp,
         "#define ATL_GetPartSYMV(A_, lda_, mb_, nb_) *(mb_) = *(nb_) = NB \n");
      }
      else if (ATL_MVIsAxpy(flag))
      {
         fprintf(fp,
                 "#define ATL_GetPartSYMV(A_, lda_, mb_, nb_) \\\n{ \\\n");
         fprintf(fp, "   *(nb_) = ATL_mvSMU; \\\n");
         mb = (l1/2 - mu)/(mu+1);
         if (mb > 2*nu) mb = (mb/(2*nu))*nu;
         else mb = nu;
         fprintf(fp, "   *(mb_) = %d; \\\n}\n", mb);
      }
      else
      {
         fprintf(fp,"#define ATL_GetPartSYMV(A_, lda_, mb_, nb_) \\\n{ \\\n");
         fprintf(fp, "   for(*(nb_)=ATL_mvNNU; (*(nb_))*(*(nb_))+((*(nb_))<<2) + ATL_mvTMU*(*(nb_)+2) < ATL_L1mvelts; *(nb_) += ATL_mvNNU); \\\n");
         fprintf(fp,
               "   *(mb_) = Mmax(ATL_mvTNU, ATL_mvNMU) + *(nb_); \\\n");
         fprintf(fp,
         "   *(mb_) = ( ATL_L1mvelts - ((*(mb_))<<1) ) / (*(mb_) + 2); \\\n");
         fprintf(fp,
   "   if (*(mb_) > ATL_mvTNU) *(mb_) = (*(mb_)/ATL_mvTNU)*ATL_mvTNU; \\\n");
         fprintf(fp, "   else *(mb_) = ATL_mvTNU; \\\n");
         fprintf(fp, "}\n");
      }
   }
   else if (TA == 'N')
   {
      if (ATL_MVIsAxpy(flag)) fprintf(fp, "#define ATL_AXPYMV\n");
      else if (nu < 32)
      {
         fprintf(fp, "#ifndef ATL_mvpagesize\n");
         fprintf(fp, "   #define ATL_mvpagesize ATL_DivBySize(4096)\n");
         fprintf(fp, "#endif\n");
         fprintf(fp, "#ifndef ATL_mvntlb\n");
         fprintf(fp, "   #define ATL_mvntlb 56\n");
         fprintf(fp, "#endif\n");
      }

      fprintf(fp, "\n#define ATL_GetPartMVN(A_, lda_, mb_, nb_) \\\n{ \\\n");
      if (ATL_MVIsAxpy(flag))
      {

         fprintf(fp,
"   *(mb_) = (ATL_L1mvelts - (ATL_mvNNU<<1)) / ((ATL_mvNNU<<1)+1); \\\n");
         fprintf(fp,
"   if (*(mb_) > ATL_mvNMU) *(mb_) = ATL_mvNMU*( *(mb_)/ATL_mvNMU ); \\\n");
         fprintf(fp, "   else *(mb_) = ATL_mvNMU; \\\n");
         fprintf(fp, "   *(nb_) = ATL_mvNNU; \\\n");
      }
      else if (nu >= 32)
      {
         fprintf(fp, "   *(nb_) = ATL_mvNNU; \\\n");
         fprintf(fp,
"   *(mb_) = (ATL_L1mvelts - ATL_mvNNU - (ATL_mvNNU+1)*ATL_mvNMU) / (ATL_mvNNU+1); \\\n");
         fprintf(fp,
"   if (*(mb_) > ATL_mvNMU) *(mb_) = (*(mb_) / ATL_mvNMU)*ATL_mvNMU; \\\n");
         fprintf(fp, "   else *(mb_) = ATL_mvNMU; \\\n");
      }
      else if (ATL_MVIsMM(flag)) /* gemv calls gemm */
      {
         fprintf(fp,
            "   *(mb_) = (ATL_L1mvelts - (ATL_mmMU+1)*NB) / (NB + 2); \\\n");
         fprintf(fp,
         "   if (*(mb_) > NB) *(mb_) = ATL_MulByNB(ATL_DivByNB(*(mb_))); \\\n");
         fprintf(fp,
         "   else if (*(mb_) < ATL_mmMU) *(mb_) = NB; \\\n");
         fprintf(fp,
         "   else  *(mb_) = (*(mb_) / ATL_mmMU) * ATL_mmMU; \\\n");
         fprintf(fp, "   *(nb_) = NB; \\\n");
      }
      else
      {
         assert(mu && nu);
         fprintf(fp,
"   *(nb_) = (ATL_L1mvelts - (ATL_mvNMU<<1)) / ((ATL_mvNMU<<1)+1); \\\n");
         fprintf(fp, "   if (ATL_mvpagesize > (lda_)) \\\n   { \\\n");
         fprintf(fp,
            "      *(mb_) = (ATL_mvpagesize / (lda_)) * ATL_mvntlb; \\\n");
         fprintf(fp,
            "      if ( *(mb_) < *(nb_) ) *(nb_) = *(mb_); \\\n");
         fprintf(fp, "   } \\\n");
         fprintf(fp,
            "   else if (ATL_mvntlb < *(nb_)) *(nb_) = ATL_mvntlb; \\\n");
         fprintf(fp,
"   if (*(nb_) > ATL_mvNNU) *(nb_) = (*(nb_) / ATL_mvNNU) * ATL_mvNNU; \\\n");
         fprintf(fp, "   else *(nb_) = ATL_mvNNU; \\\n");
         fprintf(fp,
"   *(mb_) = (ATL_L1mvelts - *(nb_) * (ATL_mvNMU+1)) / (*(nb_)+2); \\\n");
         fprintf(fp,
"   if (*(mb_) > ATL_mvNMU) *(mb_) = (*(mb_) / ATL_mvNMU) * ATL_mvNMU; \\\n");
         fprintf(fp, "   else *(mb_) = ATL_mvNMU; \\\n");
      }
      fprintf(fp, "}\n");
   }
   else
   {
      if (mu != 0)
      {
         fprintf(fp, "#ifndef ATL_mvNNU\n");
         fprintf(fp, "   #include \"atlas_%cmvN.h\"\n", pre);
         fprintf(fp, "#endif\n");
      }
      fprintf(fp, "\n");
      fprintf(fp, "#define ATL_GetPartMVT(A_, lda_, mb_, nb_) \\\n{ \\\n");
      if (mu == 0)
      {
         fprintf(fp,
         "   *(mb_) = (ATL_L1mvelts - NB - ATL_mmMU*(NB+1)) / (NB+1); \\\n");
         fprintf(fp, "   if (*(mb_) > NB) \\\n   { \\\n");
         fprintf(fp, "      *(mb_) = ATL_MulByNB(ATL_DivByNB(*(mb_))); \\\n");
         fprintf(fp, "      *(nb_) = NB; \\\n   } \\\n");
         fprintf(fp, "   else \\\n   { \\\n");
         fprintf(fp,
"      if (*(mb_) > ATL_mmMU) *(nb_) = (*(mb_) / ATL_mmMU)*ATL_mmMU; \\\n");
         fprintf(fp, "      else *(nb_) = ATL_mmMU; \\\n");
         fprintf(fp, "      *(mb_) = NB; \\\n   } \\\n");
         fprintf(fp, "}\n");
      }
      else
      {
         fprintf(fp,
"   *(mb_) = (ATL_L1mvelts - (ATL_mvTMU<<1)) / ((ATL_mvTMU<<1)+1); \\\n");
         fprintf(fp,
"   if (*(mb_) > ATL_mvTNU) *(mb_) = (*(mb_)/ATL_mvTNU)*ATL_mvTNU; \\\n");
         fprintf(fp, "   else (*mb_) = ATL_mvTNU; \\\n");
         fprintf(fp, "   *(nb_) = ATL_mvTMU; \\\n");
         fprintf(fp, "}\n");
      }
   }

   fprintf(fp, "\n#endif\n");
   fclose(fp);
}

void PrintUsage(char *nam)
{
   fprintf(stderr,
      "USAGE: %s -p <pre> -l <l1mul> -N/T/S [-x <xu> -y <yu> -f <flag> ]\n",
           nam);
   exit(-1);
}

void GetFlags(int nargs, char **args, char *pre, int *l1mul, int *flagN,
              int *muN, int *nuN, int *flagT, int *muT, int *nuT,
              int *flagS, int *muS, int *nuS)
{
   int NoTrans=1, i;

   *pre = 'd';
   *l1mul = 80;
   *flagS = *flagN = *flagT = 0;
   *muT = *muN = 16;
   *nuN = *nuT = 4;
   *muS = *nuS = -1;

   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'p' :
         i++;
         *pre = Mlowcase(args[i][0]);
         break;
      case 'T':
         NoTrans = 0;
         break;
      case 'N':
         NoTrans = 1;
         break;
      case 'S':
         NoTrans = 2;
         break;
      case 'y':
         if (NoTrans == 1) *muN = atoi(args[++i]);
         else if (NoTrans == 2) *muS = atoi(args[++i]);
         else *muT = atoi(args[++i]);
         break;
      case 'x':
         if (NoTrans == 1) *nuN = atoi(args[++i]);
         else if (NoTrans == 2) *nuS = atoi(args[++i]);
         else *nuT = atoi(args[++i]);
         break;
      case 'f':
         if (NoTrans == 1) *flagN = atoi(args[++i]);
         else if (NoTrans == 2) *flagS = atoi(args[++i]);
         else *flagT = atoi(args[++i]);
         break;
      case 'l':
         *l1mul = atoi(args[++i]);
         break;
      default:
         fprintf(stderr, "Unknown flag = %s\n", args[i]);
         PrintUsage(args[0]);
      }
   }
   if (*muS == -1) *muS = *muT;
   if (*nuS == -1) *nuS = *nuT;
}

main(int nargs, char **args)
{
   int flagN, muN, nuN, flagT, muT, nuT, flagS, muS, nuS, l1mul;
   char pre;
   GetFlags(nargs, args, &pre, &l1mul, &flagN, &muN, &nuN, &flagT, &muT, &nuT,
            &flagS, &muS, &nuS);

   emit_mvhead(pre, l1mul);
   emit_transhead(pre, 'N', flagN, muN, nuN, l1mul);
   emit_transhead(pre, 'T', flagT, muT, nuT, l1mul);
   emit_transhead(pre, 'S', flagS, muS, nuS, l1mul);
   exit(0);
}
