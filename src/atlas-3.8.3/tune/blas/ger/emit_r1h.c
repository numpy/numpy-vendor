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
#define ATL_NoBlock(iflag_) ( ((iflag_) | 32) == (iflag_) )
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

void emit_r1head(char pre, double l1mul, int flag, int mu, int nu)
{
   char fnam[64];
   int l1;
   FILE *fp;

   l1 = 1024 * GetL1CacheSize();
   if (pre == 's') l1 /= ATL_ssize;
   else if (pre == 'd') l1 /= ATL_dsize;
   else if (pre == 'c') l1 /= ATL_csize;
   else if (pre == 'z') l1 /= ATL_zsize;
   l1 = l1mul * l1;

   sprintf(fnam, "atlas_%cr1.h", pre);
   fp = fopen(fnam, "w");
   fprintf(fp, "#ifndef ATLAS_%cR1_H\n", toupper(pre));
   fprintf(fp, "#define ATLAS_%cR1_H\n\n", toupper(pre));

   if ( ATL_NoBlock(flag) ) fprintf(fp, "#define ATL_R1NOBLOCK\n");
   fprintf(fp, "#define ATL_L1r1elts %d\n", l1);
   fprintf(fp, "#define ATL_r1MU %d\n", mu);
   fprintf(fp, "#define ATL_r1NU %d\n", nu);
   fprintf(fp, "\n");

   fprintf(fp, "#define ATL_GetPartR1(A_, lda_, mb_, nb_) \\\n{ \\\n");
   fprintf(fp, "   (mb_) = (ATL_L1r1elts - (ATL_r1NU+ATL_r1NU)) / (ATL_r1NU+ATL_r1NU+1); \\\n");
   fprintf(fp,
      "   if ((mb_) > ATL_r1MU) (mb_) = ATL_r1MU*((mb_)/ATL_r1MU); \\\n");
   fprintf(fp, "   else (mb_) = ATL_r1MU; \\\n");
   fprintf(fp, "   (nb_) = ATL_r1NU; \\\n");
   fprintf(fp, "}\n");

   fprintf(fp, "\n#endif\n");
   fclose(fp);
}

void PrintUsage(char *nam)
{
   fprintf(stderr, "USAGE: %s -p <pre> -l <l1mul> -x <xu> -y <yu> -f <flag>\n",
           nam);
   exit(-1);
}

void GetFlags(int nargs, char **args, char *pre, double *l1mul,
              int *xu, int *yu, int *iflag)
{
   int i, j;

   *pre = 'd';
   *l1mul = 0.75;
   *xu = 16;
   *yu = 2;
   *iflag = 0;

   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'p' :
         i++;
         *pre = Mlowcase(args[i][0]);
         break;
      case 'l':
         *l1mul = atoi(args[++i])/100.0;
         break;
      case 'x':
         *xu = atoi(args[++i]);
         break;
      case 'y':
         *yu = atoi(args[++i]);
         break;
      case 'f':
         *iflag = atoi(args[++i]);
         break;
      }
   }
}

main(int nargs, char **args)
{
   double l1mul;
   int xu, yu, iflag;
   char pre;
   GetFlags(nargs, args, &pre, &l1mul, &xu, &yu, &iflag);
   emit_r1head(pre, l1mul, iflag, xu, yu);
   exit(0);
}

