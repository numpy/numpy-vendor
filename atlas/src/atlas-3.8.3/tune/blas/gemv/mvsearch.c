/*
 *             Automatically Tuned Linear Algebra Software v3.8.3
 *                    (C) Copyright 1999 R. Clint Whaley
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
#include <ctype.h>
#include <assert.h>
#include "atlas_fopen.h"

#define ATL_MVIsMM(iflag_) ( ((iflag_) | 8) == (iflag_) )
#define ATL_MVIsAxpy(iflag_) ( ((iflag_) | 16) == (iflag_) )
#define ATL_MVIsDot(iflag_) ( !ATL_MVIsMM(iflag_) && !ATL_MVIsAxpy(iflag_) )
#define ATL_MVNoBlock(iflag_) ( ((iflag_) | 32) == (iflag_) )

double GetAvg(int n, double tolerance, double *mflop)
{
   int i, j;
   double t0, tavg;
/*
 * Sort results, largest first
 */
   for (i=0; i != n; i++)
   {
      for (j=i+1; j < n; j++)
      {
         if (mflop[i] < mflop[j])
         {
            t0 = mflop[i];
            mflop[i] = mflop[j];
            mflop[j] = t0;
         }
      }
   }
/*
 * Throw out result if it is outside tolerance; rerun if two mflop not within
 * tolerance;  this code assumes n == 3
 */
   if (tolerance*mflop[1] < mflop[0])  /* too big a range in results */
   {
      if (tolerance*mflop[2] < mflop[1]) return(-1.0);
      tavg = (mflop[1] + mflop[2]) / 2.0;
   }
   else if (tolerance*mflop[2] < mflop[0]) tavg = (mflop[0] + mflop[1]) / 2.0;
   else tavg = (mflop[0] + mflop[1] + mflop[2]) / 3.0;

   return(tavg);
}

int GetL1CacheSize()
{
   FILE *L1f;
   int L1Size;

   if (!FileExists("res/L1CacheSize"))
   {
      assert(system("make res/L1CacheSize\n") == 0);
   }
   L1f = fopen("res/L1CacheSize", "r");
   assert(L1f != NULL);
   fscanf(L1f, "%d", &L1Size);
   fclose(L1f);
   fprintf(stderr, "\n      Read in L1 Cache size as = %dKB.\n",L1Size);
   return(L1Size);
}

int mvtstcase(char pre, char TA, char *mvnam, char *cc, char *ccflags)
{
   char ln[256];
   if (cc[0] != '\0' && ccflags[0] != '\0')
      sprintf(ln,
              "make %cmvtstcase%c mvrout=%s %cMVCC=\"%s\" %cMVFLAGS=\"%s\"\n",
              pre, TA, mvnam, pre, cc, pre, ccflags);
   else if (cc[0] == '\0' && ccflags[0] != '\0')
      sprintf(ln, "make %cmvtstcase%c mvrout=%s EXTFLAGS=\"%s\"\n",
              pre, TA, mvnam, ccflags);
   else sprintf(ln, "make %cmvtstcase%c mvrout=%s\n", pre, TA, mvnam);
   return(system(ln));
}

double svcase
   (char pre, int cas, int l1mul, int mu, int nu, int flag,
    char *mvnamN, char *ccN, char *ccflagsN,
    char *mvnamT, char *ccT, char *ccflagsT)
{
   char ln[1024], fnam[64];
   const int imul = l1mul;
   int i, mb;
   double mfs[3], mf;
   FILE *fp;

   sprintf(fnam, "res/%csymvS_%d_0", pre, cas);
   if (!FileExists(fnam))
   {
      i = sprintf(ln, "make %csvcase mvrout=%s svrout=%s cas=%d xu=%d yu=%d l1mul=%d iflag=%d gmvout=\"-o %s\"",
              pre, mvnamT, mvnamN, cas, nu, mu, imul, flag, fnam);
      if (ccN[0] != '\0' && ccflagsN[0] != '\0')
         i += sprintf(ln+i, " %cSVCC=\"%s\" %cSVFLAGS=\"%s\"",
                      pre, ccN, pre, ccflagsN);
      else if (ccN[0] == '\0' && ccflagsN[0] != '\0')
      {
         i += sprintf(ln+i,
           " SEXTFLAGS=\"%s -DATL_NOL1PREFETCH -DATL_NOL2PREFETCH\"", ccflagsN);
      }
      else i += sprintf(ln+i,
           " SEXTFLAGS=\"-DATL_NOL1PREFETCH -DATL_NOL2PREFETCH\"");
      if (ccT[0] != '\0' && ccflagsT[0] != '\0')
         i += sprintf(ln+i, " %cMVCC=\"%s\" %cMVFLAGS=\"%s\"",
                      pre, ccT, pre, ccflagsT);
      else if (ccT[0] == '\0' && ccflagsT[0] != '\0')
         i += sprintf(ln+i, " SEXTFLAGS=\"%s\"", ccflagsT);
      sprintf(ln+i, "\n");
      fprintf(stderr, "%s", ln);
      if (system(ln)) return(-1.0);  /* won't compile here */
   }
   fp = fopen(fnam, "r");
   assert(fp);
   assert(fscanf(fp, " %lf %lf %lf", mfs, mfs+1, mfs+2) == 3);
   fclose(fp);
   mf = GetAvg(3, 1.20, mfs);
   if (mf == -1.0)
   {
      fprintf(stderr,
"\n\n%s : VARIATION EXCEEDS TOLERENCE, RERUN WITH HIGHER REPS.\n\n", fnam);
      sprintf(ln, "rm -f %s\n", fnam);
      system(ln);
      exit(-1);
   }
   return(mf);
}

double mvcase(int SY, char pre, char *mvnam, char TA, int flag, int mu, int nu,
              int cas, char *cc, char *ccflags, int l1mul)
{
   char nTA;
   char ln[1024], fnam[64];
   const int imul = l1mul;
   int i, mb;
   double mfs[3], mf;
   FILE *fp;

   if (TA == 'n' || TA == 'N') nTA = 'T';
   else nTA = 'N';
   if (mu == 0) mu = 1;
   if (nu == 0) nu = 1;

   if (SY) sprintf(fnam, "res/%csymv%c_%d_0", pre, TA, cas);
   else
   {
      if (ATL_MVNoBlock(flag)) sprintf(fnam, "res/%cgemv%c_%d_0", pre, TA, cas);
      else sprintf(fnam, "res/%cgemv%c_%d_%d", pre, TA, cas, imul);
   }
   if (!FileExists(fnam))
   {
      i = sprintf(ln, "make %cmvcase%c mvrout=%s cas=%d xu=%d yu=%d l1mul=%d iflag=%d gmvout=\"-o %s\"",
              pre, TA, mvnam, cas, nu, mu, imul, flag, fnam);
      if (cc[0] != '\0' && ccflags[0] != '\0')
         i += sprintf(ln+i, " %cMVCC=\"%s\" %cMVFLAGS=\"%s\"",
                      pre, cc, pre, ccflags);
      else if (cc[0] == '\0' && ccflags[0] != '\0')
      {
         if (SY) i += sprintf(ln+i,
            " EXTFLAGS=\"%s -DATL_NOL1PREFETCH -DATL_NOL2PREFETCH\"", ccflags);
         else i += sprintf(ln+i, " EXTFLAGS=\"%s\"", ccflags);
      }
      if (SY)
      {
         if (pre == 'd' || pre == 'c') mb = 1024 / 8;
         else if (pre == 's') mb = 1024 / 4;
         else if (pre == 'z') mb = 1024 / 16;
         mb = 0.01*l1mul*mb*GetL1CacheSize();
         mb = (mb-2*nu)/(nu+2);
         if (mb > mu) mb = (mb/mu)*mu;
         if (!(cc[0] == '\0' && ccflags[0] != '\0'))
            i += sprintf(ln+i,
                         " EXTFLAGS=\"-DATL_NOL1PREFETCH -DATL_NOL2PREFETCH\"");
         i += sprintf(ln+i, " opt=\"-2 1 -L %d\" M=%d N=%d", mb, mb, nu);
      }
      sprintf(ln+i, "\n");
      fprintf(stderr, "%s", ln);
      if (system(ln)) return(-1.0);  /* won't compile here */
   }
   fp = fopen(fnam, "r");
   assert(fp);
   assert(fscanf(fp, " %lf %lf %lf", mfs, mfs+1, mfs+2) == 3);
   fclose(fp);
   mf = GetAvg(3, 1.20, mfs);
   if (mf == -1.0)
   {
      fprintf(stderr,
"\n\n%s : VARIATION EXCEEDS TOLERENCE, RERUN WITH HIGHER REPS.\n\n", fnam);
      sprintf(ln, "rm -f %s\n", fnam);
      system(ln);
      exit(-1);
   }
   return(mf);
}

int FindL1Mul(char pre, int cas, char *mvnam, char *cc, char *ccflags,
              char TA, int flag, int mu, int nu)
{
   double low = .5, high = 1.0;
   double mflow, mfhigh;
   int ilow, ihigh;

   if (ATL_MVNoBlock(flag)) flag -= 32;  /* always actually block these times */
   do
   {
      ilow = (low  * 100.0);
      ihigh = (high * 100.0);
      mflow  = mvcase(0, pre, mvnam, TA, flag, mu, nu, cas, cc, ccflags, ilow);
      mfhigh = mvcase(0, pre, mvnam, TA, flag, mu, nu, cas, cc, ccflags, ihigh);
      fprintf(stdout, "      %d%% %.2fMFLOP  ---  %d%% %.2fMFLOP\n",
              ilow, mflow, ihigh, mfhigh);
      if (mflow < 1.005*mfhigh) low += 0.5*(high-low);
      else high -= 0.5 * (high-low);
   }
   while (ihigh-ilow);
   fprintf(stdout, "\n\nBEST %% of L1 cache: %d\n", ilow);
   return(ilow);
}

int ConfirmBlock(char pre, char *mvnam, char *cc, char *ccflags, char TA,
                 int flag, int mu, int nu, int cas, int l1mul)
{
   int bflag;
   double mfblock, mfnoblock;

   if ( ATL_MVNoBlock(flag) )
   {
      bflag = flag - 32;
      mfblock   = mvcase(0, pre, mvnam, TA, bflag, mu, nu, cas, cc, ccflags,
                         l1mul);
      mfnoblock = mvcase(0, pre, mvnam, TA,  flag, mu, nu, cas, cc, ccflags,
                         l1mul);
      fprintf(stdout, "\nWith blocking=%lf, without=%lf\n\n",
              mfblock, mfnoblock);
      if (mfblock >= mfnoblock) return(bflag);
   }
   return(flag);
}

void NoEndLineWhiteSpace(char *ln)
{
   int i;

   for (i=0; ln[i]; i++);
   if (i)
      for (i--; isspace(ln[i]); i--) ln[i] = '\0';
}

int LineIsCont(char *ln)
{
   int i, iret=0;
   for(i=0; ln[i]; i++);
   if (i)
   {
      for(i--; isspace(ln[i]); i--);
      if (ln[i] == '\\') iret = 1;
   }
   return(iret);
}

int iGetIDCase(FILE *fp, int ID, char *fnam, char *auth, char *cc,
               char *ccflags, int *flag, int *mu, int *nu)
{
   char ln[256];
   int i, n, id;

   assert(fgets(ln, 128, fp));
   assert(sscanf(ln, " %d", &n) == 1);
   for (i=0; i < n; i++)
   {
      assert(fgets(ln, 256, fp));
      assert(sscanf(ln, " %d %d %d %d %s \"%[^\"]",
                    &id, flag, mu, nu, fnam, auth) == 6);
      if (LineIsCont(ln))
      {
         assert(fgets(cc, 256, fp));
         assert(fgets(ccflags, 256, fp));
         NoEndLineWhiteSpace(cc);
         NoEndLineWhiteSpace(ccflags);
      }
      else cc[0] = ccflags[0] = '\0';
      if (id == ID) break;
   }
   return(id);
}

void GetIDCase(char pre, int ID, char *fnam, char *auth, char *cc,
               char *ccflags, int *flag, int *mu, int *nu)
{
   char ln[256];
   int i, n, id;
   FILE *fp;

   sprintf(ln, "%ccases.dsc", pre);
   fp = fopen(ln, "r");
   assert(fp);

   id = iGetIDCase(fp, ID, fnam, auth, cc, ccflags, flag, mu, nu);
   if (id != ID) id = iGetIDCase(fp, ID, fnam, auth, cc, ccflags, flag, mu, nu);

   fclose(fp);
   assert(id == ID);
}

void GetCases(FILE *fp, int *N, int **IDs, char ***fnams, char ***auths,
              char ***CC, char ***CCFLAGS, int **flags, int **mus, int **nus)
{
   int i, j, n;
   int *mu, *nu, *flag, *id;
   char **fnam, **auth, **cc, **ccflags;
   char ln[256];

   assert(fgets(ln, 128, fp));
   assert(sscanf(ln, " %d", &n) == 1);
   assert(n < 100 && n > 0);
   fnam    = malloc(n * sizeof(char*));
   auth    = malloc(n * sizeof(char*));
   cc      = malloc(n * sizeof(char*));
   ccflags = malloc(n * sizeof(char*));
   assert(fnam && auth && cc && ccflags);
   for (i=0; i < n; i++)
   {
      assert(fnam[i] = malloc(64*sizeof(char)));
      assert(auth[i] = malloc(64*sizeof(char)));
   }
   id = malloc(n * sizeof(int));
   mu = malloc(n * sizeof(int));
   nu = malloc(n * sizeof(int));
   flag = malloc(n * sizeof(int));
   assert(id && mu && nu && flag);
   for (i=0; i < n; i++)
   {
      assert(fgets(ln, 256, fp));
      assert(sscanf(ln, " %d %d %d %d %s \"%[^\"]",
                    id+i, flag+i, mu+i, nu+i, fnam[i], auth[i]) == 6);
      assert(mu[i] >= 0 && nu[i] >= 0 && fnam[i][0] != '\0' && id[i] > 0);
      if (LineIsCont(ln))
      {
         fgets(ln, 256, fp);
         NoEndLineWhiteSpace(ln);
         if (ln[0] == '+' && (ln[1] == '\0' || ln[1] == ' ')) /* adding flags */
            cc[i] = NULL;
         else
         {
            j = strlen(ln);
            cc[i] = malloc((j+1)*sizeof(char));
            assert(cc[i]);
            strcpy(cc[i], ln);
         }
         fgets(ln, 256, fp);
         NoEndLineWhiteSpace(ln);
         j = strlen(ln);
         ccflags[i] = malloc((j+1)*sizeof(char));
         assert(ccflags[i]);
         strcpy(ccflags[i], ln);
      }
      else cc[i] = ccflags[i] = NULL;
   }
   for (i=0; i < n; i++) /* fix cc & ccflags */
   {
      if (cc[i] == NULL)
      {
         cc[i] = malloc(2*sizeof(char));
         assert(cc[i]);
         cc[i][0] = cc[i][1] = '\0';
      }
      if (ccflags[i] == NULL)
      {
         ccflags[i] = malloc(2*sizeof(char));
         assert(ccflags[i]);
         ccflags[i][0] = ccflags[i][1] = '\0';
      }
   }

   *N = n;
   *IDs = id;
   *fnams = fnam;
   *auths = auth;
   *CC = cc;
   *CCFLAGS = ccflags;
   *flags = flag;
   *mus = mu;
   *nus = nu;
}

int RunTransCases(int SY, char pre, char TA, int ncases, int *ids, char **fnams,
                  char **cc, char **ccflags, int *flags, int *mus, int *nus)
{
   int i, imax=0;
   double mf, mfmax=0.0;

   for (i=0; i < ncases; i++)
   {
      mf = mvcase(SY, pre, fnams[i], TA, flags[i], mus[i], nus[i], ids[i],
                  cc[i], ccflags[i], 75);
      fprintf(stdout, "%s : %.2f\n", fnams[i], mf);
      if (mf > mfmax)
      {
         if (mvtstcase(pre, TA, fnams[i], cc[i], ccflags[i]) == 0)
         { /* ensure it passes test */
            mfmax = mf;
            imax = i+1;
         }
         else fprintf(stderr, "\n\nROUTINE %s FAILED TESTS!!!\n\n", fnams[i]);
      }
   }
   assert(imax);
   fprintf(stdout,
           "\nbest %cgemv%c : ID=%d, mu=%d, nu=%d at %.2f MFLOPS\n\n",
           pre, TA, ids[imax-1], mus[imax-1], nus[imax-1], mfmax);
   return(imax-1);
}

void ReadSum(char pre, int *l1mul,
int *idN, char *fnamN, char *authN, int *flagN, int *muN, int *nuN, double *mfN,
int *idT, char *fnamT, char *authT, int *flagT, int *muT, int *nuT, double *mfT,
int *idS, char *fnamS, char *authS, int *flagS, int *muS, int *nuS, double *mfS)
{
   char ln[256];
   int i;
   FILE *fp;
   double t0;

   sprintf(ln, "res/%cMVRES", pre);
   fp = fopen(ln, "r");
   assert(fp);

   assert(fgets(ln, 256, fp));
   assert(sscanf(ln, " %d %d %d %d %d %lf %s \"%[^\"]",
                 idN, flagN, muN, nuN, l1mul, mfN, fnamN, authN) == 8);
   assert(fgets(ln, 256, fp));
   assert(sscanf(ln, " %d %d %d %d %d %lf %s \"%[^\"]",
                 idT, flagT, muT, nuT, l1mul, mfT, fnamT, authT) == 8);
   assert(fgets(ln, 256, fp));
   assert(sscanf(ln, " %d %d %d %d %d %lf %s \"%[^\"]",
                 idS, flagS, muS, nuS, &i, mfS, fnamS, authS) == 8);
   fclose(fp);
}

void CreateSum
   (char pre, int l1mul,
    int idN, char *fnamN, char *authN, int flagN, int muN, int nuN, double mfN,
    int idT, char *fnamT, char *authT, int flagT, int muT, int nuT, double mfT,
    int idS, char *fnamS, char *authS, int flagS, int muS, int nuS, double mfS)
{
   char fnam[32];
   FILE *fp;

   sprintf(fnam, "res/%cMVRES", pre);
   fp = fopen(fnam, "w");
   assert(fp);
   fprintf(fp, "%d %d %d %d %d %.2f %s \"%s\"\n",
           idN, flagN, muN, nuN, l1mul, mfN, fnamN, authN);
   fprintf(fp, "%d %d %d %d %d %.2f %s \"%s\"\n",
           idT, flagT, muT, nuT, l1mul, mfT, fnamT, authT);
   fprintf(fp, "%d %d %d %d %d %.2f %s \"%s\"\n",
           idS, flagS, muS, nuS, l1mul, mfS, fnamS, authS);
   fclose(fp);
}

void GenMakefile(char pre, char *ccN, char *ccflagsN, char *ccT, char *ccflagsT,
                 char *ccS, char *ccflagsS)
{
   char ln[256];
   char renam[128];
   char *bet[4] = {"0", "1", "X", "XI0"};
   char *tts[6] = {"N", "T", "S", "N",  "T"};
   char *TAs[6] = {"N", "T", "S", "Nc", "C"};
   char *typ, *ct;
   char UPR;
   FILE *fp;
   int i, j, nbet=3, ntran=3;

   if (pre == 's') typ = "SREAL";
   else if (pre == 'd') typ = "DREAL";
   else if (pre == 'c') typ = "SCPLX";
   else if (pre == 'z') typ = "DCPLX";
   else assert(0);
   UPR = (pre == 's' || pre == 'c') ? 'S' : 'D';
   sprintf(ln, "Make_%cmv", pre);
   fp = fopen(ln, "w");
   assert(fp);
   fprintf(fp, "include Make.inc\n\n");

   if (ccN[0] == '\0') fprintf(fp, "MVNCC = $(%cKC)\n", UPR);
   else fprintf(fp, "MVNCC = %s\n", ccN);
   if (ccN[0] == '\0' && ccflagsN[0] != '\0')
      fprintf(fp, "MVNFLAGS = $(CDEFS) $(%cKCFLAGS) %s -D%s\n",
              UPR, ccflagsN, typ);
   else if (ccflagsN[0] != '\0')
      fprintf(fp, "MVNFLAGS = $(CDEFS) %s -D%s\n", ccflagsN, typ);
   else fprintf(fp, "MVNFLAGS = $(CDEFS) $(%cKCFLAGS) -D%s\n", UPR, typ);
   if (pre == 'c' || pre == 'z')
   {
      fprintf(fp, "MVNcCC = $(MVNCC)\n");
      fprintf(fp, "MVNcFLAGS = $(MVNFLAGS) -DConj_\n\n");
   }

   if (ccT[0] == '\0') fprintf(fp, "MVTCC = $(%cKC)\n", UPR);
   else fprintf(fp, "MVTCC = %s\n", ccT);
   if (ccT[0] == '\0' && ccflagsT[0] != '\0')
      fprintf(fp, "MVTFLAGS = $(CDEFS) $(%cKCFLAGS) %s -D%s\n",
              UPR, ccflagsT, typ);
   else if (ccflagsT[0] != '\0')
      fprintf(fp, "MVTFLAGS = $(CDEFS) %s -D%s\n", ccflagsT, typ);
   else fprintf(fp, "MVTFLAGS = $(CDEFS) $(%cKCFLAGS) -D%s\n", UPR, typ);
   if (pre == 'c' || pre == 'z')
   {
      fprintf(fp, "MVCCC = $(MVTCC)\n");
      fprintf(fp, "MVCFLAGS = $(MVTFLAGS) -DConj_\n\n");
   }

   if (ccS[0] == '\0') fprintf(fp, "MVSCC = $(%cKC)\n", UPR);
   else fprintf(fp, "MVSCC = %s\n", ccS);
   if (ccS[0] == '\0' && ccflagsS[0] != '\0')
      fprintf(fp, "MVSFLAGS = $(CDEFS) $(%cKCFLAGS) %s -D%s",
              UPR, ccflagsS, typ);
   else if (ccflagsS[0] != '\0')
      fprintf(fp, "MVSFLAGS = $(CDEFS) %s -D%s", ccflagsS, typ);
   else fprintf(fp, "MVSFLAGS = $(CDEFS) $(%cKCFLAGS) -D%s", UPR, typ);
   fprintf(fp, " -DATL_NOL1PREFETCH -DATL_NOL2PREFETCH\n\n");

   if (pre == 's' || pre == 'd')
   {
      fprintf(fp, "\nobj = ");
      for (i=0; i < ntran; i++)
      {
         fprintf(fp,
         "\\\n      ATL_%cgemv%s_b0.o ATL_%cgemv%s_b1.o ATL_%cgemv%s_bX.o ",
                 pre, TAs[i], pre, TAs[i], pre, TAs[i]);
      }
   }
   else
   {
      nbet = 4;
      ntran = 5;
      fprintf(fp, "\nobj = ");
      for (i=0; i < ntran; i++)
      {
         fprintf(fp,
"\\\n      ATL_%cgemv%s_b0.o ATL_%cgemv%s_b1.o ATL_%cgemv%s_bX.o ATL_%cgemv%s_bXI0.o ",
                 pre, TAs[i], pre, TAs[i], pre, TAs[i], pre, TAs[i]);
      }
   }
   fprintf(fp, "\n\n");
   fprintf(fp, "%clib : lib\n", pre);
   fprintf(fp, "lib : %cmvk.grd\n", pre);
   fprintf(fp, "%cmvk.grd : $(obj)\n", pre);
   fprintf(fp, "\t$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(obj)\n");
   fprintf(fp, "\t$(RANLIB) $(ATLASlib)\n");
   fprintf(fp, "\ttouch %cmvk.grd\n\n", pre);

   fprintf(fp, "clean :\n\t rm -f $(obj) %cmvk.grd\n", pre);
   fprintf(fp, "%cclean : clean\n\n", pre);

   for (j=0; j < ntran; j++)
   {
      for (i=0; i < nbet; i++)
      {
         if (TAs[j][0] == 'S')
         {
            if (i == 3)  /* hack to get cap correct without rewriting */
               sprintf(renam,
                  " -DATL_%cgemvN_a1_x1_bXi0_y1=ATL_%cgemvS_a1_x1_bXi0_y1 ",
                       pre, pre);
            else
               sprintf(renam,
                       " -DATL_%cgemvN_a1_x1_b%s_y1=ATL_%cgemvS_a1_x1_b%s_y1 ",
                       pre, bet[i], pre, bet[i]);
            fprintf(fp, "ATL_%cgemvS_b%s.o : ATL_%cgemvS.c\n",
                    pre, bet[i], pre);
            fprintf(fp,
                    "\t$(MVSCC) -o $@ -c $(MVSFLAGS) -DBETA%s -DALPHA1 %s \\\n",
                    bet[i], renam);
            fprintf(fp, "                 ATL_%cgemvS.c\n", pre);
         }
         else
         {
            fprintf(fp, "ATL_%cgemv%s_b%s.o : ATL_%cgemv%s.c\n",
                    pre, TAs[j], bet[i], pre, tts[j]);
            fprintf(fp,
"\t$(MV%sCC) -o $@ -c $(MV%sFLAGS) -DBETA%s -DALPHA1 ATL_%cgemv%s.c\n",
                    TAs[j], TAs[j], bet[i], pre, tts[j]);
         }
      }
   }
   fclose(fp);
}

void mvinstall(char pre, int l1mul,
               int idN, char *fnamN, char *authN, char *ccN, char *ccflagsN,
               int flagN, int muN, int nuN,
               int idT, char *fnamT, char *authT, char *ccT, char *ccflagsT,
               int flagT, int muT, int nuT,
               int idS, char *fnamS, char *authS, char *ccS, char *ccflagsS,
               int flagS, int muS, int nuS)
{
   char ln[256];
   double mfN, mfT, mfS;

   GenMakefile(pre, ccN, ccflagsN, ccT, ccflagsT, ccS, ccflagsS);
   mfS = mvcase(1,pre, fnamS, 'N', flagS, muS, nuS, 5000, ccT, ccflagsT, l1mul);
   mfN = mvcase(0,pre, fnamN, 'N', flagN, muN, nuN, 5000, ccN, ccflagsN, l1mul);
   mfT = mvcase(0,pre, fnamT, 'T', flagT, muT, nuT, 5000, ccT, ccflagsT, l1mul);
   sprintf(ln, "make xemit_head \n");
   assert(system(ln) == 0);
   sprintf(ln, "./xemit_head -p %c -l %d -N -f %d -y %d -x %d -T -f %d -y %d -x %d -S -f %d -y %d -x %d\n", pre, l1mul, flagN, muN, nuN, flagT, muT, nuT, flagS, muS, nuS);
   fprintf(stderr, "%s", ln);
   assert(system(ln) == 0);
   sprintf(ln, "make %cinstall mvNrout=%s mvTrout=%s mvSrout=%s\n",
           pre, fnamN, fnamT, fnamS);
   fprintf(stderr, "%s", ln);
   assert(system(ln) == 0);
   CreateSum(pre, l1mul, idN, fnamN, authN, flagN, muN, nuN, mfN,
             idT, fnamT, authT, flagT, muT, nuT, mfT,
             idS, fnamS, authS, flagS, muS, nuS, mfS);
}

static int ilcm(const int m, const int n)
{
   int ia, iq, ir, ilcm;

   if (m == n) return(n);
   if (m >= n)
   {
      ia = m;
      ilcm = n;
   }
   else
   {
      ia = n;
      ilcm = m;
   }
   do
   {
      iq = ia / ilcm;
      ir = ia - iq * ilcm;
      if (ir == 0) return( (m*n)/ilcm );
      ia = ilcm;
      ilcm = ir;
   }
   while(1);
}

int RunSymmCases(char pre, char *fnamT, char *ccT, char *ccflagsT,
                 int muT, int nuT, int flagT, int ncases,
                 int *ids, char **fnams, char **cc, char **ccflags,
                 int *flags, int *mus, int *nus, int *muS, int *nuS)
{
   int i, mu, nu;
   int imax=(-1);
   double mf, mfB=0.0;

   if (ATL_MVIsMM(flagT)) /* if using gemm for trans, use it for notrans too */
   {
      *muS = *nuS = 80;
      for (i=0; i < ncases; i++) if (ATL_MVIsMM(flags[i])) return(i);
      assert(0);
   }
   for (i=0; i < ncases; i++)
   {
      if (ATL_MVIsAxpy(flags[i]))  /* only try axpy-based routs */
      {
         mu = ilcm(muT, nus[i]);
         nu = ilcm(nuT, mus[i]);
         mf = svcase(pre, ids[i], 75, mu, nu, flags[i], fnams[i],
                     cc[i], ccflags[i], fnamT, ccT, ccflagsT);
         if (mf > mfB)
         {
            mfB = mf;
            imax = i;
            *muS = mu;
            *nuS = nu;
         }
      }
   }
   assert(imax != -1);
   return(imax);
}

void RunCases(char pre)
{
   char fnam[128];
   char Nfnam[256], Tfnam[256], Nauth[256], Ncc[256], Tauth[256], Tcc[256];
   char Nccflags[512], Tccflags[512];
   char Sfnam[256], Sauth[256], Scc[256], Sccflags[512];
   char **fnamN, **fnamT, **authN, **authT, **ccN, **ccflagsN,
        **ccT, **ccflagsT;
   int i, nNTcases, nTcases, Nbest, Tbest, Sbest;
   int Nid, Nflag, Nmu, Nnu, Tid, Tflag, Tmu, Tnu;
   int Sid, Sflag, Smu, Snu;
   int *flagN, *muN, *nuN, *flagT, *muT, *nuT, *idN, *idT;
   int l1mul;
   FILE *fp;

/*
 * Read in cases to try
 */
   sprintf(fnam, "%ccases.dsc", pre);
   fp = fopen(fnam, "r");
   assert(fp);
   GetCases(fp, &nNTcases, &idN, &fnamN, &authN, &ccN, &ccflagsN, &flagN,
            &muN, &nuN);
   GetCases(fp, &nTcases,  &idT, &fnamT, &authT, &ccT, &ccflagsT, &flagT,
            &muT, &nuT);
   fclose(fp);
/*
 * Try all cases for each trans case
 */
   Nbest = RunTransCases(0, pre, 'N', nNTcases, idN, fnamN, ccN, ccflagsN,
                         flagN, muN, nuN);
   Tbest = RunTransCases(0, pre, 'T', nTcases, idT, fnamT, ccT, ccflagsT,
                         flagT, muT, nuT);
   Sbest = RunSymmCases(pre, fnamT[Tbest], ccT[Tbest], ccflagsT[Tbest],
                        muT[Tbest], nuT[Tbest], flagT[Tbest],
                        nNTcases, idN, fnamN, ccN, ccflagsN, flagN, muN, nuN,
                        &Smu, &Snu);
   Nid = idN[Nbest]; Nflag = flagN[Nbest]; Tflag = flagT[Tbest];
   Nmu = muN[Nbest]; Nnu = nuN[Nbest];
   strcpy(Nfnam, fnamN[Nbest]); strcpy(Nauth, authN[Nbest]);
   strcpy(Ncc, ccN[Nbest]);
   strcpy(Nccflags, ccflagsN[Nbest]);
   Tid = idT[Tbest]; Tmu = muT[Tbest]; Tnu = nuT[Tbest];
   strcpy(Tfnam, fnamT[Tbest]); strcpy(Tauth, authT[Tbest]);
   strcpy(Tcc, ccT[Tbest]);
   strcpy(Tccflags, ccflagsT[Tbest]);
   if (Sbest == -1) Sbest = Nbest;
   Sid = idN[Sbest]; Sflag = flagN[Sbest];
   strcpy(Sfnam, fnamN[Sbest]); strcpy(Sauth, authN[Sbest]);
   strcpy(Scc, ccN[Sbest]);
   strcpy(Sccflags, ccflagsN[Sbest]);

   free(flagN); free(flagT); free(muN); free(muT); free(nuN); free(nuT);
   for (i=0; i < nNTcases; i++)
   {
      free(fnamN[i]);
      free(authN[i]);
      free(ccN[i]);
      free(ccflagsN[i]);
   }
   free(fnamN); free(authN); free(ccN); free(ccflagsN);
   for (i=0; i < nTcases; i++)
   {
      free(fnamT[i]);
      free(authT[i]);
      free(ccT[i]);
      free(ccflagsT[i]);
   }
   free(fnamT); free(authT); free(ccT); free(ccflagsT);
   l1mul = FindL1Mul(pre, Tid, Tfnam, Tcc, Tccflags, 'T', Tflag, Tmu, Tnu);
   Tflag = ConfirmBlock(pre, Tfnam, Tcc, Tccflags, 'T', Tflag, Tmu, Tnu,
                        Tid, l1mul);
   Nflag = ConfirmBlock(pre, Nfnam, Ncc, Nccflags, 'N', Nflag, Nmu, Nnu,
                        Nid, l1mul);
   CreateSum(pre, l1mul, Nid, Nfnam, Nauth, Nflag, Nmu, Nnu, -1.0,
             Tid, Tfnam, Tauth, Tflag, Tmu, Tnu, -1.0,
             Sid, Sfnam, Sauth, Sflag, Smu, Snu, -1.0);
}

void GoToTown(pre)
{
   char fnamN[256], authN[256], fnamT[256], authT[256], ln[128];
   char ccN[256], ccflagsN[512], ccT[256], ccflagsT[512];
   char fnamS[256], authS[256], ccS[256], ccflagsS[512];
   int idN, flagN, muN, nuN, idT, flagT, muT, nuT, l1mul;
   int idS, flagS, muS, nuS, i, j;
   double mfN, mfT, mfS;

   sprintf(ln, "res/%cMVRES", pre);
   if (!FileExists(ln)) RunCases(pre);
   else  /* if default does not pass tester, rerun cases */
   {
      ReadSum(pre, &l1mul, &idN, fnamN, authN, &flagN, &muN, &nuN, &mfN,
              &idT, fnamT, authT, &flagT, &muT, &nuT, &mfT,
              &idS, fnamS, authS, &flagS, &muS, &nuS, &mfS);
      GetIDCase(pre, idN, fnamN, authN, ccN, ccflagsN, &flagN, &muN, &nuN);
      GetIDCase(pre, idT, fnamT, authT, ccT, ccflagsT, &flagT, &muT, &nuT);
      GetIDCase(pre, idS, fnamS, authS, ccS, ccflagsS, &flagS, &i, &j);
      if ( (mvtstcase(pre, 'N', fnamN, ccN, ccflagsN) != 0) ||
           (mvtstcase(pre, 'T', fnamT, ccT, ccflagsT) != 0) ||
           (mvtstcase(pre, 'N', fnamS, ccS, ccflagsS) != 0) )
      {
         assert(remove(ln) == 0);
         RunCases(pre);
      }
   }
   ReadSum(pre, &l1mul, &idN, fnamN, authN, &flagN, &muN, &nuN, &mfN,
           &idT, fnamT, authT, &flagT, &muT, &nuT, &mfT,
           &idS, fnamS, authS, &flagS, &muS, &nuS, &mfS);
   GetIDCase(pre, idN, fnamN, authN, ccN, ccflagsN, &flagN, &muN, &nuN);
   GetIDCase(pre, idT, fnamT, authT, ccT, ccflagsT, &flagT, &muT, &nuT);
   GetIDCase(pre, idS, fnamS, authS, ccS, ccflagsS, &flagS, &i, &j);
   mvinstall(pre, l1mul, idN, fnamN, authN, ccN, ccflagsN, flagN, muN, nuN,
             idT, fnamT, authT, ccT, ccflagsT, flagT, muT, nuT,
             idS, fnamS, authS, ccS, ccflagsS, flagS, muS, nuS);
}

void PrintUsage(char *fnam)
{
   fprintf(stderr, "USAGE: %s [-p <s,d,c,z>]\n", fnam);
   exit(-1);
}

void GetFlags(int nargs, char **args, char *pre, char *TA, int *mu, int *nu)
{
   char ctmp;
   int i;
   *pre = 'd';
   *mu = *nu = 0;
   *TA = ' ';
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(&args[0][0]);
      switch(args[i][1])
      {
      case 'p':
         ctmp = args[++i][0];
         ctmp = tolower(ctmp);
         if (ctmp == 's' || ctmp == 'd' || ctmp == 'c' || ctmp == 'z')
            *pre = ctmp;
         else PrintUsage(&args[0][0]);
         break;
      case 'm':
         *mu = atoi(args[++i]);
         assert(*mu > 0);
         break;
      case 'n':
         *nu = atoi(args[++i]);
         assert(*nu > 0);
         break;
      case 'A':
         ctmp = args[++i][0];
         ctmp = toupper(ctmp);
         if (ctmp == 'N' || ctmp == 'T') *TA = ctmp;
         else PrintUsage(args[0]);
         break;
      default:
         fprintf(stderr, "Unknown flag : %s\n", args[i]);
         PrintUsage(&args[0][0]);
      }
   }
   if (*mu == 0 && *nu) *mu = *nu;
   else if (*nu == 0 && *mu) *nu = *mu;
   if ( (*mu || *nu) && *TA == ' ') *TA = 'N';
}


main(int nargs, char **args)
{
   char pre, TA;
   int mu, nu;
   GetFlags(nargs, args, &pre, &TA, &mu, &nu);
   GoToTown(pre);
   exit(0);
}
