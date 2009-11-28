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

#define ATL_NoBlock(iflag_) ( ((iflag_) | 32) == (iflag_) )

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


void emit_r1head(char pre, int l1mul, int flag, int mu, int nu)
{
   char ln[256];
   sprintf(ln, "./xemit_r1h -p %c -l %d -f %d -x %d -y %d\n",
           pre, l1mul, flag, mu, nu);
   assert(system(ln) == 0);
}

int r1tstcase(char pre, char *r1nam, char *cc, char *ccflags)
{
   char ln[256];

   if (cc[0] != '\0' && ccflags[0] != '\0')
      sprintf(ln,
              "make %cr1tstcase r1rout=%s %cR1CC=\"%s\" %cR1FLAGS=\"%s\"\n",
              pre, r1nam, pre, cc, pre, ccflags);
   else if (cc[0] == '\0' && ccflags[0] != '\0')
      sprintf(ln, "make %cr1tstcase r1rout=%s EXTFLAGS=\"%s\"\n",
              pre, r1nam, ccflags);
   else sprintf(ln, "make %cr1tstcase r1rout=%s\n", pre, r1nam);


   return(system(ln));
}

double r1case(char pre, char *r1nam, int flag, int mu, int nu, int cas,
              char *cc, char *ccflags, int l1mul)
{
   char fnam[128], ln[128];
   double mf, mfs[3];
   int i;
   FILE *fp;

   if (ATL_NoBlock(flag)) sprintf(fnam, "res/%cger1_%d_0", pre, cas);
   else sprintf(fnam, "res/%cger1_%d_%d", pre, cas, l1mul);
   if (!FileExists(fnam))
   {
      i = sprintf(ln,
              "make %cr1case r1rout=%s cas=%d l1mul=%d iflag=%d xu=%d yu=%d",
              pre, r1nam, cas, l1mul, flag, mu, nu);
      if (cc[0] != '\0' && ccflags[0] != '\0')
         i += sprintf(ln+i, " %cR1CC=\"%s\" %cR1FLAGS=\"%s\"",
                      pre, cc, pre, ccflags);
      else if (cc[0] == '\0' && ccflags[0] != '\0')
         i += sprintf(ln+i, " EXTFLAGS=\"%s\"", ccflags);
      i += sprintf(ln+i, "\n");
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

void CreateSum(char pre, int l1mul, int ID, char *r1nam, char *auth,
               int flag, int mu, int nu, double mf)
{
   char fnam[64];
   FILE *fp;

   sprintf(fnam, "res/%cR1RES", pre);
   fp = fopen(fnam, "w");
   assert(fp);
   fprintf(fp, "%d %d %d %d %d %.2f %s \"%s\"\n", ID, flag, mu, nu, l1mul, mf,
           r1nam, auth);
   fclose(fp);
}

void ReadSum(char pre, int *l1mul, int *ID, char *r1nam, char *auth,
             int *flag, int *mu, int *nu, double *mf)
{
   char fnam[64];
   FILE *fp;

   sprintf(fnam, "res/%cR1RES", pre);
   fp = fopen(fnam, "r");
   assert(fp);
   assert(fscanf(fp, " %d %d %d %d %d %lf %s \"%[^\"]",
                 ID, flag, mu, nu, l1mul, mf, r1nam, auth) == 8);
   fclose(fp);
}

void GenMakefile(char pre, char *cc, char *ccflags)
{
   int i;
   FILE *fp;
   char *typ, *comp, *cflags;
   char ln[64];

   if (pre == 's') typ = "SREAL";
   else if (pre == 'd') typ = "DREAL";
   else if (pre == 'c') typ = "SCPLX";
   else if (pre == 'z') typ = "DCPLX";
   else assert(0);
   if (pre == 's' || pre == 'c')
   {
      comp = "$(SKC)";
      cflags = "$(SKCFLAGS)";
   }
   else
   {
      comp = "$(DKC)";
      cflags = "$(DKCFLAGS)";
   }
   sprintf(ln, "Make_%cr1", pre);
   fp = fopen(ln, "w");
   assert(fp);
   fprintf(fp, "include Make.inc\n\n");

   if (cc[0] == '\0') fprintf(fp, "R1CC = %s\n", comp);
   else fprintf(fp, "R1CC = %s\n", cc);
   if (cc[0] == '\0' && ccflags[0] != '\0')
      fprintf(fp, "R1FLAGS = $(CDEFS) %s %s -D%s\n", cflags, ccflags, typ);
   else if (ccflags[0] != '\0')
      fprintf(fp, "R1FLAGS = $(CDEFS) %s -D%s\n", ccflags, typ);
   else fprintf(fp, "R1FLAGS = $(CDEFS) %s -D%s\n", cflags, typ);

   if (pre == 's' || pre == 'd')
      fprintf(fp, "\nobj = ATL_%cger1_a1_x1_yX.o\n\n", pre);
   else
      fprintf(fp, "\nobj = ATL_%cger1u_a1_x1_yX.o ATL_%cger1c_a1_x1_yX.o\n\n",
              pre, pre);

   fprintf(fp, "%clib : lib\nlib : %cr1k.grd\n", pre, pre);
   fprintf(fp, "%cr1k.grd : $(obj)\n", pre);
   fprintf(fp, "\t$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(obj)\n");
   fprintf(fp, "\t$(RANLIB) $(ATLASlib)\n");
   fprintf(fp, "\ttouch %cr1k.grd\n\n", pre);

   fprintf(fp, "clean :\n\t rm -f $(obj) %cr1k.grd\n", pre);
   fprintf(fp, "%cclean : clean\n\n", pre);
   if (pre == 's' || pre == 'd')
   {
      fprintf(fp, "ATL_%cger1_a1_x1_yX.o : ATL_%cger1.c\n", pre, pre);
      fprintf(fp, "\t$(R1CC) -o $@ -c $(R1FLAGS) ATL_%cger1.c\n", pre);
   }
   else
   {
      fprintf(fp, "ATL_%cger1u_a1_x1_yX.o : ATL_%cger1.c\n", pre, pre);
      fprintf(fp, "\t$(R1CC) -o $@ -c $(R1FLAGS) -D%s ATL_%cger1.c\n",
              typ, pre);
      fprintf(fp, "ATL_%cger1c_a1_x1_yX.o : ATL_%cger1.c\n", pre, pre);
      fprintf(fp, "\t$(R1CC) -o $@ -c $(R1FLAGS) -DConj_ ATL_%cger1.c\n",
              pre);
   }
   fclose(fp);
}

void r1install(char pre, int id, char *r1nam, char *auth, int l1mul,
               int flag, char *cc, char *ccflags, int mu, int nu)
{
   char ln[128];
   double mf;

   GenMakefile(pre, cc, ccflags);
   mf = r1case(pre, r1nam, flag, mu, nu, id, cc, ccflags, l1mul);
   emit_r1head(pre, l1mul, flag, mu, nu);
   sprintf(ln, "make %cinstall r1rout=%s\n", pre, r1nam);
   fprintf(stderr, "%s", ln);
   assert(system(ln) == 0);
   CreateSum(pre, l1mul, id, r1nam, auth, flag, mu, nu, mf);
}

int FindL1Mul(char pre, int cas, char *r1nam, int flag, int mu, int nu,
              char *cc, char *ccflags)
{
   double low = .5, high = 1.0;
   double mflow, mfhigh;
   int ilow, ihigh;

   if (ATL_NoBlock(flag)) flag -= 32;  /* always actually block these times */
   do
   {
      ilow = (low  * 100.0);
      ihigh = (high * 100.0);
      mflow  = r1case(pre, r1nam, flag, mu, nu, cas, cc, ccflags, ilow);
      mfhigh = r1case(pre, r1nam, flag, mu, nu, cas, cc, ccflags, ihigh);
      fprintf(stdout, "      %d%% %.2fMFLOP  ---  %d%% %.2fMFLOP\n",
              ilow, mflow, ihigh, mfhigh);
      if (mflow < mfhigh) low += 0.5*(high-low);
      else high -= 0.5 * (high-low);
   }
   while (ihigh-ilow);
   fprintf(stdout, "\n\nBEST %% of L1 cache: %.2f\n", low*100.0);
   return(ilow);
}

int ConfirmBlock(char pre, char *r1nam, int flag, int mu, int nu, int cas,
                 char *cc, char *ccflags, int l1mul)
{
   int bflag;
   double mfblock, mfnoblock;

   if ( ATL_NoBlock(flag) )
   {
      bflag = flag - 32;
      mfblock   = r1case(pre, r1nam, bflag, mu, nu, cas, cc, ccflags, l1mul);
      mfnoblock = r1case(pre, r1nam,  flag, mu, nu, cas, cc, ccflags, l1mul);
      fprintf(stdout, "\nWith blocking=%lf, without=%lf\n\n",
              mfblock, mfnoblock);
      if (mfblock >= mfnoblock) return(bflag);
   }
   return(flag);
}

void iGetIDCase(char pre, int ID, char *r1nam, char *auth,
                char *cc, char *ccflags, int *flag, int *mu, int *nu)
{
   char ln[512];
   FILE *fp;
   int i, id, n;

   sprintf(ln, "%ccases.dsc", pre);
   fp = fopen(ln, "r");
   assert(fp);

   assert(fgets(ln, 128, fp));
   assert(sscanf(ln, " %d", &n) == 1);
   for (i=0; i < n; i++)
   {
      assert(fgets(ln, 512, fp));
      assert(sscanf(ln, " %d %d %d %d %s \"%[^\"]",
                    &id, flag, mu, nu, r1nam, auth) == 6);
      if (LineIsCont(ln))
      {
         assert(fgets(cc, 256, fp));
         assert(fgets(ccflags, 512, fp));
         NoEndLineWhiteSpace(cc);
         NoEndLineWhiteSpace(ccflags);
      }
      else cc[0] = ccflags[0] = '\0';
      if (id == ID) break;
   }
   fclose(fp);
}

void GetCases(char pre, int *ncases, int **IDs, char ***r1nams, char ***auths,
              char ***CC, char ***CCFLAGS, int **flags, int **mus, int **nus)
{
   int n, i, j;
   char fnam[64], ln[512];
   char **r1nam, **auth, **cc, **ccflags;
   int *flag, *mu, *nu, *id;
   FILE *fp;

   sprintf(fnam, "%ccases.dsc", pre);
   fp = fopen(fnam, "r");
   assert(fp);

   assert(fgets(ln, 128, fp));
   assert(sscanf(ln, " %d", &n) == 1);
   assert(n < 100 && n > 0);

   r1nam   = malloc(n * sizeof(char*));
   auth    = malloc(n * sizeof(char*));
   cc      = malloc(n * sizeof(char*));
   ccflags = malloc(n * sizeof(char*));
   assert(r1nam && auth && cc && ccflags);
   for (i=0; i < n; i++)
   {
      assert(r1nam[i] = malloc(64*sizeof(char)));
      assert(auth[i] = malloc(64*sizeof(char)));
   }
   mu = malloc(n * sizeof(int));
   nu = malloc(n * sizeof(int));
   id = malloc(n * sizeof(int));
   flag = malloc(n * sizeof(int));
   assert(mu && nu && flag && id);
   for (i=0; i < n; i++)
   {
      assert(fgets(ln, 512, fp));
      assert(sscanf(ln, " %d %d %d %d %s \"%[^\"]",
                    id+i, flag+i, mu+i, nu+i, r1nam[i], auth[i]) == 6);
      assert(mu[i] >= 0 && nu[i] >= 0 && r1nam[i][0] != '\0');
      if (LineIsCont(ln))
      {
         fgets(ln, 512, fp);
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
   fclose(fp);
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
   *IDs = id;
   *ncases = n;
   *r1nams = r1nam;
   *auths = auth;
   *CC = cc;
   *CCFLAGS = ccflags;
   *flags = flag;
   *mus = mu;
   *nus = nu;
}

void RunCases(char pre)
{
   char **r1nams, **auths, **cc, **ccflags;
   int ncases, *ids, *flags, *mus, *nus;
   int i, imax=0;
   int l1mul, mf, mfmax=0.0;

   GetCases(pre, &ncases, &ids, &r1nams, &auths, &cc, &ccflags, &flags,
            &mus, &nus);

   for (i=0; i < ncases; i++)
   {
      mf = r1case(pre, r1nams[i], flags[i], mus[i], nus[i], ids[i],
                  cc[i], ccflags[i], 75);
      fprintf(stdout, "%s : %.2f\n", r1nams[i], mf);
      if (mf > mfmax)
      {
         if (r1tstcase(pre, r1nams[i], cc[i], ccflags[i]) == 0)
         {    /* make it pass tester */
            mfmax = mf;
            imax = i+1;
         }
      }
   }

   assert(imax);
   imax--;
   l1mul = FindL1Mul(pre, imax+1, r1nams[imax], flags[imax], mus[imax],
                     nus[imax], cc[imax], ccflags[imax]);
   flags[imax] = ConfirmBlock(pre, r1nams[imax], flags[imax], mus[imax],
                              nus[imax], ids[imax], cc[imax],
                              ccflags[imax], l1mul);
   fprintf(stdout, "\nBEST: %s, case %d, mu=%d, nu=%d; at %.2f MFLOPS\n\n",
           r1nams[imax], imax+1, mus[imax], nus[imax], mfmax);
   CreateSum(pre, l1mul, ids[imax], r1nams[imax], auths[imax], flags[imax],
             mus[imax], nus[imax], mfmax);

   for (i=0; i < ncases; i++)
   {
      free(r1nams[i]);
      free(auths[i]);
      free(cc[i]);
      free(ccflags[i]);
   }
   free(cc);
   free(ccflags);
   free(r1nams);
   free(auths);
   free(flags);
   free(mus);
   free(nus);
}

void GoToTown(pre)
{
   char r1nam[128], auth[128], ln[128], cc[256], ccflags[512];
   int flag, mu, nu, id, l1mul;
   double mf;
   FILE *fp;

   sprintf(ln, "res/%cR1RES", pre);
   if (!FileExists(ln))
   {
      RunCases(pre);
      fp = fopen(ln, "r");
      assert(fp);
   }
   else /* if default does not pass tester, rerun cases */
   {
      fp = fopen(ln, "r");
      ReadSum(pre, &l1mul, &id, r1nam, auth, &flag, &mu, &nu, &mf);
      iGetIDCase(pre, id, r1nam, auth, cc, ccflags, &flag, &mu, &nu);
      if (r1tstcase(pre, r1nam, cc, ccflags) != 0)
      {
         fclose(fp);
         assert(remove(ln) == 0);
         RunCases(pre);
         fp = fopen(ln, "r");
         assert(fp);
      }
   }
   fclose(fp);
   ReadSum(pre, &l1mul, &id, r1nam, auth, &flag, &mu, &nu, &mf);
   iGetIDCase(pre, id, r1nam, auth, cc, ccflags, &flag, &mu, &nu);
   fprintf(stdout, "\nBEST: %s, mu=%d, nu=%d; at %.2f MFLOPS\n\n",
           r1nam, mu, nu, mf);
   r1install(pre, id, r1nam, auth, l1mul, flag, cc, ccflags, mu, nu);
}

void PrintUsage(char *fnam)
{
   fprintf(stderr, "USAGE: %s [-p <d/s/z/c>]\n", fnam);
   exit(-1);
}

void GetFlags(int nargs, char **args, char *pre)
{
   char ctmp;
   int i;
   *pre = 'd';
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
      default:
         fprintf(stderr, "Unknown flag : %s\n", args[i]);
         PrintUsage(&args[0][0]);
      }
   }
}

main(int nargs, char **args)
{
   char pre;
   GetFlags(nargs, args, &pre);
   GoToTown(pre);
   exit(0);
}
