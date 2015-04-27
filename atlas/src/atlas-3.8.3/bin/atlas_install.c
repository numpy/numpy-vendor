/*
 *             Automatically Tuned Linear Algebra Software v3.8.3
 *                    (C) Copyright 1998 R. Clint Whaley
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
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "atlas_fopen.h"

char *redir="2>&1 | ./xatlas_tee";
char *fmake="make -f Makefile";
int DOSEARCH=1, REGS=32;
int QUERY=0;
FILE *fpI, *fparch;

#define Mciswspace(C) ( (((C) > 8) && ((C) < 14)) || ((C) == 32) )
#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )

#include <time.h>
void GetDate(int *month, int *day, int *year, int *hour, int *min)
{
   time_t tv;
   struct tm *tp;

   tv = time(NULL);
   tp = localtime(&tv);
   *month = tp->tm_mon + 1;
   *day = tp->tm_mday;
   *year = tp->tm_year + 1900;
   *hour = tp->tm_hour;
   *min = tp->tm_min;
}

long GetInt(FILE *fpin, long Default, char *spc, char *expstr)
/*
 * Gets a signed integral type from fpin.  If nothing or garbage is entered,
 * Default is returned.
 */
{
   char str[64];
   long iin;
   if (expstr) fprintf(stdout, "%sEnter %s [%d]: ", spc, expstr, Default);
   if (fgets(str, 64, fpin) == NULL) return(Default);
   if (sscanf(str, " %ld ", &iin) != 1) return(Default);
   return(iin);
}

long GetIntRange(long Default, long Min, long Max, char *spc, char *expstr)
{
   long i;
   int keepOn=0;
   do
   {
      i = GetInt(stdin, Default, spc, expstr);
      if (i > Max)
      {
         keepOn = 1;
         fprintf(stderr, "\n%d larger than max value of %d.  Try again.\n\n",
                 i, Max);
      }
      else if (i < Min)
      {
         keepOn = 1;
         fprintf(stderr, "\n%d smaller than min value of %d.  Try again.\n\n",
                 i, Min);
      }
      else keepOn = 0;
   }
   while (keepOn);
   return(i);
}

long GetIntVer(long Default, long Min, long Max, char *spc, char *expstr)
{
   long i, j;

   do
   {
      i = GetIntRange(Default, Min, Max, spc, expstr);
      fprintf(stdout, "%s   You entered: %d\n", spc, i);
      j = GetIntRange(0, 0, 1, spc, "1 to reenter, 0 accepts");
   }
   while(j);
   return(i);
}


void GetString(FILE *fpin, char *Default, char *spc, char *expstr,
               int len, char *str0)
/*
 * Get a string of length len, not including NULL terminator; pads
 * any extra len with NULLs
 */
{
   char str[512], *sp;
   int i;

   assert(len+1 <= 512);
   if (expstr)
   {
      if (Default) fprintf(stdout, "%sEnter %s [%s]: ", spc, expstr, Default);
      else fprintf(stdout, "%sEnter %s:", spc, expstr);
   }
   sp = fgets(str, 512, fpin);
   if ( (sp == NULL) || (str[0] == '\0') || (str[0] == '\n') )
   {
      if (Default) strcpy(str0, Default);
      else str0[0] = '\0';
      return;
   }
   str[len] = '\0';
   for (i=0; str0[i] = str[i]; i++);
   if (i) i--;
   while (Mciswspace(str0[i])) i--;
   while (++i < len) str0[i] = '\0';
   str0[i] = '\0';
}

void GetStrVer(char *def, char *spc, char *expstr, int len, char *str)
{
   int i;

   do
   {
      GetString(stdin, def, spc, expstr, len, str);
      fprintf(stdout, "%sYou have entered '%s'\n", spc, str);
      i = GetIntRange(0, 0, 1, spc, "1 to reenter, 0 to accept");
   }
   while(i);
}

int IsYes(char def, char *spc, char *expstr)
{
   char ch, ln[256];
   fprintf(stdout, "%s%s [%c]: ", spc, expstr, def);
   if (fgets(ln, 256, stdin) == NULL) ch=def;
   else if (ln[0] == '\0' || ln[0] == '\n') ch=def;
   else ch = ln[0];
   return( ((ch == 'y') || (ch == 'Y')) );
}

char GetChar(char def, char *spc, char *expstr)
{
   char ch, ln[256];
   fprintf(stdout, "%s%s [%c]: ", spc, expstr, def);
   if (fgets(ln, 256, stdin) == NULL) ch=def;
   else if (ln[0] == '\0' || ln[0] == '\n') ch=def;
   else ch = ln[0];
   return(ch);
}

int FileIsThere(char *nam)
{
   FILE *fp;

   fp = fopen(nam, "r");
   if (fp == NULL) return(0);
   fclose(fp);
   return(1);
}

#include <stdarg.h>
#define ATL_UseStringVarArgs
void ATL_mprintf(int np, ...)
/*
 * Prints same message to np output file streams
 */
{
   va_list argptr;
   FILE *fp[8];
   char *form;
   int i;
   #ifdef ATL_UseStringVarArgs
      char ln[1024];
   #endif

   if (np > 0)
   {
      va_start(argptr, np);
      assert(np <= 8);
      for (i=0; i < np; i++) fp[i] = va_arg(argptr, FILE *);
      form = va_arg(argptr, char *);
      #ifdef ATL_UseStringVarArgs
         vsprintf(ln, form, argptr);
         assert(strlen(ln) < 1024);/* sanity test only, will not stop overrun */
         va_end(argptr);
         for (i=0; i < np; i++) if (fp[i]) fprintf(fp[i], ln);
      #else
         for (i=0; i < np; i++) if (fp[i]) vfprintf(fp[i], form, argptr);
         va_end(argptr);
      #endif
   }
}

int GetFirstInt(char *ln)
{
   int i, iret=0;
   for (i=0; ln[i]; i++)
   {
      if (isdigit(ln[i]))
      {
         sscanf(ln+i, "%d", &iret);
         break;
      }
   }
   return(iret);
}

long long GetFirstLong(char *ln)
{
   int i;
   long long iret=0;
   for (i=0; ln[i]; i++)
   {
      if (isdigit(ln[i]))
      {
         sscanf(ln+i, "%ld", &iret);
         break;
      }
   }
   return(iret);
}

double GetFirstDouble(char *ln)
/*
 * Gets a double, which begins wt digit or "." (i.e., won't get e10)
 */
{
   int i;
   double dret=0;
   for (i=0; ln[i]; i++)
   {
      if (isdigit(ln[i]))
      {
         if (i > 0 && ln[i-1] == '.') i--;
         sscanf(ln+i, "%lf", &dret);
         break;
      }
   }
   return(dret);
}

int GetLastInt(char *ln)
{
   int i, iret=0;
   for (i=0; ln[i]; i++);
   if (i > 0) for (i--; i > 0 && !isdigit(ln[i]); i--);
   if (i > 0 || (i == 0 && isdigit(ln[0])))
   {
      while(isdigit(ln[i]) && i > 0) i--;
      if (!isdigit(ln[i])) i++;
      sscanf(ln+i, "%d", &iret);
   }
   return(iret);
}

long long GetLastLong(char *ln)
{
   int i;
   long iret=0;
   for (i=0; ln[i]; i++);
   if (i > 0) for (i--; i > 0 && !isdigit(ln[i]); i--);
   if (i > 0 || (i == 0 && isdigit(ln[0])))
   {
      while(isdigit(ln[i]) && i > 0) i--;
      if (!isdigit(ln[i])) i++;
      sscanf(ln+i, "%ld", &iret);
   }
   return(iret);
}

void PrintBanner(char *fnam, int START, int sec, int subsec, int subsubsec);
static void ATL_Cassert0(int cond, char *exp, char *logfile, int line)
{
   FILE *fperr;

   if (cond) return;

   fperr = fopen("INSTALL_LOG/ERROR.LOG", "a");
   if (logfile)
      ATL_mprintf(2, stderr, fperr,
                  "ERROR %d DURING %s!!.  CHECK %s FOR DETAILS.\n",
                  line, exp, logfile);
   else ATL_mprintf(2, stderr, fperr, "ERROR %d DURING %s!!.\n", line, exp);
   if (system("make error_report") == 0)
   {
      ATL_mprintf(2, stderr, fperr,
"Error report error_<ARCH>.tgz has been created in your top-level ATLAS\n");
      ATL_mprintf(2, stderr, fperr,
         "directory.  Be sure to include this file in any help request.\n");
   }
   if (fperr) fclose(fperr);
   system("cat ../../CONFIG/error.txt >> INSTALL_LOG/ERROR.LOG");
   system("cat ../../CONFIG/error.txt");
   exit(-1);
}
#define ATL_Cassert(cond_, exp_, logfile_) \
    ATL_Cassert0((int)(cond_), exp_, logfile_, __LINE__)

void GetInstLogFile(char *nam, int pre, int *muladd, int *pf, int *lat,
                    int *nb, int *mu, int *nu, int *ku, int *ForceFetch,
                    int *ifetch, int *nfetch, double *mflop)
{
   char ln[128];
   FILE *fp;

   fp = fopen(nam, "r");
   if (fp == NULL) fprintf(stderr, "file %s not found!!\n\n", nam);
   assert(fp);
   fgets(ln, 128, fp);
   fscanf(fp, " %d %d %d %d %d %d %d %d %d %d %lf\n",
          muladd, lat, pf, nb, mu, nu, ku, ForceFetch, ifetch, nfetch, mflop);
   fclose(fp);
}

void PrintBanner(char *fnam, int START, int sec, int subsec, int subsubsec)
{
   int month, day, year, hour, min;
   char *sep="*******************************************************************************\n";
   char ln[256];
   FILE *fp;

   sprintf(ln, "ERROR OPENING FILE %s\n", fnam);
   fp = fopen(fnam, "a");
   if (fp == NULL) return;
   GetDate(&month, &day, &year, &hour, &min);
   fprintf(fp, "\n%s%s%s", sep, sep, sep);
   if (START)
      fprintf(fp, "*       BEGAN ATLAS3.8.3  INSTALL OF SECTION %1d-%1d-%1d ON %02d/%02d/%04d AT %02d:%02d     *\n",
              sec, subsec, subsubsec, month, day, year, hour, min);
   else
      fprintf(fp, "*      FINISHED ATLAS3.8.3  INSTALL OF SECTION %1d-%1d-%1d ON %02d/%02d/%04d AT %02d:%02d   *\n",
                      sec, subsec, subsubsec, month, day, year, hour, min);
   fprintf(fp, "%s%s%s\n\n\n", sep, sep, sep);
   fclose(fp);
}

void GetMMRES(char pre, int *muladd, int *lat, int *nb, int *pref,
              int *mu, int *nu, int *ku, int *ff, int *iff, int *nf,
              double *mf, int *icase, char *ufile, char *auth, double *umf)
{
   int h, i, j, k;
   char ln[256];
   FILE *fp;

   sprintf(ln, "INSTALL_LOG/%cMMRES", pre);
   fp = fopen(ln, "r");
   assert(fp);
   fgets(ln, 256, fp);
   fgets(ln, 256, fp);
   assert(sscanf(ln, " %d %d %d %d %d %d %d %d %d %d %lf", muladd, lat, pref,
                 nb, mu, nu, ku, ff, iff, nf, mf) == 11);
   fgets(ln, 256, fp);
   if ( fgets(ln, 256, fp) )  /* user-supplied GEMM was best */
   {
      assert(fscanf(fp, " %d %d %lf \"%[^\"]\" \"%[^\"]",
                    &icase, &i, umf, ufile, auth) == 5);
   }
   else
   {
      ufile[0] = auth[0] = '\0';
      *umf = 0.0;
      *icase = -1;
   }
   fclose(fp);
}

void GoToTown(int ARCHDEF, int L1DEF)
{
   const char TR[2] = {'N','T'};
   char prec[4] = {'d', 's', 'z', 'c'}, pre, upre, *typ;
   char ln[512], tnam[256], ln2[512], ln3[512], fnam[128];
   char *mulinst, *peakstr, *peakstr2;
   int nprec=4;
   int iL1, lat, muladd, maused, latuse, lbnreg;
   int len, i, j, ip, ia, ib, ncnb, pf;
   int FoundCE=0;
   int maxreg, latus, nb, mu, nu, ku;
   int ffetch, ifetch, nfetch;
   int DefInstall=0;
   long imf;
   int maU, latU, muU, nuU, kuU, il1mul, pfA;
   double mfU, mf4x1, mf4x4, mf, mfp, mmmf, mfpeak[2], l1mul;
   FILE *fp, *fpsum, *fpabr;

   PrintBanner("INSTALL_LOG/SUMMARY.LOG", 1, 0, 0, 0);
   fpsum = fopen("INSTALL_LOG/SUMMARY.LOG", "a");
   ATL_Cassert(fpsum, "OPENING INSTALL_LOG/SUMMARY.LOG", NULL);

   ATL_Cassert(tmpnam(tnam), "GETTING TEMPFILE", NULL);

   if (L1DEF)
   {
      ATL_Cassert(system("make IBozoL1.grd\n")==0,
                  "USING BOZO L1 DEFAULTS", NULL);
   }
   if (ARCHDEF)
      DefInstall = !system("make IArchDef.grd\n");

   ATL_mprintf(2, stdout, fpsum,
               "\n\nIN STAGE 1 INSTALL:  SYSTEM PROBE/AUX COMPILE\n");

   sprintf(ln2, "INSTALL_LOG/Stage1.log");
   PrintBanner("INSTALL_LOG/Stage1.log", 1, 1, 0, 0);
   sprintf(ln, "%s IStage1 %s INSTALL_LOG/Stage1.log\n", fmake, redir);
   ATL_Cassert(system(ln)==0, "Stage 1 install", ln2);
   fp = fopen("INSTALL_LOG/L1CacheSize", "r");

   ATL_Cassert(fp, "CACHESIZE SEARCH", ln2);
   ATL_Cassert(fscanf(fp, "%d", &i) == 1, "CACHESIZE SEARCH", ln2);
   fclose(fp);
   fprintf(stdout, "\n\n   Level 1 cache size calculated as %dKB\n", i);
   fprintf(fpsum, "   Level 1 cache size calculated as %dKB.\n\n", i);
   iL1 = i;

   for (ip=0; ip < 2; ip++)
   {
      sprintf(ln, "INSTALL_LOG/%cMULADD", prec[ip]);
      fp = fopen(ln, "r");
      ATL_Cassert(fp, "FPU PROBE", NULL);
      ATL_Cassert(fscanf(fp, " %d", &muladd)==1, "FPU PROBE", ln2);
      ATL_Cassert(fscanf(fp, " %d", &lat)==1, "FPU PROBE", ln2);
      ATL_Cassert(fscanf(fp, " %lf", &mfpeak[ip])==1, "FPU PROBE", ln2);
      ATL_Cassert(fscanf(fp, " %d", &lbnreg)==1, "FPU PROBE", ln2);
      fclose(fp);
      if (muladd) mulinst = "Combined muladd instruction";
      else mulinst = "Separate multiply and add instructions";
      ATL_mprintf(2, stdout, fpsum, "   %cFPU: %s with %d cycle pipeline.\n",
                  prec[ip], mulinst, lat);
      ATL_mprintf(2, stdout, fpsum,
                  "         Apparent number of registers : %d\n", lbnreg);
      ATL_mprintf(2, stdout, fpsum,
                  "         Register-register performance=%.2fMFLOPS\n",
                  mfpeak[ip]);
   }

   fprintf(stdout, "\n\nIN STAGE 2 INSTALL:  TYPE-DEPENDENT TUNING\n");
   fprintf(fpsum , "\n\nIN STAGE 2 INSTALL:  TYPE-DEPENDENT TUNING\n");

   for (ip=0; ip < nprec; ip++)
   {
      sprintf(ln, "INSTALL_LOG/%cPerfSumm.txt", prec[ip]);
      fpabr = fopen(ln, "w");
      ATL_Cassert((fpabr != NULL), "Unable to open abbreviation file", NULL);
      pre = prec[ip];
      switch(pre)
      {
      case 's':
         len = sizeof(float);
         typ = "SREAL";
         upre = 's';
         break;
      case 'd':
         len = sizeof(double);
         typ = "DREAL";
         upre = 'd';
         break;
      case 'c':
         len = sizeof(float);
         typ = "SCPLX";
         upre = 's';
         break;
      case 'z':
         len = sizeof(double);
         typ = "DCPLX";
         upre = 'd';
         break;
      }

      fprintf(stdout,
              "\n\nSTAGE 2-%d: TUNING PREC=\'%c\' (precision %d of %d)\n",
              ip+1, pre, ip+1, nprec);
      fprintf(fpsum ,
              "\n\nSTAGE 2-%d: TUNING PREC=\'%c\' (precision %d of %d)\n",
              ip+1, pre, ip+1, nprec);
      fprintf(stdout, "\n\n   STAGE 2-%d-1 : BUILDING BLOCK MATMUL TUNE\n",
              ip+1);
      fprintf(fpsum , "\n\n   STAGE 2-%d-1 : BUILDING BLOCK MATMUL TUNE\n",
              ip+1);
/*
 *    If necessary, install matmul for this type
 */
      sprintf(fnam, "INSTALL_LOG/%cMMRES", pre);
      if (!FileExists(fnam))  /* need to run search or make link */
      {
         sprintf(ln2, "INSTALL_LOG/%cMMSEARCH.LOG", pre);
         PrintBanner(ln2, 1, 2, ip+1, 1);
         if (DefInstall)
         {
            sprintf(ln, "%s IRunMMDef pre=%c %s %s\n", fmake, pre, redir, ln2);
            fprintf(stdout, ln);
            ATL_Cassert(system(ln)==0, "BUILDING BLOCK MATMUL TUNE", ln2);
         }
         sprintf(ln, "%s %s pre=%c %s %s\n", fmake, fnam, pre, redir, ln2);
         fprintf(stdout, ln);
         ATL_Cassert(system(ln)==0, "BUILDING BLOCK MATMUL TUNE", ln2);
         PrintBanner(ln2, 0, 2, ip+1, 1);
      }
      GetMMRES(pre, &muladd, &lat, &nb, &pfA, &mu, &nu, &ku, &ffetch, &ifetch,
               &nfetch, &mf, &i, fnam, ln, &mfU);
      #ifdef ATL_CPUMHZ
         fprintf(fpabr, "Clock_rate=%d Mhz\n", ATL_CPUMHZ);
         mfp = ATL_CPUMHZ;
         peakstr = "of detected clock rate";
         peakstr2 = "clock";
      #else
         fprintf(fpabr, "Clock_rate=0 Mhz\n");
         if (pre == 'd' || pre == 'z') mfp = mfpeak[0];
         else mfp = mfpeak[1];
         peakstr = "of apparent peak";
         peakstr2 = "peak";
      #endif
      fprintf(fpabr, "%% clock      MFLOP  ROUTINE/PROBLEM\n");
      fprintf(fpabr, "=======  =========  ===============\n");
      if (mfU > mf)
      {
         ATL_mprintf(2, fpsum, stdout,
            "      The best matmul kernel was %s, NB=%d, written by %s\n",
                    fnam, nb, ln);
         ATL_mprintf(2, fpsum, stdout,
            "      Performance: %.2fMFLOPS (%5.2f percent of %s)\n",
                     mfU, (mfU/mfp)*100.0, peakstr);
         ATL_mprintf(2, fpsum, stdout,
                     "        (Gen case got %.2fMFLOPS)\n", mf);
         mmmf = mfU;
      }
      else
      {
         ATL_mprintf(2, fpsum, stdout,
"      %cL1MATMUL: lat=%d, nb=%d, pf=%d, mu=%d, nu=%d, ku=%d, if=%d, nf=%d;\n",
                 pre, lat, nb, pfA, mu, nu, ku, ifetch, nfetch);
         ATL_mprintf(2, fpsum, stdout,
            "                 Performance: %.2f (%5.2f percent of %s)\n",
                     mf, (mf/mfp)*100.0, peakstr);
         mmmf = mf;
      }
      fprintf(fpabr, "%7.1f %10.1f  %s\n", (mmmf/mfp)*100.0, mmmf,
              "Chosen kgemm");
      fprintf(fpabr, "%7.1f %10.1f  %s\n", (mf/mfp)*100.0, mf,
              "Generated kgemm");

      sprintf(fnam, "INSTALL_LOG/%cNCNB", pre);
      if (!FileExists(fnam))
      {
         sprintf(ln, "%s %s pre=%c %s %s", fmake, fnam, pre, redir, ln2);
         fprintf(stdout, ln);
         ATL_Cassert(system(ln)==0, "BUILDING BLOCK MATMUL TUNE", ln2);
      }
      fp = fopen(fnam, "r");
      ATL_Cassert(fp, "OPENING NCNB", NULL);
      ATL_Cassert(fscanf(fp, " %d", &ncnb) == 1, "READING NCNB", NULL);
      fclose(fp);

      for (ia=0; ia < 2; ia++)
      {
         for (ib=0; ib < 2; ib++)
         {
            sprintf(fnam, "INSTALL_LOG/%cbest%c%c_%dx%dx%d", pre, TR[ia], TR[ib],
                    ncnb, ncnb, ncnb);
            if (!FileExists(fnam))
            {
               sprintf(ln, "%s %s pre=%c nb=%d %s %s",
                       fmake, fnam, pre, ncnb, redir, ln2);
               fprintf(stdout, ln);
               ATL_Cassert(system(ln)==0, "BUILDING BLOCK MATMUL TUNE", ln2);
            }
            GetInstLogFile(fnam, pre, &muladd, &pf, &lat, &nb, &mu, &nu, &ku,
                           &ffetch, &ifetch, &nfetch, &mf);
            fprintf(stdout,
   "      NCgemm%c%c : muladd=%d, lat=%d, pf=%d, nb=%d, mu=%d, nu=%d ku=%d,\n",
                    TR[ia], TR[ib], muladd, lat, pf, nb, mu, nu, ku);
            fprintf(stdout,"                 ForceFetch=%d, ifetch=%d nfetch=%d\n",
                    ffetch, ifetch, nfetch);
            fprintf(stdout,
"                 Performance = %.2f (%5.2f of copy matmul, %5.2f of %s)\n",
                    mf, (mf/mmmf)*100.0, (mf / mfp)*100.0, peakstr2);
            fprintf(fpsum,
"      mm%c%c   : ma=%d, lat=%d, nb=%d, mu=%d, nu=%d ku=%d, ff=%d, if=%d, nf=%d\n",
                    TR[ia], TR[ib], muladd, lat, nb, mu, nu, ku,
                    ffetch, ifetch, nfetch);
            fprintf(fpsum,
"               Performance = %.2f (%5.2f of copy matmul, %5.2f of %s)\n",
                    mf, (mf/mmmf)*100.0, (mf / mfp)*100.0, peakstr2);
            if (ia != ib)
               fprintf(fpabr, "%7.1f %10.1f  kgemm%c%c\n", (mf/mfp)*100.0, mf,
                       TR[ia], TR[ib]);
         }
      }

      sprintf(ln, "%s MMinstall pre=%c %s %s\n", fmake, pre, redir, ln2);
      fprintf(stdout, ln);
      ATL_Cassert(system(ln)==0, "BUILDING BLOCK MATMUL TUNE", ln2);

      fprintf(fpsum, "\n");
/*
 *    If necessary, find cacheedge
 */
      ATL_mprintf(2, fpsum, stdout,
                 "\n\n   STAGE 2-%d-2: CacheEdge DETECTION\n", ip+1);
      sprintf(ln2, "INSTALL_LOG/%cMMCACHEEDGE.LOG", pre);
      if (!FileExists("INSTALL_LOG/atlas_cacheedge.h"))
      {
         PrintBanner(ln2, 1, 2, ip+1, 2);
         sprintf(ln, "%s INSTALL_LOG/atlas_cacheedge.h pre=%c %s %s\n",
                 fmake, pre, redir, ln2);
         fprintf(stdout, ln);
         ATL_Cassert(system(ln)==0, "CACHEEDGE DETECTION", ln2);
         PrintBanner(ln2, 0, 2, ip+1, 2);
      }
      fp = fopen("INSTALL_LOG/atlas_cacheedge.h", "r");
      ATL_Cassert(fp, "CACHE EDGE DETECTION", NULL);
      ATL_Cassert(fgets(ln, 256, fp), "CACHE EDGE DETECTION", NULL);
      ATL_Cassert(fgets(ln, 256, fp), "CACHE EDGE DETECTION", NULL);
      ATL_Cassert(fgets(ln, 256, fp), "CACHE EDGE DETECTION", NULL);
      if (fgets(ln3, 256, fp))
      {
         ATL_Cassert(sscanf(ln+21, " %d", &i)==1, "CACHE EDGE DETECTION", NULL);
      }
      else i = 0;
      fprintf(fpsum, "      CacheEdge set to %d bytes\n", i);
      fclose(fp);
/*
 *    Determine [ZD,CS]NKB, if necessary
 */
      if (pre == 'z' || pre == 'c')
      {
         sprintf(ln3, "INSTALL_LOG/atlas_%c%cNKB.h", pre, upre);
         if (!FileExists(ln3))
         {
            sprintf(ln, "%s %s pre=%c %s %s\n",
                    fmake, ln3, pre, redir, ln2);
            fprintf(stdout, ln);
            ATL_Cassert(system(ln)==0, "CACHEEDGE DETECTION", ln2);
         }
         fp = fopen(ln3, "r");
         ATL_Cassert(fp, "CACHE EDGE DETECTION", NULL);
         ATL_Cassert(fgets(ln, 256, fp), "CACHE EDGE DETECTION", NULL);
         ATL_Cassert(fgets(ln, 256, fp), "CACHE EDGE DETECTION", NULL);
         ATL_Cassert(fgets(ln, 256, fp), "CACHE EDGE DETECTION", NULL);
         if (fgets(ln3, 256, fp))
         {
            ATL_Cassert(sscanf(ln+21, " %d", &i)==1,
                               "CACHE EDGE DETECTION", NULL);
         }
         else i = 0;
         fprintf(fpsum, "      %c%cNKB set to %d bytes\n", pre, upre, i);
         fclose(fp);
      }
/*
 *    If necessary, determine Xover for this data type
 */
      ATL_mprintf(2, fpsum, stdout,
                 "\n\n   STAGE 2-%d-3: LARGE/SMALL CASE CROSSOVER DETECTION\n",
                 ip+1);
      sprintf(fnam, "INSTALL_LOG/%cXover.h", pre);
      if (!FileExists(fnam))  /* need to run Xover tests */
      {
         sprintf(ln2, "INSTALL_LOG/%cMMCROSSOVER.LOG", pre);
         PrintBanner(ln2, 1, 2, ip+1, 3);
            fprintf(stdout,
              "\n\n   STAGE 2-%d-3: COPY/NO-COPY CROSSOVER DETECTION\n", ip+1);
         fprintf(fpsum,
              "\n\n   STAGE 2-%d-3: COPY/NO-COPY CROSSOVER DETECTION\n", ip+1);

         sprintf(ln, "%s %s pre=%c %s %s\n", fmake, fnam, pre, redir, ln2);
         fprintf(stdout, ln);
         ATL_Cassert(system(ln)==0, "COPY/NO-COPY CROSSOVER DETECTION", ln2);
         PrintBanner(ln2, 0, 2, ip+1, 3);
         fprintf(stdout, "      done.\n");
         fprintf(fpsum , "      done.\n");
      }

      sprintf(ln2, "INSTALL_LOG/%cL3TUNE.LOG", pre);
      PrintBanner(ln2, 1, 2, ip+1, 5);
      ATL_mprintf(2, stdout, fpsum,
                 "\n\n   STAGE 2-%d-4: LEVEL 3 BLAS TUNE\n", ip+1);
      if (pre == 's' || pre == 'd')
      {
         sprintf(ln, "%s INSTALL_LOG/atlas_%ctrsmXover.h pre=%c %s %s\n",
                 fmake, pre, pre, redir, ln2);
         fprintf(stdout, ln);
         ATL_Cassert(system(ln)==0, "L3BLAS TUNING", ln2);
      }
      else
      {
         sprintf(ln, "%s Il3lib pre=%c %s %s\n", fmake, pre, redir, ln2);
         fprintf(stdout, ln);
         ATL_Cassert(system(ln)==0, "L3BLAS TUNING", ln2);
      }
      sprintf(ln, "%s %ccblaslib %s %s\n", fmake, pre, redir, ln2); /* cblas */
      fprintf(stdout, ln);
      ATL_Cassert(system(ln)==0, "L3BLAS TUNING", ln2);
      PrintBanner(ln2, 0, 2, ip+1, 5);
      ATL_mprintf(2, fpsum, stdout, "      done.\n");


      ATL_mprintf(2, fpsum, stdout,"\n\n   STAGE 2-%d-5: GEMV TUNE\n", ip+1);
      sprintf(fnam, "INSTALL_LOG/%cMVRES", pre);
      if (!FileExists(fnam))
      {
         sprintf(ln2, "INSTALL_LOG/%cMVTUNE.LOG", pre);
         PrintBanner(ln2, 1, 2, ip+1, 7);
         sprintf(ln, "%s %s pre=%c %s %s\n", fmake, fnam, pre, redir, ln2);
         fprintf(stdout, ln);
         ATL_Cassert(system(ln)==0, "MVTUNE", ln2);
         ATL_Cassert(FileIsThere(fnam), "MVTUNE", ln2);
         PrintBanner(ln2, 0, 2, ip+1, 7);
      }
      fp = fopen(fnam, "r");
      ATL_Cassert(fp, "MVTUNE", NULL);
      ATL_Cassert(fgets(ln, 256, fp), "MVTUNE", NULL);
      ATL_Cassert(sscanf(ln, " %d %d %d %d %d %lf %s \"%[^\"]", &j, &i,
                  &mu, &nu, &il1mul, &mf, ln3, tnam) == 8, "MVTUNE", NULL);
      ATL_mprintf(2, fpsum, stdout,
                  "      gemvN : chose routine %d:%s written by %s\n",
                  j, ln3, tnam);
      ATL_mprintf(2, fpsum, stdout,
         "              Yunroll=%d, Xunroll=%d, using %d percent of L1\n",
                  mu, nu, il1mul);
      ATL_mprintf(2, fpsum, stdout,
"              Performance = %.2f (%5.2f of copy matmul, %5.2f of %s)\n",
              mf, (mf/mmmf)*100.0, (mf / mfp)*100.0, peakstr2);
      fprintf(fpabr, "%7.1f %10.1f  %s\n", (mf/mfp)*100.0, mf, "kgemvN");
      ATL_Cassert(fgets(ln, 256, fp), "MVTUNE", NULL);
      ATL_Cassert(sscanf(ln, " %d %d %d %d %d %lf %s \"%[^\"]", &j, &i,
                  &mu, &nu, &il1mul, &mf, ln3, tnam) == 8, "MVTUNE", NULL);
      fclose(fp);
      ATL_mprintf(2, fpsum, stdout,
                  "      gemvT : chose routine %d:%s written by %s\n",
                  j, ln3, tnam);
      ATL_mprintf(2, fpsum, stdout,
         "              Yunroll=%d, Xunroll=%d, using %d percent of L1\n",
                  mu, nu, il1mul);
      ATL_mprintf(2, fpsum, stdout,
"              Performance = %.2f (%5.2f of copy matmul, %5.2f of %s)\n",
              mf, (mf/mmmf)*100.0, (mf / mfp)*100.0, peakstr2);
      fprintf(fpabr, "%7.1f %10.1f  %s\n", (mf/mfp)*100.0, mf, "kgemvT");

      ATL_mprintf(2, fpsum, stdout,"\n\n   STAGE 2-%d-6: GER TUNE\n", ip+1);
      sprintf(fnam, "INSTALL_LOG/%cR1RES", pre);
      if (!FileExists(fnam))
      {
         sprintf(ln2, "INSTALL_LOG/%cR1TUNE.LOG", pre);
         PrintBanner(ln2, 1, 2, ip+1, 7);
         sprintf(ln, "%s %s pre=%c %s %s\n", fmake, fnam, pre, redir, ln2);
         fprintf(stdout, ln);
         ATL_Cassert(system(ln)==0, "R1TUNE", ln2);
         ATL_Cassert(FileIsThere(fnam), "R1TUNE", ln2);
         PrintBanner(ln2, 0, 2, ip+1, 7);
      }
      fp = fopen(fnam, "r");
      ATL_Cassert(fp, "R1TUNE", NULL);
      ATL_Cassert(fgets(ln, 256, fp), "R1TUNE", NULL);
      ATL_Cassert(sscanf(ln, " %d %d %d %d %d %lf %s \"%[^\"]", &j, &i, &mu,
                         &nu, &il1mul, &mf, ln3, tnam) == 8, "R1TUNE", NULL);
      fclose(fp);
      ATL_mprintf(2, fpsum, stdout,
              "      ger : chose routine %d:%s written by %s\n", j, ln3, tnam);
      ATL_mprintf(2, fpsum, stdout,
              "            mu=%d, nu=%d, using %5.2f percent of L1 Cache\n",
              mu, nu, il1mul*0.01);
      ATL_mprintf(2, fpsum, stdout,
"              Performance = %.2f (%5.2f of copy matmul, %5.2f of %s)\n",
              mf, (mf/mmmf)*100.0, (mf / mfp)*100.0, peakstr2);
      fprintf(fpabr, "%7.1f %10.1f  %s\n", (mf/mfp)*100.0, mf, "kger");
      fclose(fpabr);
   }
   ATL_mprintf(2, fpsum, stdout,"\n\nSTAGE 3: GENERAL LIBRARY BUILD\n");

   sprintf(ln2, "INSTALL_LOG/LIBBUILD.LOG");
   PrintBanner(ln2, 1, 3, 1, 1);
   sprintf(ln, "%s IBuildLibs %s %s\n", fmake, redir, ln2);
   fprintf(stdout, ln);
   ATL_Cassert(system(ln)==0, "LIBRARY BUILD", ln2);
   ATL_Cassert(FileIsThere(fnam), "LIBRARY BUILD", ln2);
   PrintBanner(ln2, 0, 3, 1, 1);

   ATL_mprintf(2, fpsum, stdout,"\n\nSTAGE 4: POST-BUILD TUNING\n");
   sprintf(ln2, "INSTALL_LOG/POSTTUNE.LOG");
   PrintBanner(ln2, 1, 4, 1, 1);
   sprintf(ln, "%s IPostTune %s %s\n", fmake, redir, ln2);
   fprintf(stdout, ln);
   ATL_Cassert(system(ln)==0, "POST-BUILD TUNE", ln2);
   ATL_Cassert(FileIsThere(fnam), "POST-BUILD TUNE", ln2);
   PrintBanner(ln2, 0, 4, 1, 1);

   fprintf(stdout, "   done.\n");
   fprintf(fpsum,  "   done.\n");

#ifdef ATL_NCPU
   ATL_mprintf(2, stdout, fpsum, "\n\nSTAGE 4: Threading install\n");

   sprintf(ln2, "INSTALL_LOG/LIBPTBUILD.LOG");
   PrintBanner(ln2, 1, 4, 2, 1);
   sprintf(ln, "%s IBuildPtlibs %s %s\n", fmake, redir, ln2);
   fprintf(stdout, ln);
   ATL_Cassert(system(ln)==0, "PTLIBRARY BUILD", ln2);
   PrintBanner(ln2, 0, 4, 2, 1);
#endif

   fprintf(stdout, "\n\n\n\n");
   fprintf(stdout, "ATLAS install complete.  Examine \n");
   fprintf(stdout, "ATLAS/bin/<arch>/INSTALL_LOG/SUMMARY.LOG for details.\n");
   fclose(fpsum);
   PrintBanner("INSTALL_LOG/SUMMARY.LOG", 0, 0, 0, 0);
}

void PrintUsage(char *nam)
{
   fprintf(stderr, "\n\nUSAGE: %s [-a <#(use archdef)> -1 <#(bozoL1)>]\n\n",
           nam);
   exit(-1);
}

void GetFlags(int nargs, char *args[], int *ARCHDEF, int *L1DEF)
{
   int i;

   *L1DEF = 0;
   *ARCHDEF = 1;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'a':
         if (++i > nargs)
            PrintUsage(args[0]);
         *ARCHDEF = atoi(args[i]);
         break;
      case '1':
         if (++i > nargs)
            PrintUsage(args[0]);
         *L1DEF = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0]);
      }
   }
}

main(int nargs, char *args[])
{
   int L1DEF, ARCHDEF;
   GetFlags(nargs, args, &ARCHDEF, &L1DEF);
   GoToTown(ARCHDEF, L1DEF);
   exit(0);
}
