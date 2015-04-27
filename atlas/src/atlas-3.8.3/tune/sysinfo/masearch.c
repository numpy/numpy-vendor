/*
 *             Automatically Tuned Linear Algebra Software v3.8.3
 *                    (C) Copyright 1997 R. Clint Whaley
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
#include <assert.h>

#define NTIM 3
#define TOLERANCE 1.2
#define BETTER 1.04  /* require at least 4% improvement, see below */

/*
 * We require at least 4% improvement before declaring a longer latency to
 * be better.  This is because longer latencies have less loop overhead, which
 * gives them a possible advantage, depending on the unrollings the compiler
 * does.  When you actually have too small a latency, this event should be
 * relatively catastrophic, with drops greater > 10%, so this is relatively
 * safe assumption.  Also, the smaller latencies require less registers to
 * satisfy, so they tend to be preferable.  - Clint
 */

void emit_muladd(char *type, char pre, int MULADD, int lat,
                 unsigned long nlma,  /* # muladds in loop */
                 unsigned long mflop)    /* n mflops to time */
{
   int i, j, nama;
   char *ma;
   unsigned long reps;
   double mf;
   FILE *fpout;

   nama = (nlma + lat-1) / lat;  /* # of actual muladds in loop */
   mflop *= 1000000;
   reps = mflop / (2*nama);

   fpout = fopen("muladd.c", "w");
   assert(fpout != NULL);
   if (MULADD) ma = "Combined MULADD";
   else ma = "Separate multiply and add";
   fprintf(fpout,
           "#include <stdio.h>\n#include <stdlib.h>\n#include <assert.h>\n");
   fprintf(fpout, "#if defined(PentiumCPS) || defined(WALL)\n");
   fprintf(fpout, "   #define time00 ATL_walltime\n");
   fprintf(fpout, "#else\n   #define time00 ATL_cputime\n#endif\n");
   fprintf(fpout, "double time00();\n");
/*
 * function macase()
 */
   fprintf(fpout, "\nstatic double macase(long nreps, int PRINT)\n{\n");
   fprintf(fpout, "   long i = nreps;\n");
   fprintf(fpout, "   double t0, tim, mf;\n");
   fprintf(fpout, "   register %s c0", type);
   for (i=1; i < lat; i++) fprintf(fpout, ", c%d", i);
   if (!MULADD) for (i=0; i < lat; i++) fprintf(fpout, ", m%d", i);
   fprintf(fpout, ";\n\n");
   fprintf(fpout,
      "   if (!(rand()|rand()|rand()|rand()|rand()|rand())) nreps = rand();\n");

/*
 * Try to make sure compiler cannot figure out ops can be reduced to asg
 */
   for (i=0; i < lat; i++)
      fprintf(fpout,
      "   if (rand()|rand()|rand()) c%d = 0.0;\n   else c%d = 1.%d*rand();\n",
              i, i, i);
   if (!MULADD)
      for (i=0; i < lat; i++)
         fprintf(fpout,
      "   if (rand()|rand()|rand()) m%d = 0.0;\n   else m%d = 0.%d*rand();\n",
                 i, i, i+1);
   fprintf(fpout, "   t0 = time00();\n");
   fprintf(fpout, "   do\n   {\n");
   for (j=0; j < nama; j++)
   {
      for (i=0; i < lat; i++)
      {
         if (MULADD) fprintf(fpout, "      c%d += c%d * c%d;\n", i, i, i);
         else fprintf(fpout,"      c%d += m%d;\n      m%d = m%d * m%d;\n",
                      i, i, i, i, i);
      }
   }
   fprintf(fpout, "   }\n   while(--i);\n");
   fprintf(fpout, "   tim = time00() - t0;\n");
   fprintf(fpout, "   c0 = c0");
   if (!MULADD)
   {
      fprintf(fpout, "*m0");
      for (i=1; i < lat; i++) fprintf(fpout, " + c%d*m%d", i, i);
   }
   else
      for (i=1; i < lat; i++) fprintf(fpout, " + c%d", i);
   fprintf(fpout, ";\n");
   mf = (2.0 * nama)*lat ;
   fprintf(fpout, "   if (tim < 0.0) mf = tim = 0.0;\n");
   fprintf(fpout, "   else mf = (nreps*%lf) / (1000000.0 * tim);\n", mf);
   fprintf(fpout, "   if (PRINT) printf(\"%%.1f:   %s, lat=%d, time=%%f, mflop=%%f\\n\", (float) c0, tim, mf);\n", ma, lat);
   fprintf(fpout, "   else printf(\"      %%.0f: NFLOP=%%.0f, tim=%%f\\n\", (float) c0, nreps*%f, tim);\n", mf);
   fprintf(fpout, "   return(tim);\n");
   fprintf(fpout, "}\n");
/*
 * main()
 */
   fprintf(fpout, "\nmain(int nargs, char **args[])\n{\n");
   fprintf(fpout, "   long nreps = 16000000/%d;\n", (2*nama));
   fprintf(fpout, "   int i, k;\n");
   fprintf(fpout, "   double t0, tim, mf;\n");
   fprintf(fpout, "   FILE *fp;\n");
   fprintf(fpout, "   fp = fopen(\"res/%cmuladd%d_%d\", \"w\");\n", pre, MULADD, lat);
   fprintf(fpout, "   assert(fp != NULL);\n");

   fprintf(fpout, "   fprintf(stdout, \"Finding granularity of timer:\\n\");");
   fprintf(fpout, "   while(macase(nreps, 0) < 0.75) nreps *= 4;\n");
   fprintf(fpout, "   fprintf(stdout, \"Done.\\n\");");

   fprintf(fpout, "   for(k=0; k < 3; k++)\n   {\n");
   fprintf(fpout, "   tim = macase(nreps, 1);\n");
   fprintf(fpout, "   if (tim < 0.0) mf = tim = 0.0;\n");
   fprintf(fpout, "   else mf = (nreps*%lf) / (1000000.0 * tim);\n", mf);
   fprintf(fpout, "   if (fp) fprintf(fp, \"%%f\\n\", mf);\n");
   fprintf(fpout, "   }\n");
   fprintf(fpout, "   fclose(fp);\n");
   fprintf(fpout, "   exit(0);\n");
   fprintf(fpout, "}\n");
   fclose(fpout);
}

double macase(char *type, char pre, int MULADD, int lat, int nlma, int mflops)
{
   char fnam[80], ln[80];
   int i, j;
   double mflop[NTIM], t0, tavg;
   FILE *fp;

   lat++;
   sprintf(fnam, "res/%cmuladd%d_%d", pre, MULADD, lat);
   fp = fopen(fnam, "r");
   if (fp == NULL)
   {
      emit_muladd(type, pre, MULADD, lat, nlma, mflops);
      sprintf(ln, "make x%cmuladd\n", pre);
      assert(!system(ln));
      assert( (fp = fopen(fnam, "r")) != NULL );
   }
   t0 = 0.0;
   j = 0;
   for (i=0; i != NTIM; i++)
   {
      assert( fscanf(fp, "%lf", &mflop[i]) );
   }
   fclose(fp);
/*
 * Sort results, largest first
 */
   for (i=0; i != NTIM; i++)
   {
      for (j=i+1; j < NTIM; j++)
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
 * For walltime, return min time, else return median
 */
   #if defined(PentiumCPS) || defined(WALL)
      return(mflop[0]);
   #else
      return(mflop[NTIM/2]);
   #endif
}

int FindSafeNumRegs(char pre, char *type, int nlma, int mflops,
                    int muladd, int lat)
{
   int i, k, lat2, nreg, nrmax, nrmin;
   double mf, mf0;

   if (muladd) nreg = lat;
   else nreg = 2*lat;
   mf0 = macase(type, pre, muladd, lat-1, nlma, mflops);
   for (k=1; k <= nreg; k *= 2);
   k /= 2;

   do
   {
      k *= 2;
      if (muladd) lat2 = k;
      else { lat2 = k / 2; k = lat2 * 2; }

      mf = macase(type, pre, muladd, lat2-1, nlma, mflops);
      if (k > 256)
      {
         fprintf(stdout,
"      Either number of registers is > 256, or this test is not working.\n");
         fprintf(stdout, "      Unable to find NREG, setting to 32\n");
         return(32);
         break;
      }
      fprintf(stdout, "      nreg=%d, mflop = %.2lf (peak %.2lf)\n",k, mf, mf0);
   }
   while (mf*BETTER >= mf0);
   fprintf(stdout, "\n  nreg < %d (drop to %.2f%%)\n\n", k, mf/mf0);

   nrmax = k;
   nrmin = k / 2;
   do
   {
      k = nrmin + (nrmax - nrmin) / 2;
      if (muladd) lat2 = k;
      else lat2 = k / 2;
      mf = macase(type, pre, muladd, lat2-1, nlma, mflops);
      if (mf*BETTER >= mf0) nrmin = k;
      else nrmax = k;
      fprintf(stdout, "      nreg=%d, mflop = %.2lf (peak %.2lf)\n",k, mf, mf0);
   }
   while (nrmax-nrmin > 1);
   return(nrmin);
}

void RunTime(char pre, char *type, int nlma, int mflop)
/*
 * If necessary, fills in mflop field by timing previously chosen case
 */
{
   char fnam[64];
   int muladd, lat, nreg;
   double mf;
   FILE *fp;

   sprintf(fnam, "res/%cMULADD", pre);
   fp = fopen(fnam, "r");
   if (fp == NULL) return;
   fscanf(fp, "%d", &muladd);
   fscanf(fp, "%d", &lat);
   fscanf(fp, "%lf", &mf);
   fscanf(fp, "%d", &nreg);
   fclose(fp);
   if (mf < 0.0)
   {
      mf = macase(type, pre, muladd, lat-1, nlma, mflop);
      fp = fopen(fnam, "w");
      assert(fp != NULL);
      fprintf(fp, "%d\n", muladd);
      fprintf(fp, "%d\n", lat);
      fprintf(fp, "%lf\n", mf);
      fprintf(fp, "%d\n", nreg);
      fclose(fp);
   }
}

static int Mylcm(const int M, const int N)
/*
 * Returns least common multiple (LCM) of two positive integers M & N by
 * computing greatest common divisor (GCD) and using the property that
 * M*N = GCD*LCM.
 */
{
   register int tmp, max, min, gcd=0;

   if (M != N)
   {
      if (M > N) { max = M; min = N; }
      else { max = N; min = M; }
      if (min > 0)  /* undefined for negative numbers */
      {
         do  /* while (min) */
         {
            if ( !(min & 1) ) /* min is even */
            {
               if ( !(max & 1) ) /* max is also even */
               {
                  do
                  {
                     min >>= 1;
                     max >>= 1;
                     gcd++;
                     if (min & 1) goto MinIsOdd;
                  }
                  while ( !(max & 1) );
               }
               do min >>=1 ; while ( !(min & 1) );
            }
/*
 *          Once min is odd, halve max until it too is odd.  Then, use
 *          property that gcd(max, min) = gcd(max, (max-min)/2)
 *          for odd max & min
 */
MinIsOdd:
            if (min != 1)
            {
               do  /* while (max >= min */
               {
                  max -= (max & 1) ? min : 0;
                  max >>= 1;
               }
               while (max >= min);
            }
            else return( (M*N) / (1<<gcd) );
            tmp = max;
            max = min;
            min = tmp;
         }
         while(tmp);
      }
      return( (M*N) / (max<<gcd) );
   }
   else return(M);
}

main(int nargs, char *args[])
{
   char pre, *type, fnam[64];
   int k, kma, kno, kbest, maxlat, mflop, nlma=1;
   int MULADD=1;
   double mf0, mf, mfma=0, mfno=0;
   FILE *fp;

   if (nargs != 4 && nargs != 5)
   {
      fprintf(stderr, "usage: %s <pre> <maxlat> <mflop> [<nlma>]\n", args[0]);
      exit(-1);
   }
   pre = *args[1];
   if (pre == 'd') type = "double";
   else if (pre == 's') type = "float";
   else if (pre == 'q') type = "long double";
   else
   {
      fprintf(stderr, "\n\nUnknown prefix \'%c\'; aborting.\n\n", pre);
      exit(-1);
   }
   maxlat = atoi(args[2]);
   mflop = atoi(args[3]);  /* this is ignored now */
   if (nargs == 5) nlma = atoi(args[4]);
   else
   {
      nlma = 2;
      for (k=3; k < maxlat; k++)
         nlma = Mylcm(nlma, k);
      nlma = ((512+nlma-1)/nlma)*nlma;
   }
   if (maxlat <= 0)
   {
      RunTime(pre, type, nlma, mflop);
      exit(0);
   }

   kma = kno = kbest = -1;
   mf = mfma = mfno = 0.0;

   maxlat = 16;
   for (k=0; k < maxlat; k++)
   {
      mf0 = macase(type, pre, 1, k, nlma, mflop);
      if (mf0 > BETTER*mfma)
      {
         mfma = mf0;
         kma = k;
      }
      fprintf(stdout, "      MULADD=%d, lat=%d, mf=%.2lf\n", 1, k+1, mf0);
/*
 *    If we've timed the short lengths and longer is not improving perf, quit
 */
      if (k > 6 && mf0 < mfma)
         break;
   }
   for (k=0; k < maxlat; k++)
   {
      mf0 = macase(type, pre, 0, k, nlma, mflop);
      if (mf0 > BETTER*mfno)
      {
         mfno = mf0;
         kno = k;
      }
      fprintf(stdout, "      MULADD=%d, lat=%d, mf=%.2lf\n", 0, k+1, mf0);
/*
 *    If we've timed the short lengths and longer is not improving perf, quit
 */
      if (k > 6 && mf0 < mfma)
         break;
   }
/*
 * Award ties to separate multiply/add, since it may just be the compiler
 * making the machine behave like it has muladd
 */
   if (BETTER*mfno >= mfma)
   {
      MULADD = 0;
      mf = mfno;
      kbest = kno;
   }
   else
   {
      mf = mfma;
      kbest = kma;
   }
   k = FindSafeNumRegs(pre, type, nlma, mflop, MULADD, kbest+1);
   fprintf(stdout, "\nLower bound on number of registers = %d\n", k);
   sprintf(fnam, "res/%cMULADD", pre);
   fp = fopen(fnam, "w");
   assert(fp != NULL);
   fprintf(fp, "%d\n", MULADD);
   fprintf(fp, "%d\n", kbest+1);
   fprintf(fp, "%lf\n", mf);
   fprintf(fp, "%d\n", k);
   fclose(fp);
   fprintf(stdout,
"\n\nAccording to these timings, I think you need a %d cycle latency,\n", kbest+1);
   if (MULADD) fprintf(stdout, "with a combined multiply/add instruction.\n");
   else fprintf(stdout, "with separate multiply and add instructions.\n\n");
   fprintf(stdout, "These parameters gave a performance of %f MFLOPS.\n\n", mf);
   exit(0);
}

