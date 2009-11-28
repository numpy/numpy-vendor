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

#include "atlas_misc.h"
#include "atlas_lvl2.h"
#include "atlas_fopen.h"

int FAx=0, MAx=0, FAy=0, MAy=0, FAa=0, MAa=0;
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

struct FA_allocs
{
   void *mem, *memA;
   struct FA_allocs *next;
} *allocQ=NULL;

struct FA_allocs *NewAlloc(size_t size, struct FA_allocs *next,
                           int align, int misalign)
/*
 * Allocates size allocation that is aligned to [align], but not aligned
 * to [misalign].  Therefore, misalign > align.  Align must minimally be sizeof
 * while misalign may be 0 if we don't need to avoid a particular alignment.
 */
{
   void *vp;
   char *cp;
   struct FA_allocs *ap;
   int n, i;
   const int malign = align >= misalign ? align : misalign;

   n = size + align + malign;
   i = (n >> 3)<<3;
   if (n != i)
      n += n - i;
   cp = malloc(n + sizeof(struct FA_allocs));
   assert(cp);
   ap = (struct FA_allocs *) (cp + n);
   ap->mem = cp;
/*
 * Align to min alignment
 */
   ap->memA = align ? (void*) ((((size_t) cp)/align)*align + align) : cp;
/*
 * Misalign to misalign
 */
   if (misalign)
   {
      if (((size_t)ap->memA)%misalign == 0)
         ap->memA = ((char*)ap->memA) + align;
   }
   ap->next = next;
   return(ap);
}

/*
 * no-align malloc free retaining system default behavior
 */
void *NA_malloc(size_t size)
{
   return(malloc(size));
}
void *NA_calloc(size_t n, size_t size)
{
   return(calloc(n, size));
}
void NA_free(void *ptr)
{
   free(ptr);
}


/*
 * malloc/free pair that aligns data to align, but not to misalign
 */
void *FA_malloc(size_t size, int align, int misalign)
{
   if ((!misalign && align <= 8) || !size)
      return(malloc(size));
   else
   {
      allocQ = NewAlloc(size, allocQ, align, misalign);
      return(allocQ->memA);
   }
}
void *FA_calloc(size_t n, size_t size, int align, int misalign)
{
   char *cp;
   int *ip;
   double *dp;
   size_t i;
   size_t tsize;
   tsize = n * size;
   cp = FA_malloc(tsize, align, misalign);
   if (size == sizeof(int))
      for (ip=(int*)cp,i=0; i < n; i++)
        ip[i] = 0;
   else if (size == sizeof(double))
      for (dp=(double*)cp,i=0; i < n; i++)
        dp[i] = 0.0;
   else
      for (i=0; i < tsize; i++)
        cp[i] = 0;
   return(cp);
}

void FA_free(void *ptr, int align, int misalign)
/*
 * Part of malloc/free pair that aligns data to FALIGN
 */
{
   struct FA_allocs *ap, *prev;
   if (ptr)
   {
      if ((!misalign && align <= 8))
         free(ptr);
      else
      {
         for (ap=allocQ; ap && ap->memA != ptr; ap = ap->next) prev = ap;
         if (!ap)
         {
            fprintf(stderr, "Couldn't find mem=%ld\nmemQ=\n", ptr);
            for (ap=allocQ; ap; ap = ap->next)
               fprintf(stderr, "   %ld, %ld\n", ap->memA, ap->mem);
         }
         assert(ap);
         if (ap == allocQ)
            allocQ = allocQ->next;
         else
            prev->next = ap->next;
         free(ap->mem);
      }
   }
}

#define ATL_NoBlock(iflag_) ( ((iflag_) | 32) == (iflag_) )

double time00();
#ifdef TREAL
   #define test_ger(M, N, alpha, X, incX, Y, incY, A, lda) \
      Mjoin(PATL,ger)(M, N, alpha, X, incX, Y, incY, A, lda)
#else
   #define test_ger(M, N, alpha, X, incX, Y, incY, A, lda) \
      Mjoin(PATL,geru)(M, N, alpha, X, incX, Y, incY, A, lda)
#endif

double gercase(const int MFLOP, const int M, const int N, const SCALAR alpha,
               const int lda)
{
   unsigned long reps;
   const int aincX = 1, aincY = 1;
   const int incX = 1, incY = 1;
   int i, lx, ly, la;
   #ifdef TREAL
      const double flops = 2.0 * M * N;
   #else
      const double flops = 8.0 * M * N;
   #endif
   double ttest, mftest, t0;
   const int inca = lda*N SHIFT, incx = M*incX SHIFT, incy = N*incY SHIFT;
   TYPE *a, *A, *stA, *A0, *x, *X, *X0, *stX, *y, *Y, *Y0, *stY;
   #ifdef TREAL
      const TYPE nalpha = -alpha;
      TYPE alp = alpha;
   #else
      const TYPE *alp = alpha;
      TYPE nalpha[2];
   #endif
   #ifdef TCPLX
      nalpha[0] = -alpha[0]; nalpha[1] = alpha[1];
   #endif

   i = (ATL_DivBySize(L2SIZE)+M-1)/M;
   lx = i * M * aincX;
   X0 = X = x = FA_malloc(ATL_MulBySize(lx), FAx, MAx);
   if (x == NULL) return(-1);

   i = (ATL_DivBySize(L2SIZE)+N-1)/N;
   ly = i * N * aincY;
   Y0 = Y = y = FA_malloc(ATL_MulBySize(ly), FAy, MAy);
   if (y == NULL)
   {
      FA_free(x, FAx, MAx);
      return(-1);
   }
   i = (ATL_DivBySize(L2SIZE)+M*N)/(M*N);
   la = i * lda*N;
   A0 = A = a = FA_malloc(ATL_MulBySize(la), FAa, MAa);
   if (a == NULL)
   {
      FA_free(x, FAy, MAy);
      FA_free(y, FAy, MAy);
      return(-1);
   }
   if (incX < 1)
   {
      stX = x;
      x = X = x + (lx SHIFT);
   }
   else stX = x + (lx SHIFT);
   if (incY < 1)
   {
      stY = y;
      y = Y = y + (ly SHIFT);
   }
   else stY = y + (ly SHIFT);
   stA = a + (la SHIFT);

   reps = (MFLOP * 1000000.0) / flops;
   if (reps < 1) reps = 1;
   Mjoin(PATL,gegen)(ly, 1, Y0, ly, M*incY);
   Mjoin(PATL,gegen)(lx, 1, X0, lx, N*incY+127*50+77);
   Mjoin(PATL,gegen)(la, 1, A0, la, N*M+513*7+90);

   t0 = time00();
   for (i=reps; i; i--)
   {
      test_ger(M, N, alp, x, incX, y, incY, A0, lda);
      x += incx;
      y += incy;
      a += inca;
      if (x == stX) x = X;
      if (y == stY) y = Y;
      if (a == stA)
      {
         a = A;
         if (alp == alpha) alp = nalpha;
         else alp = alpha;
      }
   }
   ttest = time00() - t0;

   if (ttest > 0.0) mftest = (reps * flops) / (1000000.0 * ttest);
   else mftest = 0.0;

   FA_free(A0, FAa, MAa);
   FA_free(X0, FAx, MAx);
   FA_free(Y0, FAy, MAy);
   return(mftest);
}

void PrintUsage(char *nam)
{
   fprintf(stderr, "USAGE: %s -C <case #> -l <l1mul> -F <mflop> -m <M> -n <N> -f <iflag> -a <alpha> -o <outfile>\n\n", nam);
   exit(-1);
}

void GetFlags(int nargs, char **args, char *pre, int *MFLOP,
              int *M, int *N, TYPE *alpha, int *lda, char *outnam)
{
   int l1mul;
   int i, k, cas, iflag=0;
   char ch;

   l1mul = 75;
   outnam[0] = '\0';
   #if defined(DREAL)
      *pre = 'd';
   #elif defined(SREAL)
      *pre = 's';
   #elif defined(SCPLX)
      *pre = 'c';
   #elif defined(DCPLX)
      *pre = 'z';
   #endif
   #ifdef ATL_nkflop
      *MFLOP = (ATL_nkflop) / 2000.0;
      if (*MFLOP < 1) *MFLOP = 1;
   #else
      *MFLOP = 75;
   #endif
   *lda = *M = *N = 1000;
   *alpha = ATL_rone;
   #ifdef TCPLX
      alpha[1] = ATL_rzero;
   #endif
   cas = 1;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'C': /* case */
         cas = atoi(args[++i]);
         break;
      case 'F': /* mflops */
         ch = args[i][2];
         k = atoi(args[++i]);
         if (ch == 'a')
         {
            if (k < 0)
              MAa = -k;
            else
              FAa = k;
         }
         else if (ch == 'y')
         {
            if (k < 0)
              MAy = -k;
            else
              FAy = k;
         }
         else if (ch == 'x')
         {
            if (k < 0)
              MAx = -k;
            else
              FAx = k;
         }
         else *MFLOP = k;
         break;
      case 'f': /* iflag */
         iflag = atoi(args[++i]);
         break;
      case 'm':
         *M = atoi(args[++i]);
         break;
      case 'n':
         *N = atoi(args[++i]);
         break;
      case 'l':
         l1mul = atoi(args[++i]);
         break;
      case 'a':
         *alpha = atof(args[++i]);
         #ifdef TCPLX
            alpha[1] = atof(args[++i]);
         #endif
         break;
      case 'o':
         strcpy(outnam, args[++i]);
         break;
      default:
         PrintUsage(args[0]);
      }
   }
   if (outnam[0] == '\0')
   {
      if (ATL_NoBlock(iflag))
         sprintf(outnam, "res/%cger1_%d_0", *pre, cas);
      else sprintf(outnam, "res/%cger1_%d_%d", *pre, cas, l1mul);
   }
}

main(int nargs, char **args)
{
   char pre, fnam[128];
   int MFLOP, M, N, lda, cas, i;
   double mf, mfs[3];
   #ifdef TREAL
      TYPE alpha;
   #else
      TYPE alpha[2];
   #endif
   FILE *fp;

   GetFlags(nargs, args, &pre, &MFLOP, &M, &N, SADD alpha, &lda, fnam);

   if (!FileExists(fnam))
   {
      fp = fopen(fnam, "w");
      ATL_assert(fp);
      for (i=0; i < 3; i++)
      {
         mf = gercase(MFLOP, M, N, alpha, lda);
         fprintf(stdout, "      %s : %f MFLOPS\n", fnam, mf);
         fprintf(fp, "%lf\n", mf);
         mfs[i] = mf;
      }
   }
   else
   {
      fp = fopen(fnam, "r");
      for (i=0; i < 3; i++) ATL_assert(fscanf(fp, " %lf", &mfs[i]) == 1);
   }
   fclose(fp);
   mf = (mfs[0] + mfs[1] + mfs[2]) / 3.0;
   fprintf(stdout, "   %s : %.2f MFLOPS\n", fnam, mf);
   exit(0);
}
