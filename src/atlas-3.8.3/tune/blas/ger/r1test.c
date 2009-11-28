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
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_lvl2.h"
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

static void dumb_ger(int Conj, int M, int N, SCALAR alpha, TYPE *X, int incX,
                     TYPE *Y, int incY, TYPE *A, int lda)
{
   #ifdef TCPLX
      TYPE tmp[2];
   #endif
   int j;
   ATL_assert(incX == 1);
   ATL_assert( SCALAR_IS_ONE(alpha) );
   for (j=0; j < N; j++)
   {
   #ifdef TREAL
      Mjoin(PATL,axpy)(M, Y[j*incY], X, 1, A+lda*j, 1);
   #else
      tmp[0] = Y[2*j*incY];
      tmp[1] = Y[2*j*incY+1];
      if (Conj) tmp[1] = -tmp[1];
      Mjoin(PATL,axpy)(M, tmp, X, 1, A+2*lda*j, 1);
   #endif
   }
}

static int CheckAns(int M, int N, TYPE *G, int ldg, TYPE *U, int ldu)
{
   TYPE diff, eps;
   int i, j, ierr=0;
   #ifdef TREAL
      const int M2 = M, emul=4;
   #else
      const int M2 = M<<1, emul=16;
      ldg <<= 1; ldu <<= 1;
   #endif

   eps = Mjoin(PATL,epsilon)();
   for (j=0; j < N; j++, G += ldg, U += ldu)
   {
      for (i=0; i < M2; i++)
      {
         diff = G[i] - U[i];
         if (diff < ATL_rzero) diff = -diff;
         if (diff > emul*eps)
         {
            fprintf(stderr, "A(%d,%d): Good=%f, Computed=%f\n",
                    i, j, G[i], U[i]);
            if (!ierr) ierr = i+j*M+1;
         }
      }
   }
   return(ierr);
}

static int RunTest(int CONJ, int M, int N, int incY, int lda)
{
   #ifdef TCPLX
      TYPE one[2] = {ATL_rone, ATL_rzero};
   #else
      TYPE one = ATL_rone;
   #endif
   TYPE *A, *A0, *X, *Y, *y;
   const int aincY = Mabs(incY);
   int ierr;

   fprintf(stdout, "   TEST CONJ=%d, M=%d, N=%d, lda=%d, incY=%d, STARTED\n",
           CONJ, M, N, lda, incY);
   A = FA_malloc(ATL_MulBySize(lda)*N, FAa, MAa);
   A0 = FA_malloc(ATL_MulBySize(M)*N, FAa, MAa);
   X = FA_malloc(ATL_MulBySize(M), FAx, MAx);
   Y = FA_malloc(ATL_MulBySize(N)*aincY, FAy, MAy);
   ATL_assert(A && A0 && X && Y);
   Mjoin(PATL,gegen)(1, N, Y, aincY, N*aincY);
   Mjoin(PATL,gegen)(M, 1, X, M, N*aincY+127*50+77);
   Mjoin(PATL,gegen)(M, N, A0, M, N*M+513*7+90);
   Mjoin(PATL,gegen)(M, N, A, lda, N*M+513*7+90);
   if (incY < 0) Y += (N-1) * (aincY SHIFT);

#ifdef TCPLX
   if (CONJ)
      Mjoin(PATL,ger1c_a1_x1_yX)(M, N, one, X, 1, Y, incY, A, lda);
   else Mjoin(PATL,ger1u_a1_x1_yX)(M, N, one, X, 1, Y, incY, A, lda);
#else
   Mjoin(PATL,ger1_a1_x1_yX)(M, N, one, X, 1, Y, incY, A, lda);
#endif
   dumb_ger(CONJ, M, N, one, X, 1, Y, incY, A0, M);

   if (incY < 0) Y -= (N-1) * (aincY SHIFT);
   FA_free(Y, FAy, MAy);
   FA_free(X, FAx, MAx);
   ierr = CheckAns(M, N, A0, M, A, lda);
   FA_free(A, FAa, MAa);
   FA_free(A0, FAa, MAa);

   if (ierr)
      fprintf(stdout, "   TEST CONJ=%d, M=%d, N=%d, lda=%d, incY=%d, FAILED\n",
              CONJ, M, N, lda, incY);
   else fprintf(stdout,
                "   TEST CONJ=%d, M=%d, N=%d, lda=%d, incY=%d, PASSED\n",
                CONJ, M, N, lda, incY);
   return(ierr);
}

static int RunTests(int M, int N, int incY, int lda)
{
   int j, k=0, ierr = 0;

   ATL_assert(incY != 0);

#ifdef TCPLX
   for (k=0; k < 2; k++)
#endif
   {
      j = RunTest(k, M, N, 1, lda);
      if (!ierr) ierr = j;
      j = RunTest(k, M, N, incY, lda);
      if (!ierr) ierr = j;
      j = RunTest(k, M, N, -incY, lda);
      if (!ierr) ierr = j;
   }
   return(ierr);
}

static void PrintUsage(char *nam)
{
   fprintf(stderr, "USAGE: %s -m <M> -n <N> -Y <incY> -l <lda>\n",
           nam);
   exit(-1);
}

void GetFlags(int nargs, char **args, int *M, int *N, int *incY, int *lda)
{
   int i, k;
   char ch;

   *M = 997;
   *N = 177;
   *incY = 3;
   *lda = -1;

   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'F':
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
         else
         {
            if (k < 0)
              MAx = -k;
            else
              FAx = k;
         }
         break;
      case 'm':
         *M = atoi(args[++i]);
         break;
      case 'n':
         *N = atoi(args[++i]);
         break;
      case 'l':
         *lda = atoi(args[++i]);
         break;
      case 'Y':
         *incY = atoi(args[++i]);
         break;
      default:
         PrintUsage(args[0]);
      }
   }
   if (*lda < *M) *lda = *M + 9;
}

main(int nargs, char **args)
{
   int M, N, incY, lda, ierr=0;

   GetFlags(nargs, args, &M, &N, &incY, &lda);
   ierr = RunTests(M, N, incY, lda);
   exit(ierr);
}
