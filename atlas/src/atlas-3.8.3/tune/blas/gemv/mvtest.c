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

int CheckAns(int M, int N, SCALAR beta, TYPE *Y, TYPE *Y0)
{
   TYPE eps, diff;
   #ifdef TCPLX
      const int M2 = M<<1, mul=8;
   #else
      const int M2 = M, mul=4;
   #endif
   int i, ierr=0;

   eps = Mjoin(PATL,epsilon)();
   for (i=0; i < M2; i++)
   {
      diff = Y0[i] - Y[i];
      if (diff < 0.0) diff = -diff;
      if (diff > N*mul*eps + 2.0*eps)
      {
         ierr = i+1;
         fprintf(stderr, "Y[%d] = %f, expected %f\n", i, Y[i], Y0[i]);
      }
   }
   return(ierr);
}

#ifdef TREAL
static void dumb_gemv(enum ATLAS_TRANS TA, int M, int N, SCALAR alpha,
                      TYPE *A, int lda, TYPE *X, SCALAR beta, TYPE *Y)
{
   register TYPE y0;
   int i, j;

   ATL_assert(alpha == ATL_rone);
   if (beta == ATL_rzero) Mjoin(PATL,zero)(M, Y, 1);
   else if (beta != ATL_rone) Mjoin(PATL,scal)(M, beta, Y, 1);

   if (TA == AtlasNoTrans)
   {
      for (j=0; j < N; j++) Mjoin(PATL,axpy)(M, X[j], A+j*lda, 1, Y, 1);
   }
   else
   {
      for (i=0; i < M; i++)
      {
         y0 = Y[i];
         for (j=0; j < N; j++) y0 += A[i*lda+j] * X[j];
         Y[i] = y0;
      }
   }
}
#else
static void dumb_gemv(enum ATLAS_TRANS TA, int M, int N, SCALAR alpha,
                      TYPE *A, int lda, TYPE *X, SCALAR beta, TYPE *Y)
{
   register TYPE ry, iy, rx, ix, ra, ia;
   int i, j;
   const int M2 = M<<1, N2 = N<<1;

   ATL_assert(SCALAR_IS_ONE(alpha));
   if ( SCALAR_IS_ZERO(beta) ) Mjoin(PATL,zero)(M, Y, 1);
   else if ( !SCALAR_IS_ONE(beta) ) Mjoin(PATL,scal)(M, beta, Y, 1);

   if (TA == AtlasNoTrans)
      for (j=0; j < N; j++) Mjoin(PATL,axpy)(M, X+2*j, A+j*2*lda, 1, Y, 1);
   else if (TA == AtlasConj)
      for (j=0; j < N; j++) Mjoin(PATL,axpyConj)(M, X+2*j, A+j*2*lda, 1, Y, 1);
   else if (TA == AtlasTrans)
   {
      for (i=0; i < M2; i += 2)
      {
         ry = Y[i];
         iy = Y[i+1];
         for (j=0; j < N2; j += 2)
         {
            ra = A[i*lda+j]; ia = A[i*lda+j+1];
            rx = X[j];       ix = X[j+1];
            ry += rx * ra - ix * ia;
            iy += rx * ia + ix * ra;
         }
         Y[i] = ry;
         Y[i+1] = iy;
      }
   }
   else if (TA == AtlasConjTrans)
   {
      for (i=0; i < M2; i += 2)
      {
         ry = Y[i];
         iy = Y[i+1];
         for (j=0; j < N2; j += 2)
         {
            ra = A[i*lda+j]; ia = -A[i*lda+j+1];
            rx = X[j];       ix = X[j+1];
            ry += rx * ra - ix * ia;
            iy += rx * ia + ix * ra;
         }
         Y[i] = ry;
         Y[i+1] = iy;
      }
   }
}
#endif

int RunTest(enum ATLAS_TRANS TA, int M, int N, int lda, TYPE *beta0)
{
   #ifdef TCPLX
      TYPE alpha[2] = {ATL_rone, ATL_rzero};
      #define beta beta0
   #else
      TYPE alpha= ATL_rone;
      TYPE beta = *beta0;
   #endif
   TYPE *A, *X, *Y, *Y0;
   int mA, nA, ierr;
   char cta, ln[32];

   if (TA == AtlasNoTrans) cta = 'N';
   else if (TA == AtlasTrans) cta = 'T';
   else if (TA == AtlasConjTrans) cta = 'C';
   else cta = '-';
   if (TA == AtlasNoTrans || TA == AtlasConj) { mA = M; nA = N; }
   else { mA = N; nA = M; }

   #ifdef TREAL
      fprintf(stdout, "   TEST TA=%c, M=%d, N=%d, lda=%d, beta=%f STARTED\n",
              cta, M, N, lda, beta);
   #else
      fprintf(stdout,
              "   TEST TA=%c, M=%d, N=%d, lda=%d, beta=(%f,%f) STARTED\n",
              cta, M, N, lda, *beta, beta[1]);
   #endif
   A = FA_malloc(lda * ATL_MulBySize(nA), FAa, MAa);
   X = FA_malloc(ATL_MulBySize(N), FAx, MAx);
   Y = FA_malloc(ATL_MulBySize(M), FAy, MAy);
   Y0 = FA_malloc(ATL_MulBySize(M), FAy, MAy);
   ATL_assert(A && X && Y && Y0);
   Mjoin(PATL,gegen)(mA, nA, A, lda, M*N+lda);
   Mjoin(PATL,gegen)(N, 1, X, N, N*1139+52);
   if ( !SCALAR_IS_ZERO(beta) )
   {
      Mjoin(PATL,gegen)(M, 1, Y, M, M*8897+112);
      Mjoin(PATL,gegen)(M, 1, Y0, M, M*8897+112);
   }

   if (TA == AtlasNoTrans)
   {
      if ( SCALAR_IS_ZERO(beta) )
         Mjoin(PATL,gemvN_a1_x1_b0_y1)(M, N, alpha, A, lda, X, 1,
                                       beta, Y, 1);
      else if ( SCALAR_IS_ONE(beta) )
         Mjoin(PATL,gemvN_a1_x1_b1_y1)(M, N, alpha, A, lda, X, 1,
                                       beta, Y, 1);
      #ifdef TCPLX
         else if (beta[1] == ATL_rzero)
            Mjoin(PATL,gemvN_a1_x1_bXi0_y1)(M, N, alpha, A, lda, X, 1,
                                            beta, Y, 1);
      #endif
      else
         Mjoin(PATL,gemvN_a1_x1_bX_y1)(M, N, alpha, A, lda, X, 1,
                                       beta, Y, 1);
   }
   else if (TA == AtlasTrans)
   {
      if ( SCALAR_IS_ZERO(beta) )
         Mjoin(PATL,gemvT_a1_x1_b0_y1)(M, N, alpha, A, lda, X, 1,
                                       beta, Y, 1);
      else if ( SCALAR_IS_ONE(beta) )
         Mjoin(PATL,gemvT_a1_x1_b1_y1)(M, N, alpha, A, lda, X, 1,
                                       beta, Y, 1);
      #ifdef TCPLX
         else if (beta[1] == ATL_rzero)
            Mjoin(PATL,gemvT_a1_x1_bXi0_y1)(M, N, alpha, A, lda, X, 1,
                                            beta, Y, 1);
      #endif
      else
         Mjoin(PATL,gemvT_a1_x1_bX_y1)(M, N, alpha, A, lda, X, 1,
                                       beta, Y, 1);
   }
   #ifdef TCPLX
      else if (TA == AtlasConj)
      {
         if (beta[1] == ATL_rzero)
         {
            if (*beta == ATL_rzero)
               Mjoin(PATL,gemvNc_a1_x1_b0_y1)(M, N, alpha, A, lda, X, 1,
                                              beta, Y, 1);
            else if (*beta == ATL_rone)
               Mjoin(PATL,gemvNc_a1_x1_b1_y1)(M, N, alpha, A, lda, X, 1,
                                              beta, Y, 1);
            else
               Mjoin(PATL,gemvNc_a1_x1_bXi0_y1)(M, N, alpha, A, lda, X, 1,
                                                beta, Y, 1);
         }
         else Mjoin(PATL,gemvNc_a1_x1_bX_y1)(M, N, alpha, A, lda, X, 1,
                                             beta, Y, 1);
      }
      else if (TA == AtlasConjTrans)
      {
         if (beta[1] == ATL_rzero)
         {
            if (*beta == ATL_rzero)
               Mjoin(PATL,gemvC_a1_x1_b0_y1)(M, N, alpha, A, lda, X, 1,
                                             beta, Y, 1);
            else if (*beta == ATL_rone)
               Mjoin(PATL,gemvC_a1_x1_b1_y1)(M, N, alpha, A, lda, X, 1,
                                             beta, Y, 1);
            else
               Mjoin(PATL,gemvC_a1_x1_bXi0_y1)(M, N, alpha, A, lda, X, 1,
                                               beta, Y, 1);
         }
         else Mjoin(PATL,gemvC_a1_x1_bX_y1)(M, N, alpha, A, lda, X, 1,
                                            beta, Y, 1);
      }
   #endif
   dumb_gemv(TA, M, N, alpha, A, lda, X, beta, Y0);

   FA_free(A, FAa, MAa);
   FA_free(X, FAx, MAx);

   ierr = CheckAns(M, N, beta, Y, Y0);

   FA_free(Y, FAy, MAy);
   FA_free(Y0, FAy, MAy);

   if (ierr) sprintf(ln, "FAILED");
   else sprintf(ln, "PASSED");
   #ifdef TREAL
      fprintf(stdout, "   TEST TA=%c, M=%d, N=%d, lda=%d, beta=%f %s\n",
              cta, M, N, lda, beta, ln);
   #else
      fprintf(stdout, "   TEST TA=%c, M=%d, N=%d, lda=%d, beta=(%f,%f) %s\n",
              cta, M, N, lda, *beta, beta[1], ln);
   #endif
   return(ierr);
}

int RunTests(enum ATLAS_TRANS TA, int M, int N, int lda, SCALAR beta)
{
   int ierr = 0, i;
   #ifdef TREAL
      int nbeta = 3, incb = 1;
      TYPE betas[3];
   #else
      int nbeta = 8, incb = 2;
      TYPE betas[8];
   #endif

   #ifdef TREAL
      betas[0] = ATL_rzero;
      betas[1] = ATL_rone;
      betas[2] = beta;
   #else
      betas[0] = betas[1] = betas[3] = betas[5] = ATL_rzero;
      betas[2] = ATL_rone;
      betas[4] = betas[6] = *beta;
      betas[7] = beta[1];
   #endif
   for (i=0; i < nbeta; i += incb)
   {
      if (TA == AtlasNoTrans)
      {
         ierr = RunTest(AtlasNoTrans, M, N, lda, &betas[i]);
         if (ierr) break;
         #ifdef TCPLX
            ierr = RunTest(AtlasConj, M, N, lda, &betas[i]);
            if (ierr) break;
         #endif
      }
      else
      {
         ierr = RunTest(AtlasTrans, M, N, lda, &betas[i]);
         if (ierr) break;
         #ifdef TCPLX
            ierr = RunTest(AtlasConjTrans, M, N, lda, &betas[i]);
            if (ierr) break;
         #endif
      }
   }
   return(ierr);
}

void PrintUsage(char *nam)
{
   fprintf(stderr, "usage: %s -A t/n -m <M> -n <N> -l <lda> -b <betaX>\n", nam);
   exit(-1);
}

void GetFlags(int nargs, char **args, enum ATLAS_TRANS *TA,
              int *M0, int *MN, int *incM, int *N0, int *NN, int *incN,
              int *lda, TYPE *beta)
{
   int i, k;
   char ch;

   *TA = AtlasNoTrans;
   *M0 = *MN = *incM = 997;
   *N0 = *NN = *incN = 177;
   *lda = 0;
   *beta = 0.8;
   #ifdef TCPLX
      beta[1] = 0.3;
   #endif
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
         else if (args[i][2] == 'y')
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
      case 'A':
         i++;
         if (args[i][0] == 'N' || args[i][0] == 'n') *TA = AtlasNoTrans;
         else if (args[i][0] == 'T' || args[i][0] == 't') *TA = AtlasTrans;
         else if (args[i][0] == 'C' || args[i][0] == 'c') *TA = AtlasConjTrans;
         else PrintUsage(args[0]);
         break;
      case 'm':
         *M0 = *MN = atoi(args[++i]);
         break;
      case 'n':
         *N0 = *NN = atoi(args[++i]);
         break;
      case 'M':
         *M0 = atoi(args[++i]);
         *MN = atoi(args[++i]);
         *incM = atoi(args[++i]);
         break;
      case 'N':
         *N0 = atoi(args[++i]);
         *NN = atoi(args[++i]);
         *incN = atoi(args[++i]);
         break;
      case 'l':
         *lda = atoi(args[++i]);
         break;
      case 'b':
         *beta = atof(args[++i]);
         #ifdef TCPLX
            beta[1] = atof(args[++i]);
         #endif
         break;
      default:
         PrintUsage(args[0]);
      }
   }
   if (*lda == 0)
   {
      if (*TA == AtlasNoTrans  || *TA == AtlasConj) *lda = -(*M0 + 7);
      else *lda = -(*N0 + 7);
   }
}

main(int nargs, char **args)
/*
 * Test the GEMV primitives
 */
{
   #ifdef TREAL
      TYPE beta;
      #define bet &beta
   #else
      TYPE beta[2];
      #define bet beta
   #endif
   enum ATLAS_TRANS TA;
   int M0, MN, incM, N0, NN, incN, m, n, lda, ldagap=0, ierr=0;

   GetFlags(nargs, args, &TA, &M0, &MN, &incM,  &N0, &NN, &incN, &lda, bet);
   if (lda < 0)
   {
      lda = -lda;
      if (TA == AtlasNoTrans || TA == AtlasConj) ldagap = lda - M0;
      else ldagap = lda - N0;
      lda -= ldagap;
   }
   fprintf(stdout, "\n\nBEGINNING GEMV PRIMITIVE TESTING:\n\n");
   for (n=N0; n <= NN; n += incN)
   {
      for (m=M0; m <= MN; m += incM)
      {
         if (ldagap)
         {
            if (TA == AtlasNoTrans || TA == AtlasConj) lda = m + ldagap;
            else lda = n + ldagap;
         }
         ierr += RunTests(TA, m, n, lda, beta);
      }
   }
/*   ierr = RunTests(TA, M, N, lda, beta); */
   if (ierr == 0) fprintf(stdout, "\n\nGEMV PRIMITIVE PASSED ALL TESTS\n\n");
   else fprintf(stdout, "\n\nGEMV PRIMITIVE FAILED SOME TESTS\n\n");
   return(ierr);
}
