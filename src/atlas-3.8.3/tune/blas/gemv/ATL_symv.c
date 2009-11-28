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
#include <assert.h>
#include Mstr(Mjoin(Mjoin(atlas_,PRE),mv.h))

#define gemvN Mjoin(PATL,gemvN_a1_x1_b1_y1)
#define gemvT Mjoin(PATL,gemvT_a1_x1_b1_y1)

static void refsymvL
   (const int N, const SCALAR alpha, const TYPE *A, const int lda,
    const TYPE *X, const int incX, const SCALAR beta, TYPE *Y, const int incY)
{
#ifdef TREAL
   register TYPE t0, t1, t2;
   const int ldap1=lda+1;
   int j, i, iaij, jaj;

   for (j=jaj=0; j < N; j++, jaj+=ldap1)
   {
      t0 = X[j]; t1 = ATL_rzero; Y[j] += t0 * A[jaj];
      for (i=j+1, iaij=1+jaj; i < N; i++, iaij++)
      {
         t2 = A[iaij];
         Y[i] += t2 * t0;
         t1 += t2 * X[i];
      }
      Y[j] += t1;
   }
#else
   register TYPE yr, yi, ar, ai, xr, xi, Xr, Xi;
   const int ldap1=(lda+1)SHIFT, N2=(N SHIFT);
   const TYPE *a, *x;
   TYPE *y=Y;
   int i, j;

   for (j=0; j < N; j++, A += ldap1, X += 2, Y += 2)
   {
      yr = yi = ATL_rzero;
      xr = *X; xi = X[1];
      ar = *A;
      *Y += ar*xr; Y[1] += ar*xi;
      for (y=Y+2, x=X+2, a=A+2, i=j+1; i < N; i++, a += 2, x += 2, y += 2)
      {
         ar = *a; ai = a[1];
         Xr = *x; Xi = x[1];
         *y   += ar*xr + ai*xi;
         y[1] += ar*xi - ai*xr;
         yr += ar * Xr - ai * Xi;
         yi += ar * Xi + ai * Xr;
      }
      *Y = yr; Y[1] = yi;
   }
#endif
}
static void RsymvL(const int N, const TYPE *A, const int lda, const TYPE *X,
                   const SCALAR beta, TYPE *Y)
{
   int mb, nb, jb, jbs, j, m;
   const TYPE *x0=X, *A0=A;
   const SCALAR alpha=beta;
   TYPE *y0=Y;

   ATL_GetPartSYMV(A, lda, &mb, &nb);

   for (j=0; j < N; j += nb)
   {
      jb = N-j;
      jb = Mmin(jb,nb);
      refsymvL(jb, alpha, A, lda, X, 1, beta, Y, 1);
      m = N-j-jb;
      if (m)
      {
         jbs = jb SHIFT; X += jbs; Y += jbs;
         gemvT(jb, m, alpha, A+jbs, lda, X, 1, beta, y0, 1);
         gemvN(m, jb, alpha, A+jbs, lda, x0, 1, beta, Y, 1);
         A += jbs*(lda+1); x0 = X; y0 = Y;
      }
   }
}
void Mjoin(PATL,symv)
   (const enum ATLAS_UPLO Uplo, const int N, const SCALAR alpha,
    const TYPE *A, const int lda, const TYPE *X, const int incX,
    const SCALAR beta, TYPE *Y, const int incY)
{
   int mb, nb, jb, mb1, incA1, incA, incXY, incXY1, n, j;
   const int lda2=(lda SHIFT);
   const TYPE *x0=X, *x1, *A0=A, *A1;
   TYPE *y1, *y0=Y;

   assert(incX==1 && incY==1 && Uplo == AtlasLower);
   #ifdef TREAL
      assert(alpha == ATL_rone && beta == ATL_rone);
   #else
      assert(*alpha == ATL_rone && *beta == ATL_rone);
      assert(alpha[1] == ATL_rzero && beta[1] == ATL_rzero);
   #endif

   ATL_GetPartSYMV(A, lda, &mb, &nb);
   mb1 = N - ( (N-1) / mb ) * mb;
   incA1 = nb * lda2;  incXY1 = (nb SHIFT);
   incA = incXY = mb SHIFT;
   n = (N-mb)SHIFT;
   A += n; X += n; Y += n;
   for (n=N-mb; n > 0; n -= mb, A -= incA, X -= incXY, Y -= incXY)
   {
      RsymvL(mb, A+n*lda2, lda, X, beta, Y);
      for (j=0, A1=A, x1=x0, y1=y0; j < n; j += nb, A1 += incA1, x1 += incXY1,
           y1 += incXY1)
      {
         jb = n - j;
         jb = Mmin(jb, nb);
         gemvT(jb, mb, alpha, A1, lda, X, 1, beta, y1, 1);
         gemvN(mb, jb, alpha, A1, lda, x1, 1, beta, Y, 1);
      }
   }
   RsymvL(mb1, A0, lda, x0, beta, y0);
}
