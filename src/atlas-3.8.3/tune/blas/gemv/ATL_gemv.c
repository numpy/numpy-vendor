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
#include "atlas_mv.h"

/*
 * Special case for NoTrans gemv based on multiple axpys instead of usual
 * multiple dot products
 */
#ifdef ATL_AXPYMV
static void gemvNaxpy
   (const enum ATLAS_TRANS TA, const int M, const int N, const SCALAR alpha,
    const TYPE *A, const int lda, const TYPE *X, const int incX,
    const SCALAR beta, TYPE *Y, const int incY)
{
   int i, mb, incYY, incy;
   void *vx=NULL, *vy=NULL;
   TYPE *x = (TYPE *) X, *y=Y;
   #ifdef TREAL
      #define ONE ATL_rone
      #define ZERO ATL_rzero
      TYPE alphaY = alpha, betaY, betamv=beta;
   #else
      const TYPE ONE[2]={ATL_rone, ATL_rzero}, ZERO[2]={ATL_rzero, ATL_rzero};
      const TYPE *alphaY = alpha, *betaY, *betamv=beta;
   #endif
   void (*putY)(const int N, const SCALAR alpha, const TYPE *X,
                const int incX, const SCALAR beta, TYPE *Y, const int incY);
   void (*gemv)(const int M, const int N, const SCALAR alpha,
                const TYPE *A, const int lda, const TYPE *X, const int incX,
                const SCALAR beta, TYPE *Y, const int incY);

   #ifdef ATL_mvNNOBLOCK
      mb = M;
   #else
      ATL_GetPartMVN(A, lda, &mb, &i);
      if (mb > M) mb = M;
   #endif
   if ( (incX != 1) || (!SCALAR_IS_ONE(alpha) && incY == 1 && M >= N) )
   {
      alphaY = ONE;
      vx = malloc(ATL_Cachelen + ATL_MulBySize(N));
      ATL_assert(vx);
      x = ATL_AlignPtr(vx);
      Mjoin(PATL,cpsc)(N, alpha, X, incX, x, 1);
   }
   if ( (incY != 1) || !SCALAR_IS_ONE(alphaY) )
   {
      betaY = beta;
      betamv = ZERO;
      #ifdef TREAL
         gemv = Mjoin(PATL,gemvN_a1_x1_b0_y1);
      #else
         if (TA == AtlasNoTrans) gemv = Mjoin(PATL,gemvN_a1_x1_b0_y1);
         else gemv = Mjoin(PATL,gemvNc_a1_x1_b0_y1);
      #endif
      vy = malloc(ATL_Cachelen + ATL_MulBySize(mb));
      ATL_assert(vy);
      y = ATL_AlignPtr(vy);
      putY = Mjoin(PATL,axpby);
   }
   else if ( SCALAR_IS_ONE(beta) )
   {
      putY = NULL;
      #ifdef TREAL
         gemv = Mjoin(PATL,gemvN_a1_x1_b1_y1);
      #else
         if (TA == AtlasNoTrans) gemv = Mjoin(PATL,gemvN_a1_x1_b1_y1);
         else gemv = Mjoin(PATL,gemvNc_a1_x1_b1_y1);
      #endif
   }
   else if ( SCALAR_IS_ZERO(beta) )
   {
      putY = NULL;
      #ifdef TREAL
         gemv = Mjoin(PATL,gemvN_a1_x1_b0_y1);
      #else
         if (TA == AtlasNoTrans) gemv = Mjoin(PATL,gemvN_a1_x1_b0_y1);
         else gemv = Mjoin(PATL,gemvNc_a1_x1_b0_y1);
      #endif
   }
   else
   {
      putY = NULL;
      #ifdef TREAL
         gemv = Mjoin(PATL,gemvN_a1_x1_bX_y1);
      #else
         if (TA == AtlasNoTrans) gemv = Mjoin(PATL,gemvN_a1_x1_bX_y1);
         else gemv = Mjoin(PATL,gemvNc_a1_x1_bX_y1);
      #endif
   }
   i = M;
   incYY = mb * incY SHIFT;
   if (putY) incy = 0;
   else incy = incYY;

   do
   {
      if (mb > i) mb = i;
      gemv(mb, N, ONE, A, lda, x, 1, betamv, y, 1);
      if (putY) putY(mb, alphaY, y, 1, betaY, Y, incY);
      A += mb SHIFT;
      Y += incYY;
      y += incy;
      i -= mb;
   }
   while(i);

   if (vx) free(vx);
   if (vy) free(vy);
}
#endif

void Mjoin(PATL,gemv)
   (const enum ATLAS_TRANS TA, const int M, const int N, const SCALAR alpha,
    const TYPE *A, const int lda, const TYPE *X, const int incX,
    const SCALAR beta, TYPE *Y, const int incY)
/*
 *  y <- alpha * A * x + beta*y   .. or ..  y <- alpha * A' * x + beta*y
 *  .. or ..  y <- alpha * conj(A') * x + beta*y
 *
 *  Note:  M = length of Y
 *         N = length of X
 */
{
   const int BetaIsZero = SCALAR_IS_ZERO(beta);
   const int BetaIsOne = SCALAR_IS_ONE(beta);
   const int AlphaIsOne = SCALAR_IS_ONE(alpha);
   const int AlphaIsZero = SCALAR_IS_ZERO(alpha);
   int i, imb, nb, incAn, incx;
   void *vx=NULL, *vy=NULL;
   TYPE *x = (TYPE *) X, *y=Y;
   #ifdef TREAL
      TYPE alphaY = alpha, beta0=beta;
      #define ONE ATL_rone
      #define ZERO ATL_rzero
   #else
      const TYPE *alphaY = alpha, *beta0 = beta;
      const TYPE ONE[2]={ATL_rone, ATL_rzero}, ZERO[2]={ATL_rzero, ATL_rzero};
   #endif
   void (*getX)(const int N, const SCALAR alpha, const TYPE *X,
                const int incX, TYPE *Y, const int incY);
   void (*putY)(const int N, const SCALAR alpha, const TYPE *X,
                const int incX, const SCALAR beta, TYPE *Y, const int incY);
   void (*gemv0)(const int M, const int N, const SCALAR alpha,
                 const TYPE *A, const int lda, const TYPE *X, const int incX,
                 const SCALAR beta, TYPE *Y, const int incY);
   void (*gemv1)(const int M, const int N, const SCALAR alpha,
                 const TYPE *A, const int lda, const TYPE *X, const int incX,
                 const SCALAR beta, TYPE *Y, const int incY);
   #define gemv1_N Mjoin(PATL,gemvN_a1_x1_b1_y1)
   #define gemv1_Nc Mjoin(PATL,gemvNc_a1_x1_b1_y1)
   #define gemv1_T Mjoin(PATL,gemvT_a1_x1_b1_y1)
   #define gemv1_C Mjoin(PATL,gemvC_a1_x1_b1_y1)

   if (!M || !N) return;
   if (AlphaIsZero)  /* y <- beta * y */
   {
      if (!BetaIsOne) Mjoin(PATL,scal)(M, beta, Y, incY);
      return;
   }

   if (TA == AtlasNoTrans || TA == AtlasConj)
   {
      #ifdef ATL_AXPYMV
         gemvNaxpy(TA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
         return;
      #else
         #ifdef ATL_mvNNOBLOCK
            nb = N;
         #else
            ATL_GetPartMVN(A, lda, &i, &nb);
            if (nb > N) nb = N;
         #endif
         incAn = nb * lda SHIFT;
         #ifdef TREAL
            gemv1 = gemv1_N;
         #else
            if (TA == AtlasNoTrans) gemv1 = gemv1_N;
            else gemv1 = gemv1_Nc;
         #endif
      #endif
   }
   else
   {
      #ifdef ATL_mvTNOBLOCK
         nb = N;
      #else
         ATL_GetPartMVT(A, lda, &nb, &i);
         if (nb > N) nb = N;
      #endif
      incAn = nb SHIFT;
      #ifdef TREAL
         gemv1 = gemv1_T;
      #else
         if (TA == AtlasTrans) gemv1 = gemv1_T;
         else gemv1 = gemv1_C;
      #endif
   }

   if (incX != 1 || (incY == 1 && !AlphaIsOne && (M>>2) >= N))
   {
      vx = malloc((ATL_Cachelen<<1) + ATL_MulBySize(nb));
      ATL_assert(vx);
      x = ATL_AlignPtr(vx);
/*
 *    In trans case, adjust X to have same misalignment as A
 */
      if (TA == AtlasTrans)
      {
         ATL_FindMatAdjust(A, lda, imb);
         if (imb)
         {
            i = (int) (x - ((TYPE*)(vx)));
            #ifdef TCPLX
               i >>= 1;
            #endif
            if (imb <= i) x -= (i-imb)SHIFT;
            else x += imb SHIFT;
         }
      }
      getX = Mjoin(PATL,cpsc);
      alphaY = ONE;
   }
   else getX = NULL;
   if (incY != 1 || !SCALAR_IS_ONE(alphaY))  /* gotta buffer Y */
   {
      vy = malloc(ATL_Cachelen + ATL_MulBySize(M));
      ATL_assert(vy);
      y = ATL_AlignPtr(vy);
      putY = Mjoin(PATL,axpby);
      if (TA == AtlasNoTrans) gemv0 = Mjoin(PATL,gemvN_a1_x1_b0_y1);
      #ifdef TCPLX
         else if (TA == AtlasConj) gemv0 = Mjoin(PATL,gemvNc_a1_x1_b0_y1);
         else if (TA == AtlasConjTrans) gemv0 = Mjoin(PATL,gemvC_a1_x1_b0_y1);
      #endif
      else gemv0 = Mjoin(PATL,gemvT_a1_x1_b0_y1);
      beta0 = ZERO;
   }
   else  /* can use user's own Y */
   {
      putY = NULL;
      if (BetaIsOne)
      {
         if (TA == AtlasNoTrans) gemv0 = Mjoin(PATL,gemvN_a1_x1_b1_y1);
         #ifdef TCPLX
            else if (TA == AtlasConj)
               gemv0 = Mjoin(PATL,gemvNc_a1_x1_b1_y1);
            else if (TA == AtlasConjTrans)
               gemv0 = Mjoin(PATL,gemvC_a1_x1_b1_y1);
         #endif
         else gemv0 = Mjoin(PATL,gemvT_a1_x1_b1_y1);
      }
      else if (BetaIsZero)
      {
         if (TA == AtlasNoTrans) gemv0 = Mjoin(PATL,gemvN_a1_x1_b0_y1);
         #ifdef TCPLX
            else if (TA == AtlasConj)
               gemv0 = Mjoin(PATL,gemvNc_a1_x1_b0_y1);
            else if (TA == AtlasConjTrans)
               gemv0 = Mjoin(PATL,gemvC_a1_x1_b0_y1);
         #endif
         else gemv0 = Mjoin(PATL,gemvT_a1_x1_b0_y1);
      }
      #ifdef TCPLX
         else if (beta[1] == ATL_rzero)
         {
            if (TA == AtlasNoTrans) gemv0 = Mjoin(PATL,gemvN_a1_x1_bXi0_y1);
            else if (TA == AtlasConj) gemv0 = Mjoin(PATL,gemvNc_a1_x1_bXi0_y1);
            else if (TA == AtlasConjTrans)
               gemv0 = Mjoin(PATL,gemvC_a1_x1_bXi0_y1);
            else gemv0 = Mjoin(PATL,gemvT_a1_x1_bXi0_y1);
         }
      #endif
      else
      {
         if (TA == AtlasNoTrans) gemv0 = Mjoin(PATL,gemvN_a1_x1_bX_y1);
         #ifdef TCPLX
            else if (TA == AtlasConj) gemv0 = Mjoin(PATL,gemvNc_a1_x1_bX_y1);
            else if (TA == AtlasConjTrans)
               gemv0 = Mjoin(PATL,gemvC_a1_x1_bX_y1);
         #endif
         else gemv0 = Mjoin(PATL,gemvT_a1_x1_bX_y1);
      }
   }
   incx = nb * incX SHIFT;
   i = N;
   do
   {
      if (nb > i) nb = i;
      if (getX) getX(nb, alpha, X, incX, x, 1);
      else x = (TYPE*) X;
      gemv0(M, nb, ONE, A, lda, x, 1, beta0, y, 1);
      beta0 = ONE;
      gemv0 = gemv1;
      i -= nb;
      X += incx;
      A += incAn;
   }
   while (i);
   if (vx) free(vx);

   if (putY)
   {
      putY(M, alphaY, y, 1, beta, Y, incY);
      free(vy);
   }
}
