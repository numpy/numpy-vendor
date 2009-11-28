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
#include "atlas_lvl3.h"
#include "atlas_r1.h"


#ifdef TREAL
   #define ATL_ger Mjoin(PATL,ger)
   #define ATL_ger1 Mjoin(PATL,ger1_a1_x1_yX)
#else
   #ifdef Conj_
      #define ATL_ger Mjoin(PATL,gerc)
   #else
      #define ATL_ger Mjoin(PATL,geru)
      #define ATL_ger1 Mjoin(PATL,ger1u_a1_x1_yX)
   #endif
#endif
void ATL_ger(const int M, const int N, const SCALAR alpha,
             const TYPE *X, const int incX, const TYPE *Y, const int incY,
             TYPE *A, const int lda)
{
   int imb, mb, mb0, m=M, i;
   int incy=incY;
   #ifdef TREAL
      #define one ATL_rone
   #else
      static TYPE one[2] = {ATL_rone, ATL_rzero};
   #endif
   void *vx=NULL;
   TYPE *x, *y = (TYPE*) Y;
   void (*getX)(const int N, const SCALAR alpha, const TYPE *X,
                const int incX, TYPE *Y, const int incY);
   #ifdef Conj_
      void (*ATL_ger1)(const int M, const int N, const SCALAR alpha,
                       const TYPE *X, const int incX, const TYPE *Y,
                       const int incY, TYPE *A, const int lda);
      ATL_ger1 = Mjoin(PATL,ger1c_a1_x1_yX);
   #endif

   if ( !M || !N || SCALAR_IS_ZERO(alpha) ) return;

   #ifdef ATL_R1NOBLOCK
      imb = mb = M;
   #else
      ATL_FindMatAdjust(A, lda, imb);
      ATL_GetPartR1(A, lda, mb, i);
   #endif

   if (incX != 1 || !SCALAR_IS_ONE(alpha))
   {
/*
 *    Apply alpha to Y if X has stride 1 & Y is MUCH smaller
 *    The LAPACK barfs if you  switch which vector alpha is applied to,
 *    since it tests tiny matrices, so make it apply to X when they are
 *    close to even
 */
      if (incX == 1 && N < (M>>4))
      {
         vx = malloc(ATL_Cachelen + ATL_MulBySize(N));
         ATL_assert(vx);
         y = ATL_AlignPtr(vx);
         #ifdef Conj_
            Mjoin(PATL,moveConj)(N, alpha, Y, incY, y, 1);
            ATL_ger1 = Mjoin(PATL,ger1u_a1_x1_yX);
         #else
            Mjoin(PATL,cpsc)(N, alpha, Y, incY, y, 1);
         #endif
         incy = 1;
         getX = NULL;
      }
      else
      {
         i = Mmax(imb,mb);
         i = Mmin(i,M);
         vx = malloc(ATL_Cachelen + ATL_MulBySize(i));
         ATL_assert(vx);
         x = ATL_AlignPtr(vx);
         getX = Mjoin(PATL,cpsc);
      }
   }
   else getX = NULL;

   if (imb) mb0 = Mmin(imb,m);
   else mb0 = Mmin(mb,m);
   do
   {
      if (getX) getX(mb0, alpha, X, incX, x, 1);
      else x = (TYPE*) X;
      ATL_ger1(mb0, N, one, x, 1, y, incy, A, lda);
      A += mb0 SHIFT;
      X += mb0*incX SHIFT;
      m -= mb0;
      mb0 = Mmin(m,mb);
   }
   while(m);
   if (vx) free(vx);
}
