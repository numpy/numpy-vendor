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
#include "atlas_level1.h"
#include "atlas_level2.h"
#include "atlas_prefetch.h"

static void gemvT_Msmall(const int M, const int N, const TYPE *A, const int lda,
                         const TYPE *X, const SCALAR beta, TYPE *Y)
{
   const TYPE *stY = Y + M;
   #ifndef BETA0
      register TYPE y0;
   #endif

   do
   {
      #ifdef BETA0
         *Y = Mjoin(PATL,dot)(N, X, 1, A, 1);
      #else
         y0 = *Y;
         #ifdef BETAX
            y0 *= beta;
         #endif
         y0 += Mjoin(PATL,dot)(N, X, 1, A, 1);
         *Y = y0;
      #endif
      Y++;
      A += lda;
   }
   while (Y != stY);
}
static void gemvT2x16(const int M, const int N, const TYPE *A, const int lda,
                      const TYPE *X, const SCALAR beta, TYPE *Y)
{
   int i;
   const int M2 = (M>>1)<<1, N16 = (N>>4)<<4;
   const int incAm = (N16 ? (lda<<1) - N16 + 16 : (lda<<1));
   const int nr = (N16 ? N-N16+16 : N);
   TYPE *stY = Y + M2;
   const TYPE *stX = X + N16-16, *x;
   const TYPE *A0 = A, *A1 = A + lda;
   register TYPE y0, y1, x0, x1;
   #ifndef ATL_PFD
      #define ATL_PFD 12
   #endif

   if (M2)
   {
      do
      {
         #ifdef BETA0
            y0 = y1 = ATL_rzero;
         #elif defined(BETAX)
            y0 = *Y * beta; y1 = beta * Y[1];
         #else
            y0 = *Y; y1 = Y[1];
         #endif
         x = X;
         if (N16 > 16)
         {
            x0 = *x; x1 = x[1];
            do
            {
               y0 += x0 * *A0; ATL_pfl1R(A0+ATL_PFD); ATL_pfl1R(A1+ATL_PFD);
                               ATL_pfl1R(A0+ATL_PFD+8); ATL_pfl1R(A1+ATL_PFD+8);
               y1 += x0 * *A1; x0 = x[2];
               y0 += x1 * A0[1];
               y1 += x1 * A1[1]; x1 = x[3];
               y0 += x0 * A0[2];
               y1 += x0 * A1[2]; x0 = x[4];
               y0 += x1 * A0[3];
               y1 += x1 * A1[3]; x1 = x[5];
               y0 += x0 * A0[4];
               y1 += x0 * A1[4]; x0 = x[6];
               y0 += x1 * A0[5];
               y1 += x1 * A1[5]; x1 = x[7];
               y0 += x0 * A0[6];
               y1 += x0 * A1[6]; x0 = x[8];
               y0 += x1 * A0[7];
               y1 += x1 * A1[7]; x1 = x[9];
               y0 += x0 * A0[8];
                  ATL_pfl1R(A0+ATL_PFD+4); ATL_pfl1R(A1+ATL_PFD+4);
                  ATL_pfl1R(A0+ATL_PFD+12); ATL_pfl1R(A1+ATL_PFD+12);
               y1 += x0 * A1[8]; x0 = x[10];
               y0 += x1 * A0[9];
               y1 += x1 * A1[9]; x1 = x[11];
               y0 += x0 * A0[10];
               y1 += x0 * A1[10]; x0 = x[12];
               y0 += x1 * A0[11];
               y1 += x1 * A1[11]; x1 = x[13];
               y0 += x0 * A0[12];
               y1 += x0 * A1[12]; x0 = x[14];
               y0 += x1 * A0[13];
               y1 += x1 * A1[13]; x1 = x[15];
               y0 += x0 * A0[14]; x += 16;
               y1 += x0 * A1[14]; x0 = *x;
               y0 += x1 * A0[15]; A0 += 16;
               y1 += x1 * A1[15]; x1 = x[1]; A1 += 16;
            }
            while(x != stX);
         }
         for (i=0; i != nr; i++)
         {
            x0 = x[i];
            y0 += x0 * A0[i];
            y1 += x0 * A1[i];
         }
         A0 += incAm;
         *Y = y0;
         A1 += incAm;
         Y[1] = y1;
         Y += 2;
      }
      while (Y != stY);
   }
   if (M-M2) gemvT_Msmall(M-M2, N, A0, lda, X, beta, Y);
}
#define gemv0 Mjoin(Mjoin(Mjoin(Mjoin(gemvT,NM),_x1),BNM),_y1)
void Mjoin(PATL,gemv0)
   (const int M, const int N, const SCALAR alpha,
    const TYPE *A, const int lda, const TYPE *X, const int incX,
    const SCALAR beta, TYPE *Y, const int incY)
{
   gemvT2x16(M, N, A, lda, X, beta, Y);
}
