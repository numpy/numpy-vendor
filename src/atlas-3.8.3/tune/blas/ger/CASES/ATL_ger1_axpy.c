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
#include "atlas_prefetch.h"

static void ger_Mle8(const int M, const int N, const TYPE *X,
                     const TYPE *Y, const int incY, TYPE *A, const int lda)
{
   const TYPE *stY = Y + incY*N;
   register TYPE y0, x0, x1, x2, x3, x4, x5, x6, x7;

   switch(M)
   {
   case 1:
      x0 = *X;
      do
      {
         *A += *Y * x0;
         A += lda;
         Y += incY;
      }
      while (Y != stY);
      break;
   case 2:
      x0 = *X;
      x1 = X[1];
      do
      {
         y0 = *Y;
         *A += y0 * x0; Y += incY;
         A[1] += y0 * x1;
         A += lda;
      }
      while (Y != stY);
      break;
   case 3:
      x0 = *X;
      x1 = X[1];
      x2 = X[2];
      do
      {
         y0 = *Y;
         *A += y0 * x0;
         A[1] += y0 * x1; Y += incY;
         A[2] += y0 * x2; A += lda;
      }
      while (Y != stY);
      break;
   case 4:
      x0 = *X;
      x1 = X[1];
      x2 = X[2];
      x3 = X[3];
      do
      {
         y0 = *Y;
         *A += y0 * x0;
         A[1] += y0 * x1; Y += incY;
         A[2] += y0 * x2;
         A[3] += y0 * x3; A += lda;
      }
      while (Y != stY);
      break;
   case 5:
      x0 = *X;
      x1 = X[1];
      x2 = X[2];
      x3 = X[3];
      x4 = X[4];
      do
      {
         y0 = *Y;
         *A += y0 * x0;
         A[1] += y0 * x1; Y += incY;
         A[2] += y0 * x2;
         A[3] += y0 * x3;
         A[4] += y0 * x4; A += lda;
      }
      while (Y != stY);
      break;
   case 6:
      x0 = *X;
      x1 = X[1];
      x2 = X[2];
      x3 = X[3];
      x4 = X[4];
      x5 = X[5];
      do
      {
         y0 = *Y;
         *A += y0 * x0;
         A[1] += y0 * x1; Y += incY;
         A[2] += y0 * x2;
         A[3] += y0 * x3;
         A[4] += y0 * x4;
         A[5] += y0 * x5; A += lda;
      }
      while (Y != stY);
      break;
   case 7:
      x0 = *X;
      x1 = X[1];
      x2 = X[2];
      x3 = X[3];
      x4 = X[4];
      x5 = X[5];
      x6 = X[6];
      do
      {
         y0 = *Y;
         *A += y0 * x0;
         A[1] += y0 * x1;
         A[2] += y0 * x2; Y += incY;
         A[3] += y0 * x3;
         A[4] += y0 * x4;
         A[5] += y0 * x5;
         A[6] += y0 * x6; A += lda;
      }
      while (Y != stY);
      break;
   case 8:
      x0 = *X;
      x1 = X[1];
      x2 = X[2];
      x3 = X[3];
      x4 = X[4];
      x5 = X[5];
      x6 = X[6];
      x7 = X[7];
      do
      {
         y0 = *Y;
         *A += y0 * x0;
         A[1] += y0 * x1;
         A[2] += y0 * x2; Y += incY;
         A[3] += y0 * x3;
         A[4] += y0 * x4;
         A[5] += y0 * x5;
         A[6] += y0 * x6;
         A[7] += y0 * x7; A += lda;
      }
      while (Y != stY);
      break;
   default:
      ATL_assert(M == 0);
   }
}
static void ger_Nle4(const int M, const int N, const TYPE *X,
                     const TYPE *Y, const int incY, TYPE *A, const int lda)
{
   register TYPE y0, y1, y2, y3, x0;
   TYPE *A0 = A, *A1 = A+lda, *A2 = A1+lda, *A3 = A2+lda;
   int i;

   switch(N)
   {
   case 1:
      y0 = *Y;
      for (i=0; i != M; i++) A0[i] += y0 * X[i];
      break;
   case 2:
      y0 = *Y; y1 = Y[incY];
      for (i=0; i != M; i++)
      {
         x0 = X[i];
         A0[i] += y0 * x0;
         A1[i] += y1 * x0;
      }
      break;
   case 3:
      y0 = *Y; y1 = Y[incY]; y2 = Y[incY<<1];
      for (i=0; i != M; i++)
      {
         x0 = X[i];
         A0[i] += y0 * x0;
         A1[i] += y1 * x0;
         A2[i] += y2 * x0;
      }
      break;
   case 4:
      y0 = *Y; y1 = Y[incY]; y2 = Y[incY+incY]; y3 = Y[(incY<<1)+incY];
      for (i=0; i != M; i++)
      {
         x0 = X[i];
         A0[i] += y0 * x0;
         A1[i] += y1 * x0;
         A2[i] += y2 * x0;
         A3[i] += y3 * x0;
      }
      break;
   default:;
   }
}
void Mjoin(PATL,ger1_a1_x1_yX)
   (const int M, const int N, const SCALAR alpha, const TYPE *X, const int incX,
    const TYPE *Y, const int incY, TYPE *A, const int lda)
{
   const TYPE *stY = Y + N*incY;
   if (M > 8)
   {
      do
      {
         Mjoin(PATL,axpy)(M, *Y, X, 1, A, 1);
         Y += incY;
         A += lda;
      }
      while (Y != stY);
   }
   else ger_Mle8(M, N, X, Y, incY, A, lda);
}
