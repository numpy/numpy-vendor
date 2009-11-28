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

static void gemvT_Nsmall(const int M, const int N, const TYPE *A, const int lda,
                         const TYPE *X, const SCALAR beta, TYPE *Y)
{
   int i;
   register TYPE x0, x1, x2, x3, x4, x5, x6, x7, y0;
   TYPE *stY = Y + M;

   switch(N)
   {
   case 1:
      #ifdef BETA0
         Mjoin(PATL,cpsc)(M, *X, A, lda, Y, 1);
      #elif defined(BETAX)
         Mjoin(PATL,axpby)(M, *X, A, lda, beta, Y, 1);
      #else
         Mjoin(PATL,axpy)(M, *X, A, lda, Y, 1);
      #endif
      break;
   case 2:
      x0 = *X;
      x1 = X[1];
      do
      {
         #ifdef BETA0
            *Y = x0 * *A + x1 * A[1];
         #elif defined(BETAX)
            y0 = *Y;
            *Y = y0 * beta + x0 * *A + x1 * A[1];
         #else
            *Y += x0 * *A + x1 * A[1];
         #endif
         A += lda;
      }
      while (++Y != stY);
      break;
   case 3:
      x0 = *X;
      x1 = X[1];
      x2 = X[2];
      do
      {
         #ifdef BETA0
            *Y = x0 * *A + x1 * A[1] + x2 * A[2];
         #elif defined(BETAX)
            y0 = *Y;
            *Y = y0 * beta + x0 * *A + x1 * A[1] + x2 * A[2];
         #else
            *Y += x0 * *A + x1 * A[1] + x2 * A[2];
         #endif
         A += lda;
      }
      while (++Y != stY);
      break;
   case 4:
      x0 = *X;
      x1 = X[1];
      x2 = X[2];
      x3 = X[3];
      do
      {
         #ifdef BETA0
            *Y = x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3];
         #elif defined(BETAX)
            y0 = *Y;
            *Y = y0 * beta + x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3];
         #else
            *Y += x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3];
         #endif
         A += lda;
      }
      while (++Y != stY);
      break;
   case 5:
      x0 = *X;
      x1 = X[1];
      x2 = X[2];
      x3 = X[3];
      x4 = X[4];
      do
      {
         #ifdef BETA0
            *Y = x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4];
         #elif defined(BETAX)
            y0 = *Y;
            *Y = y0 * beta + x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4];
         #else
            *Y += x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3] + x4 * A[4];
         #endif
         A += lda;
      }
      while (++Y != stY);
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
         #ifdef BETA0
            *Y = x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4] + x5 * A[5];
         #elif defined(BETAX)
            y0 = *Y;
            *Y = y0 * beta + x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4] + x5 * A[5];
         #else
            *Y += x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4] + x5 * A[5];
         #endif
         A += lda;
      }
      while (++Y != stY);
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
         #ifdef BETA0
            *Y = x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4] + x5 * A[5] + x6 * A[6];
         #elif defined(BETAX)
            y0 = *Y;
            *Y = y0 * beta + x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4] + x5 * A[5] + x6 * A[6];
         #else
            *Y += x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4] + x5 * A[5] + x6 * A[6];
         #endif
         A += lda;
      }
      while (++Y != stY);
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
         #ifdef BETA0
            *Y = x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4] + x5 * A[5] + x6 * A[6] + x7 * A[7];
         #elif defined(BETAX)
            y0 = *Y;
            *Y = y0 * beta + x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4] + x5 * A[5] + x6 * A[6] + x7 * A[7];
         #else
            *Y += x0 * *A + x1 * A[1] + x2 * A[2] + x3 * A[3]
                 + x4 * A[4] + x5 * A[5] + x6 * A[6] + x7 * A[7];
         #endif
         A += lda;
      }
      while (++Y != stY);
      break;
   default:
      if (M)
      {
         do
         {
            #ifdef BETA0
               y0 = ATL_rzero;
            #elif defined(BETAX)
               y0 = *Y * beta;
            #else
               y0 = *Y;
            #endif
            for (i=0; i != N; i++) y0 += A[i] * X[i];
            *Y++ = y0;
            A += lda;
         }
         while (Y != stY);
      }
   }
}

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
static void gemvT4x4(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
{
   const int M4 = (M>>2)<<2, N4 = (N>>2)<<2;
   const int nr = N-N4, incAm = (lda<<2) - N4;
   const TYPE *stX = X + N4, *x;
   const TYPE *A0 = A, *A1 = A + lda, *A2 = A1 + lda, *A3 = A2 + lda;
   TYPE *stY = Y + M4;
   register TYPE a00, a10, a20, a30, a01, a11, a21, a31;
   register TYPE a02, a12, a22, a32, a03, a13, a23, a33;
   register TYPE y0, y1, y2, y3, yy0, yy1, yy2, yy3;
   register TYPE x0, x1, x2, x3;

   if (N4 > 4)
   {
      if (M4)
      {
         do
         {
            a00 = *A0;   a01 = *A1;   a02 = *A2;   a03 = *A3;
            a10 = A0[1]; a11 = A1[1]; a12 = A2[1]; a13 = A3[1];
            a20 = A0[2]; a21 = A1[2]; a22 = A2[2]; a23 = A3[2];
            a30 = A0[3]; a31 = A1[3]; a32 = A2[3]; a33 = A3[3];
            A0 += 4;     A1 += 4;     A2 += 4;     A3 += 4;
            #ifdef BETA0
               y0 = y1 = y2 = y3 = yy0 = yy1 = yy2 = yy3 = ATL_rzero;
            #elif defined(BETAX)
               y0 = beta;
               yy0 = *Y; yy1 = Y[1]; yy2 = Y[2];  yy3 = Y[3];
               yy0 *= y0; yy1 *= y0; yy2 *= y0; yy3 *= y0;
               y0 = y1 = y2 = y3 = ATL_rzero;
            #else
               yy0 = *Y; yy1 = Y[1]; yy2 = Y[2];  yy3 = Y[3];
               y0 = y1 = y2 = y3 = ATL_rzero;
            #endif
            x0 = *X; x1 = X[1]; x2 = X[2]; x3 = X[3]; x = X + 4;
            do
            {
               y0  += a00 * x0; a00 = *A0;
               y1  += a01 * x0; a01 = *A1;
               yy0 += a10 * x1; a10 = A0[1];
               yy1 += a11 * x1; a11 = A1[1];
               y2  += a02 * x0; a02 = *A2;
               y3  += a03 * x0; a03 = *A3; x0 = *x;
               yy2 += a12 * x1; a12 = A2[1];
               yy3 += a13 * x1; a13 = A3[1]; x1 = x[1];

               y0  += a20 * x2; a20 = A0[2];
               y1  += a21 * x2; a21 = A1[2];
               yy0 += a30 * x3; a30 = A0[3]; A0 += 4;
               yy1 += a31 * x3; a31 = A1[3]; A1 += 4;
               y2  += a22 * x2; a22 = A2[2];
               y3  += a23 * x2; a23 = A3[2]; x2 = x[2];
               yy2 += a32 * x3; a32 = A2[3]; A2 += 4;
               yy3 += a33 * x3; a33 = A3[3]; A3 += 4; x3 = x[3]; x += 4;
            }
            while (x != stX);
            y0  += a00 * x0;
            y1  += a01 * x0;
            yy0 += a10 * x1;
            yy1 += a11 * x1;
            y2  += a02 * x0;
            y3  += a03 * x0;
            yy2 += a12 * x1;
            yy3 += a13 * x1;

            y0  += a20 * x2;
            y1  += a21 * x2;
            yy0 += a30 * x3;
            yy1 += a31 * x3;
            y2  += a22 * x2;
            y3  += a23 * x2;
            yy2 += a32 * x3;
            yy3 += a33 * x3;
            switch(nr)
            {
            case 1:
               x0 = *x;
               y0 += *A0 * x0;
               y1 += *A1 * x0;
               y2 += *A2 * x0;
               y3 += *A3 * x0;
               break;
            case 2:
               x0 = *x;
               x1 = x[1];
               y0  += *A0 * x0;
               y1  += *A1 * x0;
               yy0 += A0[1] * x1;
               yy1 += A1[1] * x1;
               y2  += *A2 * x0;
               y3  += *A3 * x0;
               yy2 += A2[1] * x1;
               yy3 += A3[1] * x1;
               break;
            case 3:
               x0 = *x;
               x1 = x[1];
               x2 = x[2];
               y0  += *A0 * x0;
               y1  += *A1 * x0;
               yy0 += A0[1] * x1;
               yy1 += A1[1] * x1;
               y2  += *A2 * x0;
               y3  += *A3 * x0;
               yy2 += A2[1] * x1;
               yy3 += A3[1] * x1;
               y0  += A0[2] * x2;
               y1  += A1[2] * x2;
               y2  += A2[2] * x2;
               y3  += A3[2] * x2;
               break;
            default:;
            }
            A0 += incAm;
            y0 += yy0;
            A1 += incAm;
            y1 += yy1;
            A2 += incAm;
            y2 += yy2;
            A3 += incAm;
            y3 += yy3;
            *Y = y0;
            Y[1] = y1;
            Y[2] = y2;
            Y[3] = y3;
            Y += 4;
         }
         while (Y != stY);
      }
      if (M-M4) gemvT_Msmall(M-M4, N, A0, lda, X, beta, Y);
   }
   else if (M) gemvT_Nsmall(M, N, A, lda, X, beta, Y);
}
#define gemv0 Mjoin(Mjoin(Mjoin(Mjoin(gemvT,NM),_x1),BNM),_y1)
void Mjoin(PATL,gemv0)
   (const int M, const int N, const SCALAR alpha,
    const TYPE *A, const int lda, const TYPE *X, const int incX,
    const SCALAR beta, TYPE *Y, const int incY)
{
   gemvT4x4(M, N, A, lda, X, beta, Y);
}
