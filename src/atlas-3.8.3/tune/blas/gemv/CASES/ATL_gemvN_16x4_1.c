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
#include "atlas_lvl2.h"


static void gemvMlt8(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y);
static void gemvNle4(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y);
#ifdef BETA0
   #define Yget(y_, yp_, bet_) (y_) = ATL_rzero
#elif defined BETAX
   #define Yget(y_, yp_, bet_) (y_) = (yp_) * (bet_)
#else
   #define Yget(y_, yp_, bet_) (y_) = (yp_)
#endif
static void gemvN32x4(const int M, const int N, const TYPE *A, const int lda,
                      const TYPE *x, const SCALAR beta0, TYPE *y)
/*
 * rank-4 daxpy based NoTrans gemv
 */
{
   const int M16 = (M>>4)<<4;
   TYPE *stY = y + M16 - 32;
   const TYPE *A0 = A, *A1 = A+lda, *A2 = A1 + lda, *A3 = A2 + lda;
   register TYPE z0, z1, z2, z3, z4, z5, z6, z7;
   register TYPE y0, y1, y2, y3, y4, y5, y6, y7;
   const register TYPE x0 = *x, x1 = x[1], x2 = x[2], x3 = x[3];
   #ifdef BETAX
      const register TYPE beta = beta0;
   #else
      #define beta beta0
   #endif

   ATL_assert(N == 4);
   if (M16 >= 32)
   {
      #ifdef BETA0
         y0 = y1 = y2 = y3 = y4 = y5 = y6 = y7 = ATL_rzero;
      #else
         y0 = *y;   y1 = y[1]; y2 = y[2]; y3 = y[3];
         y4 = y[4]; y5 = y[5]; y6 = y[6]; y7 = y[7];
         #ifdef BETAX
            y0 *= beta; y1 *= beta; y2 *= beta; y3 *= beta;
            y4 *= beta; y5 *= beta; y6 *= beta; y7 *= beta;
         #endif
      #endif
      y0 += x0 * *A0;   Yget(z0, y[8], beta);
      y1 += x1 * A1[1];
      y2 += x2 * A2[2]; Yget(z1, y[9], beta);
      y3 += x3 * A3[3];
      y4 += x0 * A0[4]; Yget(z2, y[10], beta);
      y5 += x1 * A1[5];
      y6 += x2 * A2[6]; Yget(z3, y[11], beta);
      y7 += x3 * A3[7];

      y0 += x1 * *A1;   Yget(z4, y[12], beta);
      y1 += x2 * A2[1];
      y2 += x3 * A3[2]; Yget(z5, y[13], beta);
      y3 += x0 * A0[3];
      y4 += x1 * A1[4]; Yget(z6, y[14], beta);
      y5 += x2 * A2[5];
      y6 += x3 * A3[6]; Yget(z7, y[15], beta);
      y7 += x0 * A0[7];

      y0 += x2 * *A2;
      y1 += x3 * A3[1];
      y2 += x0 * A0[2];
      y3 += x1 * A1[3];
      y4 += x2 * A2[4];
      y5 += x3 * A3[5];
      y6 += x0 * A0[6];
      y7 += x1 * A1[7];

      y0 += x3 * *A3;
      y1 += x0 * A0[1];
      y2 += x1 * A1[2];
      y3 += x2 * A2[3];
      y4 += x3 * A3[4];
      y5 += x0 * A0[5];
      y6 += x1 * A1[6];
      y7 += x2 * A2[7];

      z0 += x0 * A0[8]; *y = y0;
      z1 += x1 * A1[9];
      z2 += x2 * A2[10]; y[1] = y1;
      z3 += x3 * A3[11];
      z4 += x0 * A0[12]; y[2] = y2;
      z5 += x1 * A1[13];
      z6 += x2 * A2[14]; y[3] = y3;
      z7 += x3 * A3[15];

      z0 += x1 * A1[8]; y[4] = y4;
      z1 += x2 * A2[9];
      z2 += x3 * A3[10]; y[5] = y5;
      z3 += x0 * A0[11];
      z4 += x1 * A1[12]; y[6] = y6;
      z5 += x2 * A2[13];
      z6 += x3 * A3[14]; y[7] = y7;
      z7 += x0 * A0[15];

      z0 += x2 * A2[8];  Yget(y0, y[16], beta);
      z1 += x3 * A3[ 9];
      z2 += x0 * A0[10]; Yget(y1, y[17], beta);
      z3 += x1 * A1[11];
      z4 += x2 * A2[12]; Yget(y2, y[18], beta);
      z5 += x3 * A3[13];
      z6 += x0 * A0[14]; Yget(y3, y[19], beta);
      z7 += x1 * A1[15];

      z0 += x3 * A3[8];  Yget(y4, y[20], beta);
      z1 += x0 * A0[9];
      z2 += x1 * A1[10]; Yget(y5, y[21], beta);
      z3 += x2 * A2[11];
      z4 += x3 * A3[12]; Yget(y6, y[22], beta); A3 += 16;
      z5 += x0 * A0[13];                        A0 += 16;
      z6 += x1 * A1[14]; Yget(y7, y[23], beta); A1 += 16;
      z7 += x2 * A2[15];                        A2 += 16;
      if (M16 != 32)
      {
         do
         {
            y0 += x0 * *A0;   y[8] = z0;
            y1 += x1 * A1[1];
            y2 += x2 * A2[2]; y[9] = z1;
            y3 += x3 * A3[3];
            y4 += x0 * A0[4]; y[10] = z2;
            y5 += x1 * A1[5];
            y6 += x2 * A2[6]; y[11] = z3;
            y7 += x3 * A3[7];

            y0 += x1 * *A1;   y[12] = z4;
            y1 += x2 * A2[1];
            y2 += x3 * A3[2]; y[13] = z5;
            y3 += x0 * A0[3];
            y4 += x1 * A1[4]; y[14] = z6;
            y5 += x2 * A2[5];
            y6 += x3 * A3[6]; y[15] = z7; y += 16;
            y7 += x0 * A0[7];

            y0 += x2 * *A2;   Yget(z0, y[8], beta);
            y1 += x3 * A3[1];
            y2 += x0 * A0[2]; Yget(z1, y[9], beta);
            y3 += x1 * A1[3];
            y4 += x2 * A2[4]; Yget(z2, y[10], beta);
            y5 += x3 * A3[5];
            y6 += x0 * A0[6]; Yget(z3, y[11], beta);
            y7 += x1 * A1[7];

            y0 += x3 * *A3;   Yget(z4, y[12], beta);
            y1 += x0 * A0[1];
            y2 += x1 * A1[2]; Yget(z5, y[13], beta);
            y3 += x2 * A2[3];
            y4 += x3 * A3[4]; Yget(z6, y[14], beta);
            y5 += x0 * A0[5];
            y6 += x1 * A1[6]; Yget(z7, y[15], beta);
            y7 += x2 * A2[7];

            z0 += x0 * A0[8];  *y = y0;
            z1 += x1 * A1[9];
            z2 += x2 * A2[10]; y[1] = y1;
            z3 += x3 * A3[11];
            z4 += x0 * A0[12]; y[2] = y2;
            z5 += x1 * A1[13];
            z6 += x2 * A2[14]; y[3] = y3;
            z7 += x3 * A3[15];

            z0 += x1 * A1[8]; y[4] = y4;
            z1 += x2 * A2[9];
            z2 += x3 * A3[10]; y[5] = y5;
            z3 += x0 * A0[11];
            z4 += x1 * A1[12]; y[6] = y6;
            z5 += x2 * A2[13];
            z6 += x3 * A3[14]; y[7] = y7;
            z7 += x0 * A0[15];

            z0 += x2 * A2[8];  Yget(y0, y[16], beta);
            z1 += x3 * A3[ 9];
            z2 += x0 * A0[10]; Yget(y1, y[17], beta);
            z3 += x1 * A1[11];
            z4 += x2 * A2[12]; Yget(y2, y[18], beta);
            z5 += x3 * A3[13];
            z6 += x0 * A0[14]; Yget(y3, y[19], beta);
            z7 += x1 * A1[15];

            z0 += x3 * A3[8];  Yget(y4, y[20], beta);
            z1 += x0 * A0[9];
            z2 += x1 * A1[10]; Yget(y5, y[21], beta);
            z3 += x2 * A2[11];
            z4 += x3 * A3[12]; Yget(y6, y[22], beta);  A3 += 16;
            z5 += x0 * A0[13];                         A0 += 16;
            z6 += x1 * A1[14]; Yget(y7, y[23], beta);  A1 += 16;
            z7 += x2 * A2[15];                         A2 += 16;
         }
         while (y != stY);
      }
      y0 += x0 * *A0;   y[8] = z0;
      y1 += x1 * A1[1];
      y2 += x2 * A2[2]; y[9] = z1;
      y3 += x3 * A3[3];
      y4 += x0 * A0[4]; y[10] = z2;
      y5 += x1 * A1[5];
      y6 += x2 * A2[6]; y[11] = z3;
      y7 += x3 * A3[7];

      y0 += x1 * *A1;   y[12] = z4;
      y1 += x2 * A2[1];
      y2 += x3 * A3[2]; y[13] = z5;
      y3 += x0 * A0[3];
      y4 += x1 * A1[4]; y[14] = z6;
      y5 += x2 * A2[5];
      y6 += x3 * A3[6]; y[15] = z7; y += 16;
      y7 += x0 * A0[7];

      y0 += x2 * *A2;   Yget(z0, y[8], beta);
      y1 += x3 * A3[1];
      y2 += x0 * A0[2]; Yget(z1, y[9], beta);
      y3 += x1 * A1[3];
      y4 += x2 * A2[4]; Yget(z2, y[10], beta);
      y5 += x3 * A3[5];
      y6 += x0 * A0[6]; Yget(z3, y[11], beta);
      y7 += x1 * A1[7];

      y0 += x3 * *A3;   Yget(z4, y[12], beta);
      y1 += x0 * A0[1];
      y2 += x1 * A1[2]; Yget(z5, y[13], beta);
      y3 += x2 * A2[3];
      y4 += x3 * A3[4]; Yget(z6, y[14], beta);
      y5 += x0 * A0[5];
      y6 += x1 * A1[6]; Yget(z7, y[15], beta);
      y7 += x2 * A2[7];

      z0 += x0 * A0[8];  *y = y0;
      z1 += x1 * A1[9];
      z2 += x2 * A2[10]; y[1] = y1;
      z3 += x3 * A3[11];
      z4 += x0 * A0[12]; y[2] = y2;
      z5 += x1 * A1[13];
      z6 += x2 * A2[14]; y[3] = y3;
      z7 += x3 * A3[15];

      z0 += x1 * A1[8]; y[4] = y4;
      z1 += x2 * A2[9];
      z2 += x3 * A3[10]; y[5] = y5;
      z3 += x0 * A0[11];
      z4 += x1 * A1[12]; y[6] = y6;
      z5 += x2 * A2[13];
      z6 += x3 * A3[14]; y[7] = y7;
      z7 += x0 * A0[15];

      z0 += x2 * A2[8];
      z1 += x3 * A3[ 9];
      z2 += x0 * A0[10];
      z3 += x1 * A1[11];
      z4 += x2 * A2[12];
      z5 += x3 * A3[13];
      z6 += x0 * A0[14];
      z7 += x1 * A1[15];

      z0 += x3 * A3[8];
      z1 += x0 * A0[9];
      z2 += x1 * A1[10];
      z3 += x2 * A2[11];
      z4 += x3 * A3[12];
      z5 += x0 * A0[13];
      z6 += x1 * A1[14];
      z7 += x2 * A2[15];
      y[8] = z0;
      y[9] = z1;
      y[10] = z2;
      y[11] = z3;
      y[12] = z4;
      y[13] = z5;
      y[14] = z6;
      y[15] = z7;
      if (M-M16) gemvMlt8(M-M16, N, A0+16, lda, x, beta, y+16);
   }
   else if (N) gemvMlt8(M, N, A, lda, x, beta, y);
}

static void gemv32x4(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
{
   #ifdef BETA1
      int j;
   #endif
   const int incA = lda<<2;

   if (N >= 4)
   {
      if (M >= 32)
      {
         #ifdef BETA1
            for (j=(N>>2); j; j--, A += incA, X += 4)
               gemvN32x4(M, 4, A, lda, X, ATL_rone, Y);
            if ( (j = N-((N>>2)<<2)) ) gemvNle4(M, j, A, lda, X, ATL_rone, Y);
         #else
            gemvN32x4(M, 4, A, lda, X, beta, Y);
            if (N != 4)
               Mjoin(PATL,gemvN_a1_x1_b1_y1)
                  (M, N-4, ATL_rone, A+incA, lda, X+4, 1, ATL_rone, Y, 1);
         #endif
      }
      else gemvMlt8(M, N, A, lda, X, beta, Y);
   }
   else if (M) gemvNle4(M, N, A, lda, X, beta, Y);
}

static void gemvMlt8(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
{
   int i;
   register TYPE y0;
   for (i=M; i; i--)
   {
      #ifdef BETA0
         y0 = Mjoin(PATL,dot)(N, A, lda, X, 1);
      #else
         Yget(y0, *Y, beta);
         y0 += Mjoin(PATL,dot)(N, A, lda, X, 1);
      #endif
      *Y++ = y0;
      A++;
   }
}
static void gemvNle4(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
{
   int i;
   const TYPE *A0 = A, *A1 = A+lda, *A2 = A1+lda, *A3 = A2+lda;
   register TYPE x0, x1, x2, x3;
   #ifdef BETAX
      const register TYPE bet=beta;
   #endif

   switch(N)
   {
   case 1:
      #if defined(BETA0)
         Mjoin(PATL,cpsc)(M, *X, A, 1, Y, 1);
      #elif defined(BETAX)
         Mjoin(PATL,axpby)(M, *X, A, 1, beta, Y, 1);
      #else
         Mjoin(PATL,axpy)(M, *X, A, 1, Y, 1);
      #endif
      break;
   case 2:
      x0 = *X; x1 = X[1];
      for (i=0; i != M; i++)
      #ifdef BETA0
         Y[i] = A0[i] * x0 + A1[i] * x1;
      #elif defined(BETAX)
         Y[i] = Y[i]*bet + A0[i] * x0 + A1[i] * x1;
      #else
         Y[i] += A0[i] * x0 + A1[i] * x1;
      #endif
      break;
   case 3:
      x0 = *X; x1 = X[1]; x2 = X[2];
      for (i=0; i != M; i++)
      #ifdef BETA0
         Y[i] = A0[i] * x0 + A1[i] * x1 + A2[i] * x2;
      #elif defined(BETAX)
         Y[i] = Y[i]*bet + A0[i] * x0 + A1[i] * x1 + A2[i] * x2;
      #else
         Y[i] += A0[i] * x0 + A1[i] * x1 + A2[i] * x2;
      #endif
      break;
   case 4:
      if (M >= 32) gemv32x4(M, 4, A, lda, X, beta, Y);
      else
      {
         x0 = *X; x1 = X[1]; x2 = X[2]; x3 = X[3];
         for (i=0; i != M; i++)
         #ifdef BETA0
            Y[i] = A0[i] * x0 + A1[i] * x1 + A2[i] * x2 + A3[i] * x3;
         #elif defined(BETAX)
            Y[i] = Y[i]*bet + A0[i] * x0 + A1[i] * x1 + A2[i] * x2 + A3[i] * x3;
         #else
            Y[i] += A0[i] * x0 + A1[i] * x1 + A2[i] * x2 + A3[i] * x3;
         #endif
      }
      break;
   default:
      ATL_assert(!N);
   }
}

static void gemv16x4(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
/*
 * 16x4 with feeble prefetch
 */
{
   int j;
   const int M16 = (M>>4)<<4, N4 = (N>>2)<<2, nr = N-N4+4;
   const int incA = lda << 2, incAm = 16 - (N4-4)*lda, incAm0 = 16 - N*lda;
   const TYPE *stX = X + N4, *x;
   const TYPE *A0 = A, *A1 = A + lda, *A2 = A1 + lda, *A3 = A2 + lda;
   TYPE *stY = Y + M16;
   register TYPE x0, x1, x2, x3;
   register TYPE y0, y1, y2, y3, y4, y5, y6, y7;
   register TYPE y8, y9, y10, y11, y12, y13, y14, y15;
   register TYPE p0, p1, p2, p3;

   if (N > 4)
   {
      if (M16)
      {
         do
         {
            #ifdef BETA0
               y0 = y1 = y2 = y3 = y4 = y5 = y6 = y7 =
               y8 = y9 = y10 = y11 = y12 = y13 = y14 = y15 = ATL_rzero;
            #elif defined BETAX
               x0 = beta;
               y0 = *Y; y1 = Y[1]; y2 = Y[2];  y3 = Y[3];
               y8 = Y[8]; y9 = Y[9]; y10 = Y[10];  y11 = Y[11];
               y4 = Y[4]; y5 = Y[5]; y6 = Y[6];  y7 = Y[7];
               y12 = Y[12]; y13 = Y[13]; y14 = Y[14];  y15 = Y[15];
               y0 *= x0; y1 *= x0; y2 *= x0; y3 *= x0;
               y8 *= x0; y9 *= x0; y10 *= x0; y11 *= x0;
               y4 *= x0; y5 *= x0; y6 *= x0; y7 *= x0;
               y12 *= x0; y13 *= x0; y14 *= x0; y15 *= x0;
            #else
               y0 = *Y; y1 = Y[1]; y2 = Y[2];  y3 = Y[3];
               y8 = Y[8]; y9 = Y[9]; y10 = Y[10];  y11 = Y[11];
               y4 = Y[4]; y5 = Y[5]; y6 = Y[6];  y7 = Y[7];
               y12 = Y[12]; y13 = Y[13]; y14 = Y[14];  y15 = Y[15];
            #endif
            p0 = *A0; p1 = A1[1];
            p2 = A2[2]; p3 = A3[3];
            x0 = *X; x1 = X[1]; x2 = X[2]; x3 = X[3];
            x = X + 4;
            if (N4 != 4)
            {
               do
               {
                  y0  += x0 * p0; p0 = A0[incA];
                  y1  += x1 * p1; p1 = A1[incA+1];
                  y2  += x2 * p2; p2 = A2[incA+2];
                  y3  += x3 * p3; p3 = A3[incA+3];
                  y8  += x0 * A0[8];
                  y9  += x1 * A1[9];
                  y10 += x2 * A2[10];
                  y11 += x3 * A3[11];
                  y4  += x0 * A0[4];
                  y5  += x1 * A1[5];
                  y6  += x2 * A2[6];
                  y7  += x3 * A3[7];
                  y12 += x0 * A0[12];
                  y13 += x1 * A1[13];
                  y14 += x2 * A2[14];
                  y15 += x3 * A3[15];

                  y0  += x1 * *A1;
                  y1  += x0 * A0[1];
                  y2  += x0 * A0[2];
                  y3  += x0 * A0[3];
                  y8  += x1 * A1[8];
                  y9  += x0 * A0[9];
                  y10 += x0 * A0[10];
                  y11 += x0 * A0[11];
                  y4  += x1 * A1[4];
                  y5  += x0 * A0[5];
                  y6  += x0 * A0[6];
                  y7  += x0 * A0[7];
                  y12 += x1 * A1[12];
                  y13 += x0 * A0[13];
                  y14 += x0 * A0[14];
                  y15 += x0 * A0[15]; x0 = *x;

                  y0  += x2 * *A2;
                  y1  += x2 * A2[1]; A0 += incA;
                  y2  += x1 * A1[2];
                  y3  += x1 * A1[3];
                  y8  += x2 * A2[8];
                  y9  += x2 * A2[9];
                  y10 += x1 * A1[10];
                  y11 += x1 * A1[11];
                  y4  += x2 * A2[4];
                  y5  += x2 * A2[5];
                  y6  += x1 * A1[6];
                  y7  += x1 * A1[7];
                  y12 += x2 * A2[12];
                  y13 += x2 * A2[13];
                  y14 += x1 * A1[14];
                  y15 += x1 * A1[15]; x1 = x[1];

                  y0  += x3 * *A3;
                  y1  += x3 * A3[1];
                  y2  += x3 * A3[2]; A1 += incA;
                  y3  += x2 * A2[3];
                  y8  += x3 * A3[8];
                  y9  += x3 * A3[9];
                  y10 += x3 * A3[10];
                  y11 += x2 * A2[11];
                  y4  += x3 * A3[4];
                  y5  += x3 * A3[5];
                  y6  += x3 * A3[6];
                  y7  += x2 * A2[7];
                  y12 += x3 * A3[12];
                  y13 += x3 * A3[13];
                  y14 += x3 * A3[14]; x3 = x[3]; A3 += incA;
                  y15 += x2 * A2[15]; x2 = x[2]; x += 4; A2 += incA;
               }
               while (x != stX);
            }
            x -= 4;
            for (j=0; j != nr; j++, A0 += lda)
            {
               x0 = x[j];
               y0  += x0 * *A0;
               y1  += x0 * A0[1];
               y2  += x0 * A0[2];
               y3  += x0 * A0[3];
               y8  += x0 * A0[8];
               y9  += x0 * A0[9];
               y10 += x0 * A0[10];
               y11 += x0 * A0[11];
               y4  += x0 * A0[4];
               y5  += x0 * A0[5];
               y6  += x0 * A0[6];
               y7  += x0 * A0[7];
               y12 += x0 * A0[12];
               y13 += x0 * A0[13];
               y14 += x0 * A0[14];
               y15 += x0 * A0[15];
            }
            A0 += incAm0;
            *Y   = y0;
            Y[ 1] = y1 ;
            Y[ 2] = y2 ;
            Y[ 3] = y3 ;
            A1 += incAm;
            Y[ 8] = y8 ;
            Y[ 9] = y9 ;
            Y[10] = y10;
            Y[11] = y11;
            A2 += incAm;
            Y[ 4] = y4 ;
            Y[ 5] = y5 ;
            Y[ 6] = y6 ;
            Y[ 7] = y7 ;
            A3 += incAm;
            Y[12] = y12;
            Y[13] = y13;
            Y[14] = y14;
            Y[15] = y15;
            Y += 16;
         }
         while (Y != stY);
      }
      if (M-M16) gemvMlt8(M-M16, N, A0, lda, X, beta, Y);
   }
   else if (M) gemvNle4(M, N, A, lda, X, beta, Y);
}
void Mjoin(Mjoin(Mjoin(Mjoin(Mjoin(PATL,gemvN),NM),_x1),BNM),_y1)
   (const int M, const int N, const SCALAR alpha, const TYPE *A, const int lda,
    const TYPE *X, const int incX, const SCALAR beta, TYPE *Y, const int incY)
{
   gemv16x4(M, N, A, lda, X, beta, Y);
}
