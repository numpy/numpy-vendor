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

void Mjoin(Mjoin(Mjoin(Mjoin(Mjoin(PATL,gemvN),NM),_x1),BNM),_y1)
   (const int M, const int N, const SCALAR alpha, const TYPE *A, const int lda,
    const TYPE *X, const int incX, const SCALAR beta, TYPE *Y, const int incY)
{
   const int incA = lda<<1, incAm = 4 - ((N>>1)<<1)*lda;
   const int m4 = (M>>2)<<2;
   int n2, nr;
   register TYPE y0, y1, y2, y3, z0, z1, z2, z3, x0, x1, m0, m1, m2, m3;
   register TYPE a00, a10, a20, a30, a01, a11, a21, a31;
   const TYPE *x, *stX = X + ((N>>1)<<1)-2, *A0 = A, *A1 = A + lda;
   TYPE *stY = Y + m4;

   if (N > 4)
   {
      n2 = ((N-4)>>1)<<1;
      nr = N - n2;
      if (m4)
      {
         do
         {
            x = X + 2;
            #ifdef BETA0
               z0 = z1 = z2 = z3 = y0 = y1 = y2 = y3 = ATL_rzero;
            #else
               z0 = *Y;
               z1 = Y[1];
               z2 = Y[2];
               z3 = Y[3];
               #ifdef BETAX
                  y0 = beta;
                  z0 *= y0;
                  z1 *= y0;
                  z2 *= y0;
                  z3 *= y0;
               #endif
               y0 = y1 = y2 = y3 = ATL_rzero;
            #endif
            x0 = *X;
            x1 = X[1];
            a00 = *A0;
            a01 = *A1;
            a10 = A0[1];
            a11 = A1[1];
            a20 = A0[2];
            a21 = A1[2];
            a30 = A0[3];
            a31 = A1[3];
            A0 += incA;
            A1 += incA;
            m0 = x0 * a00;
            a00 = *A0;
            m1 = x1 * a01;
            a01 = *A1;
            m2 = x0 * a10;
            a10 = A0[1];
            m3 = x1 * a11;
            a11 = A1[1];
            if (n2)
            {
               do
               {
                  y0 += m0;
                  m0 = x0 * a20;
                  a20 = A0[2];
                  z0 += m1;
                  m1 = x1 * a21;
                  a21 = A1[2];
                  y1 += m2;
                  m2 = x0 * a30;
                  x0 = *x;
                  a30 = A0[3];
                  A0 += incA;
                  z1 += m3;
                  m3 = x1 * a31;
                  x1 = x[1];
                  a31 = A1[3];
                  x += 2;
                  A1 += incA;

                  y2 += m0;
                  m0 = x0 * a00;
                  a00 = *A0;
                  z2 += m1;
                  m1 = x1 * a01;
                  a01 = *A1;
                  y3 += m2;
                  m2 = x0 * a10;
                  a10 = A0[1];
                  z3 += m3;
                  m3 = x1 * a11;
                  a11 = A1[1];
               }
               while (x != stX);
            }
            if (nr == 4)
            {
               y0 += m0;
               m0 = x0 * a20;
               a20 = A0[2];
               z0 += m1;
               m1 = x1 * a21;
               a21 = A1[2];
               y1 += m2;
               m2 = x0 * a30;
               x0 = *x;
               a30 = A0[3];
               z1 += m3;
               m3 = x1 * a31;
               x1 = x[1];
               a31 = A1[3];

               y2 += m0;
               m0 = x0 * a00;
               z2 += m1;
               m1 = x1 * a01;
               y3 += m2;
               m2 = x0 * a10;
               z3 += m3;
               m3 = x1 * a11;

               y0 += m0;
               m0 = x0 * a20;
               z0 += m1;
               m1 = x1 * a21;
               y1 += m2;
               m2 = x0 * a30;
               z1 += m3;
               m3 = x1 * a31;

               y2 += m0;
               A0 += incA;
               z2 += m1;
               A1 += incA;
               y3 += m2;
               z3 += m3;
            }
            else /* nr == 5 */
            {
               y0 += m0;
               m0 = x0 * a20;
               a20 = A0[2];
               z0 += m1;
               m1 = x1 * a21;
               a21 = A1[2];
               y1 += m2;
               m2 = x0 * a30;
               x0 = *x;
               a30 = A0[3];
               A0 += incA;
               z1 += m3;
               m3 = x1 * a31;
               x1 = x[1];
               x += 2;
               a31 = A1[3];

               y2 += m0;
               m0 = x0 * a00;
               a00 = *A0;
               z2 += m1;
               m1 = x1 * a01;
               y3 += m2;
               m2 = x0 * a10;
               a10 = A0[1];
               z3 += m3;
               m3 = x1 * a11;

               y0 += m0;
               m0 = x0 * a20;
               a20 = A0[2];
               z0 += m1;
               m1 = x1 * a21;
               y1 += m2;
               m2 = x0 * a30;
               x0 = *x;
               a30 = A0[3];
               z1 += m3;
               m3 = x1 * a31;

               y2 += m0;
               m0 = x0 * a00;
               z2 += m1;
               m1 = x0 * a10;
               y3 += m2;
               m2 = x0 * a20;
               z3 += m3;
               m3 = x0 * a30;

               y0 += m0;
               A1 += incA;
               y1 += m1;
               y2 += m2;
               y3 += m3;
            }

            y0 += z0;
            A0 += incAm;
            y1 += z1;
            A1 += incAm;
            y2 += z2;
            y3 += z3;
            *Y = y0;
            Y[1] = y1;
            Y[2] = y2;
            Y[3] = y3;
            Y += 4;
         }
         while (Y != stY);
      }
      for (nr=M-m4; nr; nr--)
      {
         #ifdef BETA0
            y0 = Mjoin(PATL,dot)(N, A0, lda, X, 1);
         #else
            #if defined(BETAX)
               y0 = *Y * beta;
            #else
               y0 = *Y;
            #endif
            y0 += Mjoin(PATL,dot)(N, A0, lda, X, 1);
         #endif
         *Y++ = y0;
         A0++;
      }
   }
   else if (M) gemvNle4(M, N, A, lda, X, beta, Y);
}
