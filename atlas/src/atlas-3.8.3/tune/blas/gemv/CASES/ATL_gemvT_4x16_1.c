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
static void gemvT4x16(const int M, const int N, const TYPE *A, const int lda,
                      const TYPE *X, const SCALAR beta, TYPE *Y)
{
   int j;
   const int M4 = (M>>2)<<2, N16 = (N>>4)<<4, incAm = (lda<<2) - N;
   const int nr = ( N16 ? N-N16+16 : N );
   const TYPE *stX = X + N16-16, *x;
   const TYPE *A0 = A, *A1 = A + lda, *A2 = A1 + lda, *A3 = A2 + lda;
   TYPE *stY = Y + M4;
   register TYPE a00, a01, a02, a03, a10, a11, a12, a13;
   register TYPE b00, b01, b02, b03, b10, b11, b12, b13;
   register TYPE y0, y1, y2, y3, z0, z1, z2, z3;
   register TYPE x0, x1, x2, x3;
   #ifdef ATL_AltiVec
      int cwrd = ATL_MulBySize(N)>>4;
      if (cwrd >= 64) cwrd = ATL_GetCtrl(512, (cwrd+31)>>5, 0);
      else cwrd = ATL_GetCtrl(64, (cwrd+3)>>2, 4);
   #endif


   if (N >= 10)
   {
      if (M4)
      {
         do
         {
            #ifdef ATL_AltiVec
               ATL_pfavR(A0, cwrd, 0);
               ATL_pfavR(A1, cwrd, 1);
               ATL_pfavR(A2, cwrd, 2);
               ATL_pfavR(A3, cwrd, 3);
            #endif

            a00 = *A0;   a01 = *A1;   a02 = *A2;   a03 = *A3;
            a10 = A0[1]; a11 = A1[1]; a12 = A2[1]; a13 = A3[1];
            b00 = A0[8]; b01 = A1[8]; b02 = A2[8]; b03 = A3[8];
            b10 = A0[9]; b11 = A1[9]; b12 = A2[9]; b13 = A3[9];
            #ifdef BETA0
               y0 = y1 = y2 = y3 = z0 = z1 = z2 = z3 = ATL_rzero;
            #elif defined(BETAX)
               y0 = beta;
               z0 = *Y; z1 = Y[1]; z2 = Y[2];  z3 = Y[3];
               z0 *= y0; z1 *= y0; z2 *= y0; z3 *= y0;
               y0 = y1 = y2 = y3 = ATL_rzero;
            #else
               z0 = *Y; z1 = Y[1]; z2 = Y[2];  z3 = Y[3];
               y0 = y1 = y2 = y3 = ATL_rzero;
            #endif
            x0 = *X; x1 = X[1]; x2 = X[8]; x3 = X[9]; x = X;
            if (N16 > 16)
            {
               do
               {
                  y0  += a00 * x0; a00 = A0[2];
                  y1  += a01 * x0; a01 = A1[2];
                  y2  += a02 * x0; a02 = A2[2];
                  y3  += a03 * x0; a03 = A3[2]; x0 = x[2];
                  z0  += b00 * x2; b00 = A0[10];
                  z1  += b01 * x2; b01 = A1[10];
                  z2  += b02 * x2; b02 = A2[10];
                  z3  += b03 * x2; b03 = A3[10]; x2 = x[10];

                  y0  += a10 * x1; a10 = A0[3];
                  y1  += a11 * x1; a11 = A1[3];
                  y2  += a12 * x1; a12 = A2[3];
                  y3  += a13 * x1; a13 = A3[3]; x1 = x[3];
                  z0  += b10 * x3; b10 = A0[11];
                  z1  += b11 * x3; b11 = A1[11];
                  z2  += b12 * x3; b12 = A2[11];
                  z3  += b13 * x3; b13 = A3[11]; x3 = x[11];

                  y0  += a00 * x0; a00 = A0[4];
                  y1  += a01 * x0; a01 = A1[4];
                  y2  += a02 * x0; a02 = A2[4];
                  y3  += a03 * x0; a03 = A3[4]; x0 = x[4];
                  z0  += b00 * x2; b00 = A0[12];
                  z1  += b01 * x2; b01 = A1[12];
                  z2  += b02 * x2; b02 = A2[12];
                  z3  += b03 * x2; b03 = A3[12]; x2 = x[12];

                  y0  += a10 * x1; a10 = A0[5];
                  y1  += a11 * x1; a11 = A1[5];
                  y2  += a12 * x1; a12 = A2[5];
                  y3  += a13 * x1; a13 = A3[5]; x1 = x[5];
                  z0  += b10 * x3; b10 = A0[13];
                  z1  += b11 * x3; b11 = A1[13];
                  z2  += b12 * x3; b12 = A2[13];
                  z3  += b13 * x3; b13 = A3[13]; x3 = x[13];

                  y0  += a00 * x0; a00 = A0[6];
                  y1  += a01 * x0; a01 = A1[6];
                  y2  += a02 * x0; a02 = A2[6];
                  y3  += a03 * x0; a03 = A3[6]; x0 = x[6];
                  z0  += b00 * x2; b00 = A0[14];
                  z1  += b01 * x2; b01 = A1[14];
                  z2  += b02 * x2; b02 = A2[14];
                  z3  += b03 * x2; b03 = A3[14]; x2 = x[14];

                  y0  += a10 * x1; a10 = A0[7];
                  y1  += a11 * x1; a11 = A1[7];
                  y2  += a12 * x1; a12 = A2[7];
                  y3  += a13 * x1; a13 = A3[7]; x1 = x[7];
                  z0  += b10 * x3; b10 = A0[15]; A0 += 16;
                  z1  += b11 * x3; b11 = A1[15]; A1 += 16;
                  z2  += b12 * x3; b12 = A2[15]; A2 += 16;
                  z3  += b13 * x3; b13 = A3[15]; x3 = x[15]; x += 16; A3 += 16;

                  y0  += a00 * x0; a00 = *A0;
                  y1  += a01 * x0; a01 = *A1;
                  y2  += a02 * x0; a02 = *A2;
                  y3  += a03 * x0; a03 = *A3; x0 = *x;
                  z0  += b00 * x2; b00 = A0[8];
                  z1  += b01 * x2; b01 = A1[8];
                  z2  += b02 * x2; b02 = A2[8];
                  z3  += b03 * x2; b03 = A3[8]; x2 = x[8];

                  y0  += a10 * x1; a10 = A0[1];
                  y1  += a11 * x1; a11 = A1[1];
                  y2  += a12 * x1; a12 = A2[1];
                  y3  += a13 * x1; a13 = A3[1]; x1 = x[1];
                  z0  += b10 * x3; b10 = A0[9];
                  z1  += b11 * x3; b11 = A1[9];
                  z2  += b12 * x3; b12 = A2[9];
                  z3  += b13 * x3; b13 = A3[9]; x3 = x[9];
               }
               while (x != stX);
            }
            for (j=(nr>>1); j; j--, A0 += 2, A1 += 2, A2 += 2, A3 += 2, x += 2)
            {
               x0 = *x; x1 = x[1];
               y0 += *A0 * x0;
               y1 += *A1 * x0;
               y2 += *A2 * x0;
               y3 += *A3 * x0;
               z0 += A0[1] * x1;
               z1 += A1[1] * x1;
               z2 += A2[1] * x1;
               z3 += A3[1] * x1;
            }
            if (nr != (nr>>1)<<1)
            {
               x0 = *x;
               z0 += *A0++ * x0;
               z1 += *A1++ * x0;
               z2 += *A2++ * x0;
               z3 += *A3++ * x0;
            }
            A0 += incAm;
            y0 += z0;
            A1 += incAm;
            y1 += z1;
            A2 += incAm;
            y2 += z2;
            A3 += incAm;
            y3 += z3;
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
   gemvT4x16(M, N, A, lda, X, beta, Y);
}
