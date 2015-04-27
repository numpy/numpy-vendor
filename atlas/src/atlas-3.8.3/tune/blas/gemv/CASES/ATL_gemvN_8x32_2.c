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

#define Yass(y_) (y_) +=

#ifdef BETA1
static void gemv8x32(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
/*
 * Try to let compiler do all the work
 */
{
   const TYPE *A0=A, *A1=A+lda, *A2=A1+lda, *A3=A2+lda;
   const TYPE *A4=A3+lda, *A5=A4+lda, *A6=A5+lda, *A7=A6+lda;
   const TYPE *A8=A7+lda, *A9=A8+lda, *A10=A9+lda, *A11=A10+lda;
   const TYPE *A12=A11+lda, *A13=A12+lda, *A14=A13+lda, *A15=A14+lda;
   const TYPE *A16=A15+lda, *A17=A16+lda, *A18=A17+lda, *A19=A18+lda;
   const TYPE *A20=A19+lda, *A21=A20+lda, *A22=A21+lda, *A23=A22+lda;
   const TYPE *A24=A23+lda, *A25=A24+lda, *A26=A25+lda, *A27=A26+lda;
   const TYPE *A28=A27+lda, *A29=A28+lda, *A30=A29+lda, *A31=A30+lda;
   const TYPE x0=*X, x1=X[1], x2=X[2], x3=X[3];
   const TYPE x4=X[4], x5=X[5], x6=X[6], x7=X[7];
   const TYPE x8=X[8], x9=X[9], x10=X[10], x11=X[11];
   const TYPE x12=X[12], x13=X[13], x14=X[14], x15=X[15];
   const TYPE x16=X[16], x17=X[17], x18=X[18], x19=X[19];
   const TYPE x20=X[20], x21=X[21], x22=X[22], x23=X[23];
   const TYPE x24=X[24], x25=X[25], x26=X[26], x27=X[27];
   const TYPE x28=X[28], x29=X[29], x30=X[30], x31=X[31];
   TYPE *stY = Y + ((M>>3)<<3);

   ATL_assert(M >= 8 && (N == 32));
   do
   {
      Yass(*Y) *A0 * x0 + *A1 * x1 + *A2 * x2 + *A3 * x3 +
               *A4 * x4 + *A5 * x5 + *A6 * x6 + *A7 * x7 +
               *A8 * x8 + *A9 * x9 + *A10 * x10 + *A11 * x11 +
               *A12 * x12 + *A13 * x13 + *A14 * x14 + *A15 * x15 +
               *A16 * x16 + *A17 * x17 + *A18 * x18 + *A19 * x19 +
               *A20 * x20 + *A21 * x21 + *A22 * x22 + *A23 * x23 +
               *A24 * x24 + *A25 * x25 + *A26 * x26 + *A27 * x27 +
               *A28 * x28 + *A29 * x29 + *A30 * x30 + *A31 * x31;
      Yass(Y[1]) A0[1]*x0 + A1[1]*x1 + A2[1]*x2 + A3[1]*x3 +
                 A4[1]*x4 + A5[1]*x5 + A6[1]*x6 + A7[1]*x7 +
                 A8[1]*x8 + A9[1]*x9 + A10[1]*x10 + A11[1]*x11 +
                 A12[1]*x12 + A13[1]*x13 + A14[1]*x14 + A15[1]*x15 +
                 A16[1]*x16 + A17[1]*x17 + A18[1]*x18 + A19[1]*x19 +
                 A20[1]*x20 + A21[1]*x21 + A22[1]*x22 + A23[1]*x23 +
                 A24[1]*x24 + A25[1]*x25 + A26[1]*x26 + A27[1]*x27 +
                 A28[1]*x28 + A29[1]*x29 + A30[1]*x30 + A31[1]*x31;
      Yass(Y[2]) A0[2]*x0 + A1[2]*x1 + A2[2]*x2 + A3[2]*x3 +
                 A4[2]*x4 + A5[2]*x5 + A6[2]*x6 + A7[2]*x7 +
                 A8[2]*x8 + A9[2]*x9 + A10[2]*x10 + A11[2]*x11 +
                 A12[2]*x12 + A13[2]*x13 + A14[2]*x14 + A15[2]*x15 +
                 A16[2]*x16 + A17[2]*x17 + A18[2]*x18 + A19[2]*x19 +
                 A20[2]*x20 + A21[2]*x21 + A22[2]*x22 + A23[2]*x23 +
                 A24[2]*x24 + A25[2]*x25 + A26[2]*x26 + A27[2]*x27 +
                 A28[2]*x28 + A29[2]*x29 + A30[2]*x30 + A31[2]*x31;
      Yass(Y[3]) A0[3]*x0 + A1[3]*x1 + A2[3]*x2 + A3[3]*x3 +
                 A4[3]*x4 + A5[3]*x5 + A6[3]*x6 + A7[3]*x7 +
                 A8[3]*x8 + A9[3]*x9 + A10[3]*x10 + A11[3]*x11 +
                 A12[3]*x12 + A13[3]*x13 + A14[3]*x14 + A15[3]*x15 +
                 A16[3]*x16 + A17[3]*x17 + A18[3]*x18 + A19[3]*x19 +
                 A20[3]*x20 + A21[3]*x21 + A22[3]*x22 + A23[3]*x23 +
                 A24[3]*x24 + A25[3]*x25 + A26[3]*x26 + A27[3]*x27 +
                 A28[3]*x28 + A29[3]*x29 + A30[3]*x30 + A31[3]*x31;
      Yass(Y[4]) A0[4]*x0 + A1[4]*x1 + A2[4]*x2 + A3[4]*x3 +
                 A4[4]*x4 + A5[4]*x5 + A6[4]*x6 + A7[4]*x7 +
                 A8[4]*x8 + A9[4]*x9 + A10[4]*x10 + A11[4]*x11 +
                 A12[4]*x12 + A13[4]*x13 + A14[4]*x14 + A15[4]*x15 +
                 A16[4]*x16 + A17[4]*x17 + A18[4]*x18 + A19[4]*x19 +
                 A20[4]*x20 + A21[4]*x21 + A22[4]*x22 + A23[4]*x23 +
                 A24[4]*x24 + A25[4]*x25 + A26[4]*x26 + A27[4]*x27 +
                 A28[4]*x28 + A29[4]*x29 + A30[4]*x30 + A31[4]*x31;
      Yass(Y[5]) A0[5]*x0 + A1[5]*x1 + A2[5]*x2 + A3[5]*x3 +
                 A4[5]*x4 + A5[5]*x5 + A6[5]*x6 + A7[5]*x7 +
                 A8[5]*x8 + A9[5]*x9 + A10[5]*x10 + A11[5]*x11 +
                 A12[5]*x12 + A13[5]*x13 + A14[5]*x14 + A15[5]*x15 +
                 A16[5]*x16 + A17[5]*x17 + A18[5]*x18 + A19[5]*x19 +
                 A20[5]*x20 + A21[5]*x21 + A22[5]*x22 + A23[5]*x23 +
                 A24[5]*x24 + A25[5]*x25 + A26[5]*x26 + A27[5]*x27 +
                 A28[5]*x28 + A29[5]*x29 + A30[5]*x30 + A31[5]*x31;
      Yass(Y[6]) A0[6]*x0 + A1[6]*x1 + A2[6]*x2 + A3[6]*x3 +
                 A4[6]*x4 + A5[6]*x5 + A6[6]*x6 + A7[6]*x7 +
                 A8[6]*x8 + A9[6]*x9 + A10[6]*x10 + A11[6]*x11 +
                 A12[6]*x12 + A13[6]*x13 + A14[6]*x14 + A15[6]*x15 +
                 A16[6]*x16 + A17[6]*x17 + A18[6]*x18 + A19[6]*x19 +
                 A20[6]*x20 + A21[6]*x21 + A22[6]*x22 + A23[6]*x23 +
                 A24[6]*x24 + A25[6]*x25 + A26[6]*x26 + A27[6]*x27 +
                 A28[6]*x28 + A29[6]*x29 + A30[6]*x30 + A31[6]*x31;
      Yass(Y[7]) A0[7]*x0 + A1[7]*x1 + A2[7]*x2 + A3[7]*x3 +
                 A4[7]*x4 + A5[7]*x5 + A6[7]*x6 + A7[7]*x7 +
                 A8[7]*x8 + A9[7]*x9 + A10[7]*x10 + A11[7]*x11 +
                 A12[7]*x12 + A13[7]*x13 + A14[7]*x14 + A15[7]*x15 +
                 A16[7]*x16 + A17[7]*x17 + A18[7]*x18 + A19[7]*x19 +
                 A20[7]*x20 + A21[7]*x21 + A22[7]*x22 + A23[7]*x23 +
                 A24[7]*x24 + A25[7]*x25 + A26[7]*x26 + A27[7]*x27 +
                 A28[7]*x28 + A29[7]*x29 + A30[7]*x30 + A31[7]*x31;
      A0 += 8;
      A1 += 8;
      A2 += 8;
      A3 += 8;
      A4 += 8;
      A5 += 8;
      A6 += 8;
      A7 += 8;
      A8 += 8;
      A9 += 8;
      A10 += 8;
      A11 += 8;
      A12 += 8;
      A13 += 8;
      A14 += 8;
      A15 += 8;
      A16 += 8;
      A17 += 8;
      A18 += 8;
      A19 += 8;
      A20 += 8;
      A21 += 8;
      A22 += 8;
      A23 += 8;
      A24 += 8;
      A25 += 8;
      A26 += 8;
      A27 += 8;
      A28 += 8;
      A29 += 8;
      A30 += 8;
      A31 += 8;
      Y += 8;
   }
   while (Y != stY);
}
static void gemv1x32(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
/*
 * Try to let compiler do all the work
 */
{
   const TYPE *A0=A, *A1=A+lda, *A2=A1+lda, *A3=A2+lda;
   const TYPE *A4=A3+lda, *A5=A4+lda, *A6=A5+lda, *A7=A6+lda;
   const TYPE *A8=A7+lda, *A9=A8+lda, *A10=A9+lda, *A11=A10+lda;
   const TYPE *A12=A11+lda, *A13=A12+lda, *A14=A13+lda, *A15=A14+lda;
   const TYPE *A16=A15+lda, *A17=A16+lda, *A18=A17+lda, *A19=A18+lda;
   const TYPE *A20=A19+lda, *A21=A20+lda, *A22=A21+lda, *A23=A22+lda;
   const TYPE *A24=A23+lda, *A25=A24+lda, *A26=A25+lda, *A27=A26+lda;
   const TYPE *A28=A27+lda, *A29=A28+lda, *A30=A29+lda, *A31=A30+lda;
   const TYPE x0=*X, x1=X[1], x2=X[2], x3=X[3];
   const TYPE x4=X[4], x5=X[5], x6=X[6], x7=X[7];
   const TYPE x8=X[8], x9=X[9], x10=X[10], x11=X[11];
   const TYPE x12=X[12], x13=X[13], x14=X[14], x15=X[15];
   const TYPE x16=X[16], x17=X[17], x18=X[18], x19=X[19];
   const TYPE x20=X[20], x21=X[21], x22=X[22], x23=X[23];
   const TYPE x24=X[24], x25=X[25], x26=X[26], x27=X[27];
   const TYPE x28=X[28], x29=X[29], x30=X[30], x31=X[31];
   TYPE *stY = Y + M;

   ATL_assert(N == 32);
   do
   {
      Yass(*Y) *A0 * x0 + *A1 * x1 + *A2 * x2 + *A3 * x3 +
               *A4 * x4 + *A5 * x5 + *A6 * x6 + *A7 * x7 +
               *A8 * x8 + *A9 * x9 + *A10 * x10 + *A11 * x11 +
               *A12 * x12 + *A13 * x13 + *A14 * x14 + *A15 * x15 +
               *A16 * x16 + *A17 * x17 + *A18 * x18 + *A19 * x19 +
               *A20 * x20 + *A21 * x21 + *A22 * x22 + *A23 * x23 +
               *A24 * x24 + *A25 * x25 + *A26 * x26 + *A27 * x27 +
               *A28 * x28 + *A29 * x29 + *A30 * x30 + *A31 * x31;
      Y++;
      A0++;
      A1++;
      A2++;
      A3++;
      A4++;
      A5++;
      A6++;
      A7++;
      A8++;
      A9++;
      A10++;
      A11++;
      A12++;
      A13++;
      A14++;
      A15++;
      A16++;
      A17++;
      A18++;
      A19++;
      A20++;
      A21++;
      A22++;
      A23++;
      A24++;
      A25++;
      A26++;
      A27++;
      A28++;
      A29++;
      A30++;
      A31++;
   }
   while (Y != stY);
}
static void gemv8x16(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
/*
 * Try to let compiler do all the work
 */
{
   const TYPE *A0=A, *A1=A+lda, *A2=A1+lda, *A3=A2+lda;
   const TYPE *A4=A3+lda, *A5=A4+lda, *A6=A5+lda, *A7=A6+lda;
   const TYPE *A8=A7+lda, *A9=A8+lda, *A10=A9+lda, *A11=A10+lda;
   const TYPE *A12=A11+lda, *A13=A12+lda, *A14=A13+lda, *A15=A14+lda;
   const TYPE x0=*X, x1=X[1], x2=X[2], x3=X[3];
   const TYPE x4=X[4], x5=X[5], x6=X[6], x7=X[7];
   const TYPE x8=X[8], x9=X[9], x10=X[10], x11=X[11];
   const TYPE x12=X[12], x13=X[13], x14=X[14], x15=X[15];
   TYPE *stY = Y + ((M>>3)<<3);

   ATL_assert(M >= 8 && (N == 16));
   do
   {
      Yass(*Y) *A0 * x0 + *A1 * x1 + *A2 * x2 + *A3 * x3 +
               *A4 * x4 + *A5 * x5 + *A6 * x6 + *A7 * x7 +
               *A8 * x8 + *A9 * x9 + *A10 * x10 + *A11 * x11 +
               *A12 * x12 + *A13 * x13 + *A14 * x14 + *A15 * x15;
      Yass(Y[1]) A0[1]*x0 + A1[1]*x1 + A2[1]*x2 + A3[1]*x3 +
                 A4[1]*x4 + A5[1]*x5 + A6[1]*x6 + A7[1]*x7 +
                 A8[1]*x8 + A9[1]*x9 + A10[1]*x10 + A11[1]*x11 +
                 A12[1]*x12 + A13[1]*x13 + A14[1]*x14 + A15[1]*x15;
      Yass(Y[2]) A0[2]*x0 + A1[2]*x1 + A2[2]*x2 + A3[2]*x3 +
                 A4[2]*x4 + A5[2]*x5 + A6[2]*x6 + A7[2]*x7 +
                 A8[2]*x8 + A9[2]*x9 + A10[2]*x10 + A11[2]*x11 +
                 A12[2]*x12 + A13[2]*x13 + A14[2]*x14 + A15[2]*x15;
      Yass(Y[3]) A0[3]*x0 + A1[3]*x1 + A2[3]*x2 + A3[3]*x3 +
                 A4[3]*x4 + A5[3]*x5 + A6[3]*x6 + A7[3]*x7 +
                 A8[3]*x8 + A9[3]*x9 + A10[3]*x10 + A11[3]*x11 +
                 A12[3]*x12 + A13[3]*x13 + A14[3]*x14 + A15[3]*x15;
      Yass(Y[4]) A0[4]*x0 + A1[4]*x1 + A2[4]*x2 + A3[4]*x3 +
                 A4[4]*x4 + A5[4]*x5 + A6[4]*x6 + A7[4]*x7 +
                 A8[4]*x8 + A9[4]*x9 + A10[4]*x10 + A11[4]*x11 +
                 A12[4]*x12 + A13[4]*x13 + A14[4]*x14 + A15[4]*x15;
      Yass(Y[5]) A0[5]*x0 + A1[5]*x1 + A2[5]*x2 + A3[5]*x3 +
                 A4[5]*x4 + A5[5]*x5 + A6[5]*x6 + A7[5]*x7 +
                 A8[5]*x8 + A9[5]*x9 + A10[5]*x10 + A11[5]*x11 +
                 A12[5]*x12 + A13[5]*x13 + A14[5]*x14 + A15[5]*x15;
      Yass(Y[6]) A0[6]*x0 + A1[6]*x1 + A2[6]*x2 + A3[6]*x3 +
                 A4[6]*x4 + A5[6]*x5 + A6[6]*x6 + A7[6]*x7 +
                 A8[6]*x8 + A9[6]*x9 + A10[6]*x10 + A11[6]*x11 +
                 A12[6]*x12 + A13[6]*x13 + A14[6]*x14 + A15[6]*x15;
      Yass(Y[7]) A0[7]*x0 + A1[7]*x1 + A2[7]*x2 + A3[7]*x3 +
                 A4[7]*x4 + A5[7]*x5 + A6[7]*x6 + A7[7]*x7 +
                 A8[7]*x8 + A9[7]*x9 + A10[7]*x10 + A11[7]*x11 +
                 A12[7]*x12 + A13[7]*x13 + A14[7]*x14 + A15[7]*x15;
      A0 += 8;
      A1 += 8;
      A2 += 8;
      A3 += 8;
      A4 += 8;
      A5 += 8;
      A6 += 8;
      A7 += 8;
      A8 += 8;
      A9 += 8;
      A10 += 8;
      A11 += 8;
      A12 += 8;
      A13 += 8;
      A14 += 8;
      A15 += 8;
      Y += 8;
   }
   while (Y != stY);
}
static void gemv8x8(const int M, const int N, const TYPE *A, const int lda,
                    const TYPE *X, const SCALAR beta, TYPE *Y)
/*
 * Try to let compiler do all the work
 */
{
   const TYPE *A0=A, *A1=A+lda, *A2=A1+lda, *A3=A2+lda;
   const TYPE *A4=A3+lda, *A5=A4+lda, *A6=A5+lda, *A7=A6+lda;
   const TYPE x0=*X, x1=X[1], x2=X[2], x3=X[3];
   const TYPE x4=X[4], x5=X[5], x6=X[6], x7=X[7];
   TYPE *stY = Y + ((M>>3)<<3);

   ATL_assert(M >= 8 && (N == 8));
   do
   {
      Yass(*Y) *A0 * x0 + *A1 * x1 + *A2 * x2 + *A3 * x3 +
               *A4 * x4 + *A5 * x5 + *A6 * x6 + *A7 * x7;
      Yass(Y[1]) A0[1]*x0 + A1[1]*x1 + A2[1]*x2 + A3[1]*x3 +
                 A4[1]*x4 + A5[1]*x5 + A6[1]*x6 + A7[1]*x7;
      Yass(Y[2]) A0[2]*x0 + A1[2]*x1 + A2[2]*x2 + A3[2]*x3 +
                 A4[2]*x4 + A5[2]*x5 + A6[2]*x6 + A7[2]*x7;
      Yass(Y[3]) A0[3]*x0 + A1[3]*x1 + A2[3]*x2 + A3[3]*x3 +
                 A4[3]*x4 + A5[3]*x5 + A6[3]*x6 + A7[3]*x7;
      Yass(Y[4]) A0[4]*x0 + A1[4]*x1 + A2[4]*x2 + A3[4]*x3 +
                 A4[4]*x4 + A5[4]*x5 + A6[4]*x6 + A7[4]*x7;
      Yass(Y[5]) A0[5]*x0 + A1[5]*x1 + A2[5]*x2 + A3[5]*x3 +
                 A4[5]*x4 + A5[5]*x5 + A6[5]*x6 + A7[5]*x7;
      Yass(Y[6]) A0[6]*x0 + A1[6]*x1 + A2[6]*x2 + A3[6]*x3 +
                 A4[6]*x4 + A5[6]*x5 + A6[6]*x6 + A7[6]*x7;
      Yass(Y[7]) A0[7]*x0 + A1[7]*x1 + A2[7]*x2 + A3[7]*x3 +
                 A4[7]*x4 + A5[7]*x5 + A6[7]*x6 + A7[7]*x7;
      A0 += 8;
      A1 += 8;
      A2 += 8;
      A3 += 8;
      A4 += 8;
      A5 += 8;
      A6 += 8;
      A7 += 8;
      Y += 8;
   }
   while (Y != stY);
}

static void gemv8x4(const int M, const int N, const TYPE *A, const int lda,
                    const TYPE *X, const SCALAR beta, TYPE *Y)
/*
 * Try to let compiler do all the work
 */
{
   const TYPE *A0=A, *A1=A+lda, *A2=A1+lda, *A3=A2+lda;
   const TYPE x0=*X, x1=X[1], x2=X[2], x3=X[3];
   TYPE *stY = Y + ((M>>3)<<3);

   ATL_assert(M >= 8 && (N == 4));
   do
   {
      Yass(*Y) *A0 * x0 + *A1 * x1 + *A2 * x2 + *A3 * x3;
      Yass(Y[1]) A0[1]*x0 + A1[1]*x1 + A2[1]*x2 + A3[1]*x3;
      Yass(Y[2]) A0[2]*x0 + A1[2]*x1 + A2[2]*x2 + A3[2]*x3;
      Yass(Y[3]) A0[3]*x0 + A1[3]*x1 + A2[3]*x2 + A3[3]*x3;
      Yass(Y[4]) A0[4]*x0 + A1[4]*x1 + A2[4]*x2 + A3[4]*x3;
      Yass(Y[5]) A0[5]*x0 + A1[5]*x1 + A2[5]*x2 + A3[5]*x3;
      Yass(Y[6]) A0[6]*x0 + A1[6]*x1 + A2[6]*x2 + A3[6]*x3;
      Yass(Y[7]) A0[7]*x0 + A1[7]*x1 + A2[7]*x2 + A3[7]*x3;
      A0 += 8;
      A1 += 8;
      A2 += 8;
      A3 += 8;
      Y += 8;
   }
   while (Y != stY);
}

static void gemv1x1(const int M, const int N, const TYPE *A, const int lda,
                    const TYPE *X, const SCALAR beta, TYPE *Y)
{
   int i;
   for (i=0; i < M; i++) Y[i] += Mjoin(PATL,dot)(N, X, 1, A+i, lda);
}
#endif

static void gemv(const int M, const int N, const TYPE *A, const int lda,
                 const TYPE *X, const SCALAR beta, TYPE *Y)
{
#ifdef BETA0
   Mjoin(PATL,zero)(M, Y, 1);
   Mjoin(PATL,gemvN_a1_x1_b1_y1)(M, N, A, lda, X, beta, Y);
#elif defined(BETAX)
   Mjoin(PATL,scal)(M, beta, Y, 1);
   Mjoin(PATL,gemvN_a1_x1_b1_y1)(M, N, A, lda, X, beta, Y);
#else
   int i, n;
   const int n32 = N>>5, M8 = (M>>3)<<3, mr = M-M8, incA = lda<<5;

   #ifndef BETA1
      ATL_assert(beta == ATL_rone);
   #endif
   for (i=n32; i; i--)
   {
      if (M8) gemv8x32(M8, 32, A, lda, X, beta, Y);
      if (mr) gemv1x32(mr, 32, A+M8, lda, X, beta, Y+M8);
      A += incA;
      X += 32;
   }
   n = N - (n32<<5);
   if (n >= 16)
   {
      n -= 16;
      if (M8) gemv8x16(M8, 16, A, lda, X, beta, Y);
      if (mr) gemv1x1(mr, 16, A+M8, lda, X, beta, Y+M8);
      A += incA>>1;
      X += 16;
   }
   if (n >= 8)
   {
      n -= 8;
      if (M8) gemv8x8(M8, 8, A, lda, X, beta, Y);
      if (mr) gemv1x1(mr, 8, A+M8, lda, X, beta, Y+M8);
      A += incA>>2;
      X += 8;
   }
   if (n >= 4)
   {
      n -= 4;
      if (M8) gemv8x4(M8, 4, A, lda, X, beta, Y);
      if (mr) gemv1x1(mr, 4, A+M8, lda, X, beta, Y+M8);
      A += incA>>3;
      X += 4;
   }
   if (n) gemv1x1(M, n, A, lda, X, beta, Y);
#endif
}

void Mjoin(Mjoin(Mjoin(Mjoin(Mjoin(PATL,gemvN),NM),_x1),BNM),_y1)
   (const int M, const int N, const SCALAR alpha, const TYPE *A, const int lda,
    const TYPE *X, const int incX, const SCALAR beta, TYPE *Y, const int incY)
{
   gemv(M, N, A, lda, X, beta, Y);
}
