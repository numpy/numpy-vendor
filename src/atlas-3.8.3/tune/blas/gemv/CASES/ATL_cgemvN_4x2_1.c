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
#include "atlas_level2.h"
#include "atlas_level1.h"
#include "atlas_prefetch.h"

#ifdef Conj_
   #define PEQ -=
   #define MEQ +=
#else
   #define PEQ +=
   #define MEQ -=
#endif

#ifdef Conj_
   #define gemv0 Mjoin(Mjoin(Mjoin(Mjoin(gemvNc,NM),_x1),BNM),_y1)
#else
   #define gemv0 Mjoin(Mjoin(Mjoin(Mjoin(gemvN,NM),_x1),BNM),_y1)
#endif

static void gemvMlt4(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
{
   int i;
   #if defined(BETAX) || defined(BETA1) || defined(BETAXI0)
      register TYPE yr1, yi1;
   #endif
   #ifdef BETAXI0
      const register TYPE rbeta = *beta;
   #endif
   #ifdef BETAX
      register TYPE yr0, yi0;
      const register TYPE rbeta = *beta, ibeta = beta[1];
   #endif

   for(i=M; i; i--, A += 2, Y += 2)
   {
      #ifdef BETA0
         #ifdef Conj_
            Mjoin(PATL,dotc_sub)(N, A, lda, X, 1, Y);
         #else
            Mjoin(PATL,dotu_sub)(N, A, lda, X, 1, Y);
         #endif
      #else
         #ifdef BETAX
            yr0 = *Y;
            yi0 = Y[1];
            yr1 = yr0 * rbeta - yi0 * ibeta;
            yi1 = yr0 * ibeta + yi0 * rbeta;
         #else
            yr1 = *Y;
            yi1 = Y[1];
            #ifdef BETAXI0
               yr1 *= rbeta;
               yi1 *= rbeta;
            #endif
         #endif
         #ifdef Conj_
            Mjoin(PATL,dotc_sub)(N, A, lda, X, 1, Y);
         #else
            Mjoin(PATL,dotu_sub)(N, A, lda, X, 1, Y);
         #endif
         *Y += yr1;
         Y[1] += yi1;
      #endif
   }
}

static void gemvNlt2(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
{
   if (N)
   {
      #if defined(BETA0)
         #ifdef Conj_
            Mjoin(PATL,moveConj)(M, X, A, 1, Y, 1);
         #else
            Mjoin(PATL,cpsc)(M, X, A, 1, Y, 1);
         #endif
      #elif defined(BETAX) || defined(BETAXI0)
         #ifdef Conj_
            Mjoin(PATL,axpbyConj)(M, X, A, 1, beta, Y, 1);
         #else
            Mjoin(PATL,axpby)(M, X, A, 1, beta, Y, 1);
         #endif
      #else
         #ifdef Conj_
            Mjoin(PATL,axpyConj)(M, X, A, 1, Y, 1);
         #else
            Mjoin(PATL,axpy)(M, X, A, 1, Y, 1);
         #endif
      #endif
   }
}

static void gemvN4x2(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
{
   const int M4 = (M>>2)<<2, N2 = (N>>1)<<1, nr = N-N2;
   const int incA = (lda<<2), incAm = (4 - N2*lda)<<1;
   const TYPE *A0 = A, *A1 = A + (lda<<1);
   const TYPE *x, *stX = X + (N2<<1);
   TYPE *stY = Y + (M4<<1);
   register TYPE ar00, ar10, ar20, ar30, ai00, ai10, ai20, ai30;
   register TYPE ar01, ar11, ar21, ar31, ai01, ai11, ai21, ai31;
   register TYPE yr0, yr1, yr2, yr3, yi0, yi1, yi2, yi3;
   register TYPE xr0, xi0, xr1, xi1, p0, p1;
   #ifdef ATL_AltiVec
   const int cwrd=ATL_GetCtrl(ATL_MulBySize(lda), Mmin(N,255), ATL_sizeof>>1);
   #endif


   if (N2)
   {
      if (M4)
      {
         do
         {
            #ifdef ATL_AltiVec
               ATL_pfavR(A0, cwrd, 0);
            #endif
            #ifdef BETA0
               yr0 = yi0 = yr1 = yi1 = yr2 = yi2 = yr3 = yi3 = ATL_rzero;
            #elif defined(BETAX)
               xr0 = *beta; xi0 = beta[1];
               ar00 = *Y;   ai00 = Y[1]; ar10 = Y[2]; ai10 = Y[3];
               ar20 = Y[4]; ai20 = Y[5]; ar30 = Y[6]; ai30 = Y[7];
               yr0 = ar00 * xr0 - ai00 * xi0;
               yi0 = ar00 * xi0 + ai00 * xr0;
               yr1 = ar10 * xr0 - ai10 * xi0;
               yi1 = ar10 * xi0 + ai10 * xr0;
               yr2 = ar20 * xr0 - ai20 * xi0;
               yi2 = ar20 * xi0 + ai20 * xr0;
               yr3 = ar30 * xr0 - ai30 * xi0;
               yi3 = ar30 * xi0 + ai30 * xr0;
            #else
               yr0 = *Y;   yi0 = Y[1]; yr1 = Y[2]; yi1 = Y[3];
               yr2 = Y[4]; yi2 = Y[5]; yr3 = Y[6]; yi3 = Y[7];
               #ifdef BETAXI0
                  xr0 = *beta;
                  yr0 *= xr0; yi0 *= xr0;
                  yr1 *= xr0; yi1 *= xr0;
                  yr2 *= xr0; yi2 *= xr0;
                  yr3 *= xr0; yi3 *= xr0;
               #endif
            #endif
            xr0 = *X; xi0 = X[1]; xr1 = X[2]; xi1 = X[3]; x = X + 4;
            ar00 = *A0;   ar01 = *A1;   ai00 = A0[1]; ai01 = A1[1];
            ar10 = A0[2]; ar11 = A1[2]; ai10 = A0[3]; ai11 = A1[3];
            ar20 = A0[4]; ar21 = A1[4]; ai20 = A0[5]; ai21 = A1[5];
            ar30 = A0[6]; ar31 = A1[6]; ai30 = A0[7]; ai31 = A1[7];
            A0 += incA;  A1 += incA;
            if (N2 != 2)
            {
               do
               {
                  yr0 += ar00 * xr0;
                  yi0 += ar01 * xi1;
                  yr1 += ar10 * xr0; p0 = A1[2];
                  yi1 += ar11 * xi1;
                  yr2 += ar20 * xr0; p1 = A0[2];
                  yi2 += ar21 * xi1;
                  yr3 += ar30 * xr0;
                  yi3 += ar31 * xi1;

                  yr0 += ar01 * xr1; ar01 = *A1;
                  yi0 += ar00 * xi0; ar00 = *A0;
                  yr1 += ar11 * xr1; ar11 = p0;
                  yi1 += ar10 * xi0; ar10 = p1;
                  yr2 += ar21 * xr1; ar21 = A1[4];
                  yi2 += ar20 * xi0; ar20 = A0[4];
                  yr3 += ar31 * xr1; ar31 = A1[6];
                  yi3 += ar30 * xi0; ar30 = A0[6];

                  yr0 MEQ ai00 * xi0;
                  yi0 PEQ ai01 * xr1;
                  yr1 MEQ ai10 * xi0; p0 = *x;
                  yi1 PEQ ai11 * xr1;
                  yr2 MEQ ai20 * xi0; p1 = A1[3];
                  yi2 PEQ ai21 * xr1;
                  yr3 MEQ ai30 * xi0; xi0 = x[1];
                  yi3 PEQ ai31 * xr1; xr1 = x[2];

                  yr0 MEQ ai01 * xi1; ai01 = A1[1];
                  yi0 PEQ ai00 * xr0; ai00 = A0[1];
                  yr1 MEQ ai11 * xi1; ai11 = p1;
                  yi1 PEQ ai10 * xr0; ai10 = A0[3];
                  yr2 MEQ ai21 * xi1; ai21 = A1[5];
                  yi2 PEQ ai20 * xr0; ai20 = A0[5];
                  yr3 MEQ ai31 * xi1; ai31 = A1[7]; xi1 = x[3]; A1 += incA;
                  yi3 PEQ ai30 * xr0; xr0 = p0; ai30 = A0[7]; A0 += incA;
                  x += 4;
               }
               while (x != stX);
            }
            if (!nr)  /* 2 cols left to do, all fetching done */
            {
               yr0 += ar00 * xr0;
               yi0 += ar01 * xi1;
               yr1 += ar10 * xr0;
               yi1 += ar11 * xi1;
               yr2 += ar20 * xr0;
               yi2 += ar21 * xi1;
               yr3 += ar30 * xr0;
               yi3 += ar31 * xi1;

               yr0 += ar01 * xr1;
               yi0 += ar00 * xi0;
               yr1 += ar11 * xr1;
               yi1 += ar10 * xi0;
               yr2 += ar21 * xr1;
               yi2 += ar20 * xi0;
               yr3 += ar31 * xr1;
               yi3 += ar30 * xi0;

               yr0 MEQ ai00 * xi0;
               yi0 PEQ ai01 * xr1;
               yr1 MEQ ai10 * xi0;
               yi1 PEQ ai11 * xr1;
               yr2 MEQ ai20 * xi0;
               yi2 PEQ ai21 * xr1;
               yr3 MEQ ai30 * xi0;
               yi3 PEQ ai31 * xr1;

               yr0 MEQ ai01 * xi1;
               yi0 PEQ ai00 * xr0;
               yr1 MEQ ai11 * xi1;
               yi1 PEQ ai10 * xr0;
               yr2 MEQ ai21 * xi1;
               yi2 PEQ ai20 * xr0;
               yr3 MEQ ai31 * xi1;
               yi3 PEQ ai30 * xr0;
            }
            else  /* 3 cols left, 2 already fetched */
            {
               yr0 += ar00 * xr0;
               yi0 += ar01 * xi1;
               yr1 += ar10 * xr0; p0 = A0[2];
               yi1 += ar11 * xi1;
               yr2 += ar20 * xr0; p1 = A0[6];
               yi2 += ar21 * xi1;
               yr3 += ar30 * xr0;
               yi3 += ar31 * xi1;

               yr0 += ar01 * xr1;
               yi0 += ar00 * xi0; ar00 = *A0;
               yr1 += ar11 * xr1;
               yi1 += ar10 * xi0; ar10 = p0;
               yr2 += ar21 * xr1;
               yi2 += ar20 * xi0; ar20 = A0[4];
               yr3 += ar31 * xr1;
               yi3 += ar30 * xi0; ar30 = p1;

               yr0 MEQ ai00 * xi0;
               yi0 PEQ ai01 * xr1;
               yr1 MEQ ai10 * xi0; p0 = *x;
               yi1 PEQ ai11 * xr1;
               yr2 MEQ ai20 * xi0;
               yi2 PEQ ai21 * xr1;
               yr3 MEQ ai30 * xi0; xi0 = x[1];
               yi3 PEQ ai31 * xr1;

               yr0 MEQ ai01 * xi1;
               yi0 PEQ ai00 * xr0; ai00 = A0[1];
               yr1 MEQ ai11 * xi1;
               yi1 PEQ ai10 * xr0; ai10 = A0[3];
               yr2 MEQ ai21 * xi1;
               yi2 PEQ ai20 * xr0; ai20 = A0[5];
               yr3 MEQ ai31 * xi1;
               yi3 PEQ ai30 * xr0; xr0 = p0; ai30 = A0[7];

               yr0 += ar00 * xr0;
               yi0 += ar00 * xi0;
               yr1 += ar10 * xr0;
               yi1 += ar10 * xi0;
               yr2 += ar20 * xr0;
               yi2 += ar20 * xi0;
               yr3 += ar30 * xr0;
               yi3 += ar30 * xi0;

               yr0 MEQ ai00 * xi0;
               yi0 PEQ ai00 * xr0;
               yr1 MEQ ai10 * xi0;
               yi1 PEQ ai10 * xr0;
               yr2 MEQ ai20 * xi0;
               yi2 PEQ ai20 * xr0;
               yr3 MEQ ai30 * xi0;
               yi3 PEQ ai30 * xr0;
            }
            *Y   = yr0;
            Y[1] = yi0;
            A0 += incAm;
            Y[2] = yr1;
            Y[3] = yi1;
            A1 += incAm;
            Y[4] = yr2;
            Y[5] = yi2;
            Y[6] = yr3;
            Y[7] = yi3;
            Y += 8;
         }
         while(Y != stY);
      }
      if (M-M4) gemvMlt4(M-M4, N, A0, lda, X, beta, Y);
   }
   else if (N && M) gemvNlt2(M, N, A, lda, X, beta, Y);
}
void Mjoin(PATL,gemv0)
   (const int M, const int N, const SCALAR alpha,
    const TYPE *A, const int lda, const TYPE *X, const int incX,
    const SCALAR beta, TYPE *Y, const int incY)
{
   gemvN4x2(M, N, A, lda, X, beta, Y);
}

