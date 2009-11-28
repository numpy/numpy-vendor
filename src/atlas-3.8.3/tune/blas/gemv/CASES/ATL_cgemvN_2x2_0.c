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

void Mjoin(PATL,gemv0)
   (const int M, const int N, const SCALAR alpha,
    const TYPE *A, const int lda, const TYPE *X, const int incX,
    const SCALAR beta, TYPE *Y, const int incY)
{
   const int N2 = (N>>1)<<1, M2 = (M>>1)<<1, nr = N-N2;
   const int incA = lda<<2, incAm = 4 - N2 * (lda<<1);
   const TYPE *x, *stX = X + (N2<<1), *A0 = A, *A1 = A + (lda<<1);
   TYPE *stY = Y + (M2<<1);
   register TYPE xr0, xi0, xr1, xi1, yr0, yi0, yr1, yi1, yyr0, yyi0, yyr1, yyi1;
   register TYPE ar00, ai00, ar10, ai10, ar01, ai01, ar11, ai11;
   register TYPE m0, m1, m2, m3;

   if (N2)
   {
      if (M2)
      {
         do
         {
            #ifdef BETA0
               yr0 = yi0 = yr1 = yi1 = yyr0 = yyi0 = yyr1 = yyi1 = ATL_rzero;
            #elif defined(BETAX)
               xr0 = *beta; xi0 = beta[1];
               ar00 = *Y; ai00 = Y[1]; ar10 = Y[2]; ai10 = Y[3];
               yyr0 = ar00 * xr0 - ai00 * xi0;
               yyi0 = ar00 * xi0 + ai00 * xr0;
               yyr1 = ar10 * xr0 - ai10 * xi0;
               yyi1 = ar10 * xi0 + ai10 * xr0;
               yr0 = yi0 = yr1 = yi1 = ATL_rzero;
            #else
               yyr0 = *Y; yyi0 = Y[1]; yyr1 = Y[2]; yyi1 = Y[3];
               yr0 = yi0 = yr1 = yi1 = ATL_rzero;
               #ifdef BETAXI0
                 xr0 = *beta;
                 yyr0 *= xr0; yyi0 *= xr0;
                 yyr1 *= xr0; yyi1 *= xr0;
               #endif
            #endif
            x = X + 4;
            ar00 = *A0; ai00 = A0[1]; ar10 = A0[2]; ai10 = A0[3]; A0 += incA;
            xr0 = *X; xi0 = X[1]; xr1 = X[2]; xi1 = X[3];
            ar01 = *A1; ai01 = A1[1]; ar11 = A1[2]; ai11 = A1[3]; A1 += incA;
            m0 = ar00 * xr0;
            m1 = ar01 * xr1;
            m2 = ar00 * xi0;
            m3 = ar01 * xi1;
            if (N2 != 2)
            {
               do
               {
                  ar00 = *A0; ar01 = *A1;
                  yr0  += m0; m0 = ai00 * xi0;
                  yyr0 += m1; m1 = ai01 * xi1;
                  yi0  += m2; m2 = ai00 * xr0; ai00 = A0[1];
                  yyi0 += m3; m3 = ai01 * xr1; ai01 = A1[1];
                  #ifndef Conj_
                     yr0  -= m0; m0 = ar10 * xr0;
                     yyr0 -= m1; m1 = ar11 * xr1;
                     yi0  += m2; m2 = ar10 * xi0; ar10 = A0[2];
                     yyi0 += m3; m3 = ar11 * xi1; ar11 = A1[2];
                  #else
                     yr0  += m0; m0 = ar10 * xr0;
                     yyr0 += m1; m1 = ar11 * xr1;
                     yi0  -= m2; m2 = ar10 * xi0; ar10 = A0[2];
                     yyi0 -= m3; m3 = ar11 * xi1; ar11 = A1[2];
                  #endif

                  yr1  += m0; m0 = ai10 * xi0; xi0 = x[1];
                  yyr1 += m1; m1 = ai11 * xi1; xi1 = x[3];
                  yi1  += m2; m2 = ai10 * xr0; xr0 = *x;
                              ai10 = A0[3]; A0 += incA;
                  yyi1 += m3; m3 = ai11 * xr1; xr1 = x[2]; x += 4;
                              ai11 = A1[3]; A1 += incA;
                  #ifndef Conj_
                     yr1  -= m0; m0 = ar00 * xr0;
                     yyr1 -= m1; m1 = ar01 * xr1;
                     yi1  += m2; m2 = ar00 * xi0;
                     yyi1 += m3; m3 = ar01 * xi1;
                  #else
                     yr1  += m0; m0 = ar00 * xr0;
                     yyr1 += m1; m1 = ar01 * xr1;
                     yi1  -= m2; m2 = ar00 * xi0;
                     yyi1 -= m3; m3 = ar01 * xi1;
                  #endif
               }
               while (x != stX);
            }
            if (!nr)  /* drain pipe, no cleanup needed */
            {
               yr0  += m0; m0 = ai00 * xi0;
               yyr0 += m1; m1 = ai01 * xi1;
               yi0  += m2; m2 = ai00 * xr0;
               yyi0 += m3; m3 = ai01 * xr1;
               #ifndef Conj_
                  yr0  -= m0; m0 = ar10 * xr0;
                  yyr0 -= m1; m1 = ar11 * xr1;
                  yi0  += m2; m2 = ar10 * xi0;
                  yyi0 += m3; m3 = ar11 * xi1;
               #else
                  yr0  += m0; m0 = ar10 * xr0;
                  yyr0 += m1; m1 = ar11 * xr1;
                  yi0  -= m2; m2 = ar10 * xi0;
                  yyi0 -= m3; m3 = ar11 * xi1;
               #endif

               yr1  += m0; m0 = ai10 * xi0;
               yyr1 += m1; m1 = ai11 * xi1;
               yi1  += m2; m2 = ai10 * xr0;
               yyi1 += m3; m3 = ai11 * xr1;
               #ifndef Conj_
                  yr1  -= m0;
                  yyr1 -= m1;
                  yi1  += m2;
                  yyi1 += m3;
               #else
                  yr1  += m0;
                  yyr1 += m1;
                  yi1  -= m2;
                  yyi1 -= m3;
               #endif
            }
            else  /* drain pipe, and do 1 extra elt of X */
            {
               ar00 = *A0;
               yr0  += m0; m0 = ai00 * xi0;
               yyr0 += m1; m1 = ai01 * xi1;
               yi0  += m2; m2 = ai00 * xr0; ai00 = A0[1];
               yyi0 += m3; m3 = ai01 * xr1;
               #ifndef Conj_
                  yr0  -= m0; m0 = ar10 * xr0;
                  yyr0 -= m1; m1 = ar11 * xr1;
                  yi0  += m2; m2 = ar10 * xi0; ar10 = A0[2];
                  yyi0 += m3; m3 = ar11 * xi1;
               #else
                  yr0  += m0; m0 = ar10 * xr0;
                  yyr0 += m1; m1 = ar11 * xr1;
                  yi0  -= m2; m2 = ar10 * xi0; ar10 = A0[2];
                  yyi0 -= m3; m3 = ar11 * xi1;
               #endif

               yr1  += m0; m0 = ai10 * xi0; xi0 = x[1];
               yyr1 += m1; m1 = ai11 * xi1;
               yi1  += m2; m2 = ai10 * xr0; xr0 = *x; ai10 = A0[3];
               yyi1 += m3; m3 = ai11 * xr1;
               #ifndef Conj_
                  yr1  -= m0; m0 = ar00 * xr0;
                  yyr1 -= m1; m1 = ar00 * xi0;
                  yi1  += m2; m2 = ar10 * xr0;
                  yyi1 += m3; m3 = ar10 * xi0;
               #else
                  yr1  += m0; m0 = ar00 * xr0;
                  yyr1 += m1; m1 = ar00 * xi0;
                  yi1  -= m2; m2 = ar10 * xr0;
                  yyi1 -= m3; m3 = ar10 * xi0;
               #endif

               yr0 += m0; m0 = ai00 * xi0;
               yi0 += m1; m1 = ai00 * xr0;
               yr1 += m2; m2 = ai10 * xi0;
               yi1 += m3; m3 = ai10 * xr0;
               #ifndef Conj_
                  yr0 -= m0;
                  yi0 += m1;
                  yr1 -= m2;
                  yi1 += m3;
               #else
                  yr0 += m0;
                  yi0 -= m1;
                  yr1 += m2;
                  yi1 -= m3;
               #endif
            }
            A0 += incAm;
            A1 += incAm;
            yr0 += yyr0;
            yi0 += yyi0;
            yr1 += yyr1;
            yi1 += yyi1;
            *Y = yr0;
            Y[1] = yi0;
            Y[2] = yr1;
            Y[3] = yi1;
            Y += 4;
         }
         while(Y != stY);
      }
      if (M-M2) gemvMlt4(M-M2, N, A0, lda, X, beta, Y);
   }
   else gemvNlt2(M, N, A, lda, X, beta, Y);
}

