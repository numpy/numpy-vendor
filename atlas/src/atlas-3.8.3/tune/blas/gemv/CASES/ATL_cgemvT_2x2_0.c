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

#ifdef Conj_
   #define PEQ -=
   #define MEQ +=
#else
   #define PEQ +=
   #define MEQ -=
#endif

static void gemvMlt4(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
{
   int j;
   const int lda2 = lda<<1;
   #ifndef BETA0
      register TYPE yr0, yi0, yr1, yi1;
   #endif
   #if defined(BETAXI0)
      const register TYPE rbeta = *beta;
   #elif defined(BETAX)
      const register TYPE xr0 = *beta, xi0 = beta[1];
   #endif

   for (j=M; j; j--, A += lda2, Y += 2)
   {
   #ifdef BETA0
      #ifndef Conj_
         Mjoin(PATL,dotu_sub)(N, A, 1, X, 1, Y);
      #else
         Mjoin(PATL,dotc_sub)(N, A, 1, X, 1, Y);
      #endif
   #else
      #ifdef BETAX
         yr0 = *Y;
         yi0 = Y[1];
         yr1 = yr0 * xr0 - yi0 * xi0;
         yi1 = yr0 * xi0 + yi0 * xr0;
      #else
         yr1 = *Y;
         yi1 = Y[1];
         #ifdef BETAXI0
            yr1 *= rbeta;
            yi1 *= rbeta;
         #endif
      #endif
      #ifndef Conj_
         Mjoin(PATL,dotu_sub)(N, A, 1, X, 1, Y);
      #else
         Mjoin(PATL,dotc_sub)(N, A, 1, X, 1, Y);
      #endif
      *Y += yr1;
      Y[1] += yi1;
   #endif
   }
}

static void gemvNlt2(const int M, const int N, const TYPE *A, const int lda,
                     const TYPE *X, const SCALAR beta, TYPE *Y)
{
   #ifdef BETA0
      #ifndef Conj_
         Mjoin(PATL,cpsc)(M, X, A, lda, Y, 1);
      #else
         Mjoin(PATL,moveConj)(M, X, A, lda, Y, 1);
      #endif
   #elif defined(BETAX) || defined(BETAXI0)
      #ifndef Conj_
         Mjoin(PATL,axpby)(M, X, A, lda, beta, Y, 1);
      #else
         Mjoin(PATL,axpbyConj)(M, X, A, lda, beta, Y, 1);
      #endif
   #else
      #ifndef Conj_
         Mjoin(PATL,axpy)(M, X, A, lda, Y, 1);
      #else
         Mjoin(PATL,axpyConj)(M, X, A, lda, Y, 1);
      #endif
   #endif
}


#ifdef Conj_
   #define gemv0 Mjoin(Mjoin(Mjoin(Mjoin(gemvC,NM),_x1),BNM),_y1)
#else
   #define gemv0 Mjoin(Mjoin(Mjoin(Mjoin(gemvT,NM),_x1),BNM),_y1)
#endif
void Mjoin(PATL,gemv0)
   (const int M, const int N, const SCALAR alpha,
    const TYPE *A, const int lda, const TYPE *X, const int incX,
    const SCALAR beta, TYPE *Y, const int incY)
{
   const int N2 = (N>>1)<<1, M2 = (M>>1)<<1;
   const int nr = N-N2, incAm = ((lda<<1) - N2)<<1;
   const TYPE *x, *stX = X + (N2<<1), *A0 = A, *A1 = A + (lda<<1);
   TYPE *stY = Y + (M2<<1);
   register TYPE xr0, xi0, xr1, xi1, yr0, yi0, yr1, yi1;
   register TYPE ar00, ai00, ar10, ai10, ar01, ai01, ar11, ai11;
   register TYPE m0, m1, m2, m3;
   #ifndef ATL_PFD
      #define ATL_PFD 8
   #endif

   if (N2)
   {
      if (M2)
      {
         do
         {
            #ifdef BETA0
               yr0 = yi0 = yr1 = yi1 = ATL_rzero;
            #elif defined(BETAX)
               xr0 = *beta; xi0 = beta[1];
               ar00 = *Y; ai00 = Y[1]; ar10 = Y[2]; ai10 = Y[3];
               yr0 = ar00 * xr0 - ai00 * xi0;
               yi0 = ar00 * xi0 + ai00 * xr0;
               yr1 = ar10 * xr0 - ai10 * xi0;
               yi1 = ar10 * xi0 + ai10 * xr0;
            #else
               yr0 = *Y; yi0 = Y[1]; yr1 = Y[2]; yi1 = Y[3];
               #ifdef BETAXI0
                  xr0 = *beta;
                  yr0 *= xr0; yi0 *= xr0;
                  yr1 *= xr0; yi1 *= xr0;
               #endif
            #endif
            x = X + 4;
            ar00 = *A0; ai00 = A0[1]; ar10 = A0[2]; ai10 = A0[3]; A0 += 4;
            xr0 = *X; xi0 = X[1]; xr1 = X[2]; xi1 = X[3];
            ar01 = *A1; ai01 = A1[1]; ar11 = A1[2]; ai11 = A1[3]; A1 += 4;
            m0 = ar00 * xr0;
            m1 = ar01 * xr0;
            m2 = ar00 * xi0;
            m3 = ar01 * xi0;
            if (N2 != 2)
            {
               do
               {
                  ar00 = *A0; ar01 = *A1;
                  yr0 += m0; m0 = ai00 * xi0;
                                  ATL_pfl1R(A0+ATL_PFD); ATL_pfl1R(A1+ATL_PFD);
                  yr1 += m1; m1 = ai01 * xi0; xi0 = x[1];
                  yi0 += m2; m2 = ai00 * xr0; ai00 = A0[1];
                  yi1 += m3; m3 = ai01 * xr0; ai01 = A1[1]; xr0 = *x;
                  #ifndef Conj_
                     yr0 -= m0; m0 = ar10 * xr1;
                     yr1 -= m1; m1 = ar11 * xr1;
                     yi0 += m2; m2 = ar10 * xi1; ar10 = A0[2];
                     yi1 += m3; m3 = ar11 * xi1; ar11 = A1[2];
                  #else
                     yr0 += m0; m0 = ar10 * xr1;
                     yr1 += m1; m1 = ar11 * xr1;
                     yi0 -= m2; m2 = ar10 * xi1; ar10 = A0[2];
                     yi1 -= m3; m3 = ar11 * xi1; ar11 = A1[2];
                  #endif

                  yr0 += m0; m0 = ai10 * xi1;
                  yr1 += m1; m1 = ai11 * xi1; xi1 = x[3];
                  yi0 += m2; m2 = ai10 * xr1; ai10 = A0[3]; A0 += 4;
                  yi1 += m3; m3 = ai11 * xr1; ai11 = A1[3]; A1 += 4;
                         xr1 = x[2]; x += 4;
                  #ifndef Conj_
                     yr0 -= m0; m0 = ar00 * xr0;
                     yr1 -= m1; m1 = ar01 * xr0;
                     yi0 += m2; m2 = ar00 * xi0;
                     yi1 += m3; m3 = ar01 * xi0;
                  #else
                     yr0 += m0; m0 = ar00 * xr0;
                     yr1 += m1; m1 = ar01 * xr0;
                     yi0 -= m2; m2 = ar00 * xi0;
                     yi1 -= m3; m3 = ar01 * xi0;
                  #endif
               }
               while (x != stX);
            }
/*
 *          Drain pipes & store to Y
 */
            if (nr)
            {
               ar00 = *A0; ar01 = *A1;
               yr0 += m0; m0 = ai00 * xi0;
               yr1 += m1; m1 = ai01 * xi0; xi0 = x[1];
               yi0 += m2; m2 = ai00 * xr0; ai00 = A0[1];
               yi1 += m3; m3 = ai01 * xr0; ai01 = A1[1]; xr0 = *x;
               #ifndef Conj_
                  yr0 -= m0; m0 = ar10 * xr1;
                  yr1 -= m1; m1 = ar11 * xr1;
                  yi0 += m2; m2 = ar10 * xi1;
                  yi1 += m3; m3 = ar11 * xi1;
               #else
                  yr0 += m0; m0 = ar10 * xr1;
                  yr1 += m1; m1 = ar11 * xr1;
                  yi0 -= m2; m2 = ar10 * xi1;
                  yi1 -= m3; m3 = ar11 * xi1;
               #endif

               yr0 += m0; m0 = ai10 * xi1;
               yr1 += m1; m1 = ai11 * xi1;
               yi0 += m2; m2 = ai10 * xr1;
               yi1 += m3; m3 = ai11 * xr1;
               #ifndef Conj_
                  yr0 -= m0; m0 = ar00 * xr0;
                  yr1 -= m1; m1 = ar01 * xr0;
                  yi0 += m2; m2 = ar00 * xi0;
                  yi1 += m3; m3 = ar01 * xi0;
               #else
                  yr0 += m0; m0 = ar00 * xr0;
                  yr1 += m1; m1 = ar01 * xr0;
                  yi0 -= m2; m2 = ar00 * xi0;
                  yi1 -= m3; m3 = ar01 * xi0;
               #endif

               yr0 += m0; m0 = ai00 * xi0;
               yr1 += m1; m1 = ai01 * xi0;
               yi0 += m2; m2 = ai00 * xr0;
               yi1 += m3; m3 = ai01 * xr0;
               #ifndef Conj_
                  yr0 -= m0;
                  yr1 -= m1;
                  yi0 += m2;
                  yi1 += m3;
               #else
                  yr0 += m0;
                  yr1 += m1;
                  yi0 -= m2;
                  yi1 -= m3;
               #endif
            }
            else
            {
               yr0 += m0; m0 = ai00 * xi0;
               yr1 += m1; m1 = ai01 * xi0;
               yi0 += m2; m2 = ai00 * xr0;
               yi1 += m3; m3 = ai01 * xr0;
               #ifndef Conj_
                  yr0 -= m0; m0 = ar10 * xr1;
                  yr1 -= m1; m1 = ar11 * xr1;
                  yi0 += m2; m2 = ar10 * xi1;
                  yi1 += m3; m3 = ar11 * xi1;
               #else
                  yr0 += m0; m0 = ar10 * xr1;
                  yr1 += m1; m1 = ar11 * xr1;
                  yi0 -= m2; m2 = ar10 * xi1;
                  yi1 -= m3; m3 = ar11 * xi1;
               #endif

               yr0 += m0; m0 = ai10 * xi1;
               yr1 += m1; m1 = ai11 * xi1;
               yi0 += m2; m2 = ai10 * xr1;
               yi1 += m3; m3 = ai11 * xr1;
               #ifndef Conj_
                  yr0 -= m0;
                  yr1 -= m1;
                  yi0 += m2;
                  yi1 += m3;
               #else
                  yr0 += m0;
                  yr1 += m1;
                  yi0 -= m2;
                  yi1 -= m3;
               #endif
            }

            *Y = yr0;
            A0 += incAm;
            Y[1] = yi0;
            Y[2] = yr1;
            A1 += incAm;
            Y[3] = yi1;
            Y += 4;
         }
         while (Y != stY);
      }
      if (M-M2) gemvMlt4(M-M2, N, A0, lda, X, beta, Y);
   }
   else if (M && N) gemvNlt2(M, N, A, lda, X, beta, Y);
}
