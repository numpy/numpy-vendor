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

static void gemv2x4(const int M, const int N, const TYPE *A, const int lda,
                    const TYPE *X, const SCALAR beta, TYPE *Y)
{
   const int N2 = (N>>1)<<1, M4 = (M>>2)<<2, nr = (N-N2);
   const int incAm = ((lda<<2) - N2)<<1;
   const TYPE *A0 = A, *A1 = A+(lda<<1), *A2 = A1+(lda<<1), *A3 = A2+(lda<<1);
   const TYPE *x, *stX = X + (N2<<1);
   TYPE *stY = Y + (M4<<1);
   #ifdef BETAX
      const TYPE rbeta = *beta, ibeta = beta[1];
   #elif defined(BETAXI0)
      const TYPE rbeta = *beta;
   #endif
   register TYPE yr0, yi0, yr1, yi1, yr2, yi2, yr3, yi3;
   register TYPE xr0, xi0, xr1, xi1;
   register TYPE ar00, ar01, ar02, ar03, ar10, ar11, ar12, ar13;
   register TYPE ai00, ai01, ai02, ai03, ai10, ai11, ai12, ai13;
   #ifdef ATL_AltiVec
      int cwrd = ATL_MulBySize(N)>>4;
      if (cwrd >= 64) cwrd = ATL_GetCtrl(512, (cwrd+31)>>5, 0);
      else cwrd = ATL_GetCtrl(64, (cwrd+3)>>2, 4);
   #endif


   if (N2)
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

            #ifdef BETA0
               yr0 = yi0 = yr1 = yi1 = yr2 = yi2 = yr3 = yi3 = ATL_rzero;
            #elif defined(BETAX)
               ar00 =   *Y; ai00 = Y[1]; ar01 = Y[2]; ai01 = Y[3];
               ar02 = Y[4]; ai02 = Y[5]; ar03 = Y[6]; ai03 = Y[7];
               yr0 = rbeta * ar00 - ibeta * ai00;
               yi0 = rbeta * ai00 + ibeta * ar00;
               yr1 = rbeta * ar01 - ibeta * ai01;
               yi1 = rbeta * ai01 + ibeta * ar01;
               yr2 = rbeta * ar02 - ibeta * ai02;
               yi2 = rbeta * ai02 + ibeta * ar02;
               yr3 = rbeta * ar03 - ibeta * ai03;
               yi3 = rbeta * ai03 + ibeta * ar03;
            #else
               yr0 =   *Y; yi0 = Y[1]; yr1 = Y[2]; yi1 = Y[3];
               yr2 = Y[4]; yi2 = Y[5]; yr3 = Y[6]; yi3 = Y[7];
               #ifdef BETAXI0
                  yr0 *= rbeta; yi0 *= rbeta;
                  yr1 *= rbeta; yi1 *= rbeta;
                  yr2 *= rbeta; yi2 *= rbeta;
                  yr3 *= rbeta; yi3 *= rbeta;
               #endif
            #endif
            ar00 = *A0;   ar01 = *A1;    ar02 = *A2;   ar03 = *A3;
            ai00 = A0[1]; ai01 = A1[1];  ai02 = A2[1]; ai03 = A3[1];
            ar10 = A0[2]; ar11 = A1[2];  ar12 = A2[2]; ar13 = A3[2];
            ai10 = A0[3]; ai11 = A1[3];  ai12 = A2[3]; ai13 = A3[3];
            A0 += 4; A1 += 4; A2 += 4; A3 += 4;
            xr0 = *X; xi0 = X[1]; xr1 = X[2]; xi1 = X[3]; x = X + 4;
            if (N2 > 3)
            {
               do
               {
                  yr0 += xr0 * ar00;
                  yr1 += xr0 * ar01;
                  yr2 += xr0 * ar02;
                  yr3 += xr0 * ar03;
                  yi0 += xi0 * ar00; ar00 = *A0;
                  yi1 += xi0 * ar01; ar01 = *A1;
                  yi2 += xi0 * ar02; ar02 = *A2;
                  yi3 += xi0 * ar03; ar03 = *A3;
                  #ifndef Conj_
                     yr0 -= xi0 * ai00;
                     yr1 -= xi0 * ai01;
                     yr2 -= xi0 * ai02;
                     yr3 -= xi0 * ai03; xi0 = x[1];
                     yi0 += xr0 * ai00; ai00 = A0[1];
                     yi1 += xr0 * ai01; ai01 = A1[1];
                     yi2 += xr0 * ai02; ai02 = A2[1];
                     yi3 += xr0 * ai03; ai03 = A3[1]; xr0 = *x;
                  #else
                     yr0 += xi0 * ai00;
                     yr1 += xi0 * ai01;
                     yr2 += xi0 * ai02;
                     yr3 += xi0 * ai03; xi0 = x[1];
                     yi0 -= xr0 * ai00; ai00 = A0[1];
                     yi1 -= xr0 * ai01; ai01 = A1[1];
                     yi2 -= xr0 * ai02; ai02 = A2[1];
                     yi3 -= xr0 * ai03; ai03 = A3[1]; xr0 = *x;
                  #endif

                  yr0 += xr1 * ar10;
                  yr1 += xr1 * ar11;
                  yr2 += xr1 * ar12;
                  yr3 += xr1 * ar13;
                  yi0 += xi1 * ar10; ar10 = A0[2];
                  yi1 += xi1 * ar11; ar11 = A1[2];
                  yi2 += xi1 * ar12; ar12 = A2[2];
                  yi3 += xi1 * ar13; ar13 = A3[2];
                  #ifndef Conj_
                     yr0 -= xi1 * ai10;
                     yr1 -= xi1 * ai11;
                     yr2 -= xi1 * ai12;
                     yr3 -= xi1 * ai13; xi1 = x[3];
                     yi0 += xr1 * ai10; ai10 = A0[3]; A0 += 4;
                     yi1 += xr1 * ai11; ai11 = A1[3]; A1 += 4;
                     yi2 += xr1 * ai12; ai12 = A2[3]; A2 += 4;
                     yi3 += xr1 * ai13; ai13 = A3[3]; A3 += 4; xr1=x[2]; x += 4;
                  #else
                     yr0 += xi1 * ai10;
                     yr1 += xi1 * ai11;
                     yr2 += xi1 * ai12;
                     yr3 += xi1 * ai13; xi1 = x[3];
                     yi0 -= xr1 * ai10; ai10 = A0[3]; A0 += 4;
                     yi1 -= xr1 * ai11; ai11 = A1[3]; A1 += 4;
                     yi2 -= xr1 * ai12; ai12 = A2[3]; A2 += 4;
                     yi3 -= xr1 * ai13; ai13 = A3[3]; A3 += 4; xr1=x[2]; x += 4;
                  #endif
               }
               while (x != stX);
            }
            if (!nr) /* fetching done, just do 2 remaining ops */
            {
               yr0 += xr0 * ar00;
               yr1 += xr0 * ar01;
               yr2 += xr0 * ar02;
               yr3 += xr0 * ar03;
               yi0 += xi0 * ar00;
               yi1 += xi0 * ar01;
               yi2 += xi0 * ar02;
               yi3 += xi0 * ar03;
               #ifndef Conj_
                  yr0 -= xi0 * ai00;
                  yr1 -= xi0 * ai01;
                  yr2 -= xi0 * ai02;
                  yr3 -= xi0 * ai03;
                  yi0 += xr0 * ai00;
                  yi1 += xr0 * ai01;
                  yi2 += xr0 * ai02;
                  yi3 += xr0 * ai03;
               #else
                  yr0 += xi0 * ai00;
                  yr1 += xi0 * ai01;
                  yr2 += xi0 * ai02;
                  yr3 += xi0 * ai03;
                  yi0 -= xr0 * ai00;
                  yi1 -= xr0 * ai01;
                  yi2 -= xr0 * ai02;
                  yi3 -= xr0 * ai03;
               #endif

               yr0 += xr1 * ar10;
               yr1 += xr1 * ar11;
               yr2 += xr1 * ar12;
               yr3 += xr1 * ar13;
               yi0 += xi1 * ar10;
               yi1 += xi1 * ar11;
               yi2 += xi1 * ar12;
               yi3 += xi1 * ar13;
               #ifndef Conj_
                  yr0 -= xi1 * ai10;
                  yr1 -= xi1 * ai11;
                  yr2 -= xi1 * ai12;
                  yr3 -= xi1 * ai13;
                  yi0 += xr1 * ai10;
                  yi1 += xr1 * ai11;
                  yi2 += xr1 * ai12;
                  yi3 += xr1 * ai13;
               #else
                  yr0 += xi1 * ai10;
                  yr1 += xi1 * ai11;
                  yr2 += xi1 * ai12;
                  yr3 += xi1 * ai13;
                  yi0 -= xr1 * ai10;
                  yi1 -= xr1 * ai11;
                  yi2 -= xr1 * ai12;
                  yi3 -= xr1 * ai13;
               #endif
            }
            else  /* 3 rows left to do, 2 rows prefetched */
            {
               yr0 += xr0 * ar00;
               yr1 += xr0 * ar01;
               yr2 += xr0 * ar02;
               yr3 += xr0 * ar03;
               yi0 += xi0 * ar00; ar00 = *A0;
               yi1 += xi0 * ar01; ar01 = *A1;
               yi2 += xi0 * ar02; ar02 = *A2;
               yi3 += xi0 * ar03; ar03 = *A3;
               #ifndef Conj_
                  yr0 -= xi0 * ai00;
                  yr1 -= xi0 * ai01;
                  yr2 -= xi0 * ai02;
                  yr3 -= xi0 * ai03; xi0 = x[1];
                  yi0 += xr0 * ai00; ai00 = A0[1];
                  yi1 += xr0 * ai01; ai01 = A1[1];
                  yi2 += xr0 * ai02; ai02 = A2[1];
                  yi3 += xr0 * ai03; ai03 = A3[1]; xr0 = *x;
               #else
                  yr0 += xi0 * ai00;
                  yr1 += xi0 * ai01;
                  yr2 += xi0 * ai02;
                  yr3 += xi0 * ai03; xi0 = x[1];
                  yi0 -= xr0 * ai00; ai00 = A0[1];
                  yi1 -= xr0 * ai01; ai01 = A1[1];
                  yi2 -= xr0 * ai02; ai02 = A2[1];
                  yi3 -= xr0 * ai03; ai03 = A3[1]; xr0 = *x;
               #endif

               yr0 += xr1 * ar10;
               yr1 += xr1 * ar11;
               yr2 += xr1 * ar12;
               yr3 += xr1 * ar13;
               yi0 += xi1 * ar10;
               yi1 += xi1 * ar11;
               yi2 += xi1 * ar12;
               yi3 += xi1 * ar13;
               #ifndef Conj_
                  yr0 -= xi1 * ai10;
                  yr1 -= xi1 * ai11;
                  yr2 -= xi1 * ai12;
                  yr3 -= xi1 * ai13;
                  yi0 += xr1 * ai10;
                  yi1 += xr1 * ai11;
                  yi2 += xr1 * ai12;
                  yi3 += xr1 * ai13;
               #else
                  yr0 += xi1 * ai10;
                  yr1 += xi1 * ai11;
                  yr2 += xi1 * ai12;
                  yr3 += xi1 * ai13;
                  yi0 -= xr1 * ai10;
                  yi1 -= xr1 * ai11;
                  yi2 -= xr1 * ai12;
                  yi3 -= xr1 * ai13;
               #endif

               yr0 += xr0 * ar00;
               yr1 += xr0 * ar01;
               yr2 += xr0 * ar02;
               yr3 += xr0 * ar03;
               yi0 += xi0 * ar00;
               yi1 += xi0 * ar01;
               yi2 += xi0 * ar02;
               yi3 += xi0 * ar03;
               #ifndef Conj_
                  yr0 -= xi0 * ai00;
                  yr1 -= xi0 * ai01;
                  yr2 -= xi0 * ai02;
                  yr3 -= xi0 * ai03;
                  yi0 += xr0 * ai00;
                  yi1 += xr0 * ai01;
                  yi2 += xr0 * ai02;
                  yi3 += xr0 * ai03;
               #else
                  yr0 += xi0 * ai00;
                  yr1 += xi0 * ai01;
                  yr2 += xi0 * ai02;
                  yr3 += xi0 * ai03;
                  yi0 -= xr0 * ai00;
                  yi1 -= xr0 * ai01;
                  yi2 -= xr0 * ai02;
                  yi3 -= xr0 * ai03;
               #endif
            }
            *Y   = yr0;
            A0 += incAm;
            Y[1] = yi0;
            A1 += incAm;
            Y[2] = yr1;
            A2 += incAm;
            Y[3] = yi1;
            A3 += incAm;
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
   else if (M && N) gemvNlt2(M, N, A, lda, X, beta, Y);
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
   gemv2x4(M, N, A, lda, X, beta, Y);
}
