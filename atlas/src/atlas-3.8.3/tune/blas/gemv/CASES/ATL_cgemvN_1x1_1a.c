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
   int j;
   const int lda2 = lda<<1, N2 = N<<1;

   #ifdef Conj_
      #define axpycj Mjoin(PATL,axpyConj)
      Mjoin(PATL,axpbyConj)(M, X, A, 1, beta, Y, 1);
   #else
      #define axpycj Mjoin(PATL,axpy)
      Mjoin(PATL,axpby)(M, X, A, 1, beta, Y, 1);
   #endif
   if (N != 1)
   {
      A += lda2;
      X += 2;
      for (j=N-1; j; j--, A += lda2, X += 2) axpycj(M, X, A, 1, Y, 1);
   }
}
