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

#ifdef Conj_
   #define PEQ -=
   #define MEQ +=
#else
   #define PEQ +=
   #define MEQ -=
#endif


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
   int i, j;
   const int M2 = M<<1, N2 = N<<1;
   #ifdef BETAX
      const TYPE rbeta = *beta, ibeta = beta[1];
   #elif defined(BETAXI0)
      const TYPE rbeta = *beta;
   #endif
   register TYPE ra, ia, rx, ix, ry, iy;

   for (j=0; j != M2; j += 2)
   {
      #ifdef BETA0
         ry = iy = ATL_rzero;
      #elif defined(BETAX)
         rx = rbeta; ix = ibeta;
         ra = Y[j]; ia = Y[j+1];
         ry = ra * rx - ia * ix;
         iy = ra * ix + ia * rx;
      #else
         ry = Y[j];
         iy = Y[j+1];
         #ifdef BETAXI0
            ry *= rbeta;
            iy *= rbeta;
         #endif
      #endif
      for (i=0; i != N2; i += 2)
      {
         ra = A[i+j*lda]; ia = A[i+1+j*lda];
         rx = X[i]; ix = X[i+1];
         ry += ra * rx;
         iy += ra * ix;
         ry MEQ ia * ix;
         iy PEQ ia * rx;
      }
      Y[j] = ry;
      Y[j+1] = iy;
   }
}
