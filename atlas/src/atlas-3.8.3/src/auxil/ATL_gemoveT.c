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

#ifdef Conj_
void Mjoin(PATL,gemoveC)
#else
void Mjoin(PATL,gemoveT)
#endif
   (const int N, const int M, const SCALAR alpha, const TYPE *A, const int lda,
    TYPE *C, const int ldc)
/*
 * C <- alpha * transpose(A), C is MxN, A is NxM
 * NOTE: C is written row-wise, on assumption you are copying to a smaller
 *       matrix.  Should be optimized, say by unrolling so you write 4 elts
 *       of C at a time.
 */
#ifdef TCPLX
{
   const TYPE ralpha = *alpha, calpha=alpha[1];
   const int incA = lda+lda;
   int i;

   for (i=0; i < M; i++, A += incA, C++)
   #ifdef Conj_
      Mjoin(PATL,moveConj)(N, alpha, A, 1, C, ldc);
   #else
      Mjoin(PATL,cpsc)(N, alpha, A, 1, C, ldc);
   #endif
}
#else  /* real code */
{
   int i, j;
   const int incA = lda - M, incC = 1 - ldc*M;

   if (alpha == ATL_rone)
   {
      for (j=N; j; j--, A += incA, C += incC)
         for (i=M; i; i--, C += ldc)
            *C = *A++;
   }
   else if (alpha == ATL_rnone)
   {
      for (j=N; j; j--, A += incA, C += incC)
         for (i=M; i; i--, C += ldc)
            *C = -(*A++);
   }
   else if (alpha == ATL_rzero)
      Mjoin(PATL,gezero)(M, N, C, ldc);
   else
   {
      for (j=N; j; j--, A += incA, C += incC)
         for (i=M; i; i--, C += ldc)
            *C = alpha*(*A++);
   }
}
#endif
