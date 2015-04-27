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

#include "atlas_lvl3.h"

#ifdef Conj_
void Mjoin(PATL,Mjoin(Mjoin(Mjoin(Mjoin(gemvNc,NM),_x1),BNM),_y1))
   (const int M, const int N, const SCALAR alpha,
    const TYPE *A, const int lda, const TYPE *X, const int incX,
    const SCALAR beta, TYPE *Y, const int incY)
/*
 *  Use trick: y <- alpha * conj(A) + y ==>
 *             trans(y) <- alpha * trans(x) * ConjTrans(A) + trans(y)
 */
{
   TYPE one[2] = {ATL_rone, ATL_rzero};
   Mjoin(PATL,NCmmJIK)(AtlasTrans, AtlasConjTrans, 1, M, N, one, X, 1,
                       A, lda, beta, Y, 1);
}
#else
void Mjoin(PATL,Mjoin(Mjoin(Mjoin(Mjoin(gemvN,NM),_x1),BNM),_y1))
   (const int M, const int N, const SCALAR alpha,
    const TYPE *A, const int lda, const TYPE *X, const int incX,
    const SCALAR beta, TYPE *Y, const int incY)
{
   TYPE one[2] = {ATL_rone, ATL_rzero};
   Mjoin(PATL,NCmmJIK)(AtlasNoTrans, AtlasNoTrans, M, 1, N, one, A, lda,
                       X, N, beta, Y, M);
}
#endif
