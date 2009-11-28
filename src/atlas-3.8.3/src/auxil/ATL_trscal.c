/*
 *             Automatically Tuned Linear Algebra Software v3.8.3
 *                   (C) Copyright 1999 Antoine P. Petitet
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

void Mjoin( PATL, trscal )
(
   const enum ATLAS_UPLO      UPLO,
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   TYPE                       * A,
   const int                  LDA
)
{
/*
 * Purpose
 * =======
 *
 * ATL_trscal scales a trapezoidal m-by-n matrix A by the scalar alpha.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   int                        i, incA, j, mn;
#ifdef TCPLX
   TYPE                       tmp;
#endif
   TYPE                       * a_j;
/* ..
 * .. Executable Statements ..
 *
 */
   if( UPLO == AtlasLower )
   {
      incA = ( ( LDA - M + 1 ) SHIFT ); mn = Mmin( M, N );

      if( SCALAR_IS_ZERO( ALPHA ) )
      {
         for( j = 0; j < mn; j++ )
         {
#ifdef TREAL
            for( i = j; i < M; i++, A += (1 SHIFT) ) { *A = ATL_rzero; }
#else
            for( i = j; i < M; i++, A += (1 SHIFT) ) { *A = A[1] = ATL_rzero; }
#endif
            A += incA;
            incA += (1 SHIFT);
         }
      }
      else if( !( SCALAR_IS_ONE( ALPHA ) ) )
      {
         for( j = 0; j < mn; j++ )
         {
            for( i = j; i < M; i++, A += (1 SHIFT) )
#ifdef TREAL
            { *A *= ALPHA; }
#else
            {
               tmp  = (*A);
               *A   = tmp * ALPHA[0] - A[1] * ALPHA[1];
               A[1] = tmp * ALPHA[1] + A[1] * ALPHA[0];
            }
#endif
            A += incA;
            incA += (1 SHIFT);
         }
      }
   }
   else
   {
      incA = ( LDA SHIFT );

      if( SCALAR_IS_ZERO( ALPHA ) )
      {
         for( j = 0, mn = M - N; j < N; j++, mn++, A += incA )
         {
            a_j = A;
            for( i = 0; i <= mn; i++, a_j += (1 SHIFT) )
#ifdef TREAL
            { *a_j = ATL_rzero; }
#else
            { *a_j = a_j[1] = ATL_rzero; }
#endif
         }
      }
      else if( !( SCALAR_IS_ONE( ALPHA ) ) )
      {
         for( j = 0, mn = M - N; j < N; j++, mn++, A += incA )
         {
            a_j = A;
            for( i = 0; i <= mn; i++, a_j += (1 SHIFT) )
#ifdef TREAL
            { *a_j *= ALPHA; }
#else
            {
               tmp = (*a_j);
               (*a_j) = tmp * ALPHA[0] - a_j[1] * ALPHA[1];
               a_j[1] = tmp * ALPHA[1] + a_j[1] * ALPHA[0];
            }
#endif
         }
      }
   }
}
