/* ---------------------------------------------------------------------
 *
 * -- Automatically Tuned Linear Algebra Software (ATLAS)
 *    (C) Copyright 2000 All Rights Reserved
 *
 * -- ATLAS routine -- Version 3.2 -- December 25, 2000
 *
 * Author         : Antoine P. Petitet
 * Contributor(s) : R. Clint Whaley
 * Originally developed at the University of Tennessee,
 * Innovative Computing Laboratory, Knoxville TN, 37996-1301, USA.
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the ATLAS group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ---------------------------------------------------------------------
 */
/*
 * Include files
 */
#include "atlas_misc.h"
#include "atlas_mv.h"
#include "atlas_r1.h"
#include "atlas_level1.h"
#include "atlas_kernel2.h"
#include "atlas_reflvl2.h"
#include "atlas_lvl2.h"

void Mjoin( PATL, trmvUT )
(
   const enum ATLAS_DIAG      DIAG,
   const int                  N,         /* N > 0 assumed */
   const TYPE                 * A,
   const int                  LDA,
   TYPE                       * X
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, trmvUT ) performs the following matrix-vector operation
 *
 *    x := A'* x,
 *
 * where x is an n-element vector and  A is an n by n unit or  non-unit,
 * upper triangular matrix.
 *
 * This is a blocked version of the algorithm.  For a more detailed des-
 * cription of  the arguments of this function, see the reference imple-
 * mentation in the  ATLAS/src/blas/reference directory.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
#ifdef TREAL
#define    one                ATL_rone
#define    lda2               LDA
#else
   TYPE                       one[2] = { ATL_rone, ATL_rzero };
   const int                  lda2   = ( LDA SHIFT );
#endif
   TYPE                       * A0, * x0;
   int                        incA, incX, m1, mb, mb1, n, nb;
   void  (*trmv0)( const int, const TYPE *, const int, TYPE * );
#define    gemv0              Mjoin( PATL, gemvT_a1_x1_b1_y1  )
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_GetPartMVT( A, LDA, &mb, &nb );

   if( DIAG == AtlasNonUnit ) trmv0 = Mjoin( PATL, trmvUTN );
   else                       trmv0 = Mjoin( PATL, trmvUTU );

   mb1  = N - ( m1 = ( ( N - 1 ) / mb ) * mb );
   incA = ( incX  = (mb SHIFT) ) + mb * lda2;

   A += (m1 SHIFT) + m1 * lda2; X += (m1 SHIFT);
   A0 = (TYPE *)(A) - incX;     x0 = X;

   trmv0( mb1, A, LDA, X ); X -= incX; A -= incA;
   for( n = mb1; n < N; n += mb, A0 -= incA, A -= incA, X -= incX, x0 -= incX )
   {
      gemv0( n, mb, one, A0, LDA, X, 1, one, x0, 1 );
      trmv0( mb, A, LDA, X );
   }
/*
 * End of Mjoin( PATL, trmvUT )
 */
}
