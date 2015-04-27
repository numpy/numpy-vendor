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

void Mjoin( PATL, tpmvUT )
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
 * Mjoin( PATL, tpmvUT ) performs the following matrix-vector operation
 *
 *    x := A'* x,
 *
 * where x is an n-element vector and  A is an n by n unit or  non-unit,
 * upper triangular matrix supplied in packed form.
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
   void                       (*tpmv0)( const int, const TYPE *, const int,
                              TYPE * );
#ifdef TREAL
#define    one                ATL_rone
#else
   TYPE                       one[2] = { ATL_rone, ATL_rzero };
#endif
   TYPE                       * A0, * x0;
   int                        incX, lda, lda0, m1, mb, mb1, n, nb;
#define    gpmv0              Mjoin( PATL, gpmvUT_a1_x1_b1_y1  )
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_GetPartMVT( A, LDA, &mb, &nb );

   if( DIAG == AtlasNonUnit ) tpmv0 = Mjoin( PATL, tpmvUTN );
   else                       tpmv0 = Mjoin( PATL, tpmvUTU );

   mb1 = N - ( m1 = ( ( N - 1 ) / mb ) * mb ); incX = (mb SHIFT);
   lda = LDA; MUpnext( m1, A, lda );
   X += (m1 SHIFT); x0 = X;
   tpmv0( mb1, A, lda, X ); X -= incX;

   lda0 = lda; A0 = (TYPE *)(A) - incX;

   for( n = mb1; n < N; n += mb, X -= incX, x0 -= incX )
   {
      gpmv0( n, mb, one, A0, lda0, X, 1, one, x0, 1 );
      MUpprev( mb, A, lda ); tpmv0( mb, A, lda, X );

      lda0 = lda; A0 = (TYPE *)(A) - incX;
   }
/*
 * End of Mjoin( PATL, tpmvUT )
 */
}
