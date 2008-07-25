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

void Mjoin( PATL, tpsvUN )
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
 * Mjoin( PATL, tpsvUN ) solves the following triangular system of equations
 *
 *    A * x = b,
 *
 * where b and x are n-element vectors and  A  is an n by n unit or non-
 * unit, upper triangular matrix supplied in packed form.
 *
 * No test for  singularity  or  near-singularity  is included  in  this
 * routine. Such tests must be performed before calling this routine.
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
   void                       (*tpsv0)( const int, const TYPE *, const int,
                              TYPE * );
#ifdef TREAL
#define    none               ATL_rnone
#define    one                ATL_rone
#else
   TYPE                       none[2] = { ATL_rnone, ATL_rzero },
                              one [2] = { ATL_rone,  ATL_rzero };
#endif
#ifdef ATL_AXPYMV
   TYPE                       * x0;
   int                        incX, lda = LDA, m1, mb, mb1, n, nb;
#else
   TYPE                       * A0, * x0;
   int                        incX, lda = LDA, m, mb, nb, nb1;
#endif
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_GetPartMVN( A, N, &mb, &nb );

   if( DIAG == AtlasNonUnit ) tpsv0 = Mjoin( PATL, tpsvUNN );
   else                       tpsv0 = Mjoin( PATL, tpsvUNU );

#ifdef ATL_AXPYMV
   mb1  = N - ( m1 = ( ( N - 1 ) / mb ) * mb ); incX = (mb SHIFT);
   X += (m1 SHIFT); x0 = X; MUpnext( m1, A, lda );
   tpsv0( mb1, A, lda, X ); X -= incX;

   for( n  = mb1; n < N; n += mb, X -= incX, x0 -= incX )
   {
      Mjoin( PATL, gpmv )( AtlasUpper, AtlasNoTrans, mb, n, none,
                           A - incX, lda, x0, 1, one, X, 1 );
      MUpprev( mb, A, lda ); tpsv0( mb, A, lda, X );
   }
#else
   nb1  = N - ( ( N - 1 ) / nb ) * nb; incX = (nb SHIFT); x0 = X;
   A0 = (TYPE *)(A); MUpnext( N-nb, A, lda );

   for( m = N - nb, X += ((N-nb) SHIFT); m > 0; m -= nb, X -= incX )
   {
      tpsv0( nb, A, lda, X );
      Mjoin( PATL, gpmv )( AtlasUpper, AtlasNoTrans, m, nb, none,
                           A - (m SHIFT), lda, X, 1, one, x0, 1 );
      MUpprev( nb, A, lda );
   }
   tpsv0( nb1, A0, LDA, x0 );
#endif
/*
 * End of Mjoin( PATL, tpsvUN )
 */
}
