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

void Mjoin( PATL, tpmvLN )
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
 * Mjoin( PATL, tpmvLN ) performs the following matrix-vector operation
 *
 *    x := A * x,
 *
 * where x is an n-element vector and  A is an n by n unit or  non-unit,
 * lower triangular matrix supplied in packed form.
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
#else
   TYPE                       one[2] = { ATL_rone, ATL_rzero };
#endif
#ifdef ATL_AXPYMV
   TYPE                       * A0, * x0;
   int                        incX, lda, lda0, mb, mb1, n, nb;
#else
   TYPE                       * x0;
   int                        incX, lda, m, mb, n1, nb, nb1;
#endif
   void                       (*tpmv0)( const int, const TYPE *, const int,
                              TYPE * );
#define    gpmv0              Mjoin( PATL, gpmvLN_a1_x1_b1_y1  )
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_GetPartMVN( A, LDA, &mb, &nb );

   if( DIAG == AtlasNonUnit ) tpmv0 = Mjoin( PATL, tpmvLNN );
   else                       tpmv0 = Mjoin( PATL, tpmvLNU );

#ifdef ATL_AXPYMV
   mb1  = N - ( ( N - 1 ) / mb ) * mb; incX = (mb SHIFT); x0 = X;
   lda = lda0 = LDA; A0 = (TYPE *)(A); MLpnext( N, A, lda );

   for( n = N - mb, X += ((N-mb) SHIFT); n > 0; n -= mb, X -= incX )
   {
      MLpprev( mb, A, lda ); tpmv0( mb, A, lda, X );
      gpmv0( mb, n, one, A0 + (n SHIFT), lda0, x0, 1, one, X, 1 );
   }
   tpmv0( mb1, A0, lda0, x0 );
#else
   nb1  = N - ( n1  = ( ( N - 1 ) / nb ) * nb );
   incX = (nb SHIFT); x0 = ( X += (n1 SHIFT) );
   lda = LDA; MLpnext( n1, A, lda ); tpmv0( nb1, A, lda, X ); x0 = X; X -= incX;
   for( m = nb1; m < N; m += nb, X -= incX, x0 -= incX )
   {
      MLpprev( nb, A, lda );
      gpmv0( m, nb, one, A + (nb SHIFT), lda, X, 1, one, x0, 1 );
      tpmv0( nb, A, lda, X );
   }
#endif
/*
 * End of Mjoin( PATL, tpmvLN )
 */
}
