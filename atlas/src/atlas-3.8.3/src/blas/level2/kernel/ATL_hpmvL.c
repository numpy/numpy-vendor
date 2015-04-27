/* ---------------------------------------------------------------------
 *
 * -- Automatically Tuned Linear Algebra Software (ATLAS)
 *    (C) Copyright 2000 All Rights Reserved
 *
 * -- ATLAS routine -- Version 3.2 -- December 25, 2000
 *
 * Author         : Antoine P. Petitet
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

void Mjoin( PATL, hpmvL )
(
   const int                  N,
   const TYPE                 * A,
   const int                  LDA,
   const TYPE                 * X,
   const SCALAR               BETA,
   TYPE                       * Y
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, hpmvL ) performs the matrix-vector operation
 *
 *    y := alpha * A * x + beta * y,
 *
 * where alpha and beta are scalars, x and y are n-element vectors and A
 * is an n by n Hermitian matrix, supplied in packed form.
 *
 * This is a  recursive  version of the  algorithm.  For a more detailed
 * description of  the arguments of this function, see the reference im-
 * plementation in the  ATLAS/src/blas/reference directory.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   void                       (*gpmvT)( const int, const int, const SCALAR,
                              const TYPE *, const int, const TYPE *, const int,
                              const SCALAR, TYPE *, const int );
   void                       (*gpmvN)( const int, const int, const SCALAR,
                              const TYPE *, const int, const TYPE *, const int,
                              const SCALAR, TYPE *, const int );
#ifdef TREAL
#define    one                ATL_rone
   TYPE                       beta0;
#else
   const TYPE                 * beta0, one[2] = { ATL_rone, ATL_rzero };
#endif
   TYPE                       * A0, * x0, * y0;
   int                        j, jb, jbs, lda = LDA, m, mb, nb;
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_GetPartSPMV( A, LDA, &mb, &nb );

   beta0 = BETA; x0 = (TYPE *)(X); y0 = (TYPE *)(Y);
   gpmvT = Mjoin( PATL, gpmvLC_a1_x1_b1_y1 );
   if(      SCALAR_IS_ZERO( beta0 ) ) gpmvN = Mjoin( PATL, gpmvLN_a1_x1_b0_y1 );
   else if( SCALAR_IS_ONE ( beta0 ) ) gpmvN = Mjoin( PATL, gpmvLN_a1_x1_b1_y1 );
   else                               gpmvN = Mjoin( PATL, gpmvLN_a1_x1_bX_y1 );

   for( j = 0; j < N; j += nb )
   {
      jb = N - j; jb = Mmin( jb, nb );
      Mjoin( PATL, refhpmvL )( jb, one, A, lda, X, 1, beta0, Y, 1 );

      if( ( m = N-j-jb ) != 0 )
      {
         jbs = (jb SHIFT); A0 = (TYPE *)(A) + jbs; X += jbs; Y += jbs;
         gpmvT( jb, m, one, A0, lda, X,  1, one,   y0, 1 );
         gpmvN( m, jb, one, A0, lda, x0, 1, beta0, Y,  1 );
         beta0 = one; gpmvN = Mjoin( PATL, gpmvLN_a1_x1_b1_y1 );
         MLpnext( jb, A, lda ); x0 = (TYPE *)(X); y0 = (TYPE *)(Y);
      }
   }
/*
 * End of Mjoin( PATL, hpmvL )
 */
}
