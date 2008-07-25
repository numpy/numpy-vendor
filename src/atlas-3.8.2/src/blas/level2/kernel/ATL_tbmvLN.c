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

void Mjoin( PATL, tbmvLN )
(
   const enum ATLAS_DIAG      DIAG,
   const int                  N,         /* N > 0 assumed */
   const int                  K,
   const TYPE                 * A,
   const int                  LDA,
   TYPE                       * X
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, tbmvLN ) performs the following matrix-vector operation
 *
 *    x := A * x,
 *
 * where x is an n-element vector and  A is an n by n unit or  non-unit,
 * lower triangular band matrix.
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
   void                       (*tbmv0)( const int, const int, const TYPE *,
                              const int, TYPE * );
#ifdef TREAL
#define    lda2               LDA
#define    one                ATL_rone
#else
   TYPE                       one[2] = { ATL_rone, ATL_rzero };
   const int                  lda2   = ( LDA SHIFT );
#endif
#ifdef ATL_AXPYMV
   int                        ja, jan, kl, ku, ma, mb, mb1, na, nb;
#else
   int                        ja, jan, kl, ku, m, ma, mb, n1, na, nb, nb1;
#endif
#define    gbmv0              Mjoin( PATL, gbmvN_a1_x1_b1_y1  )
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_GetPartMVN( A, LDA, &mb, &nb );

   if( DIAG == AtlasNonUnit ) tbmv0 = Mjoin( PATL, tbmvLNN );
   else                       tbmv0 = Mjoin( PATL, tbmvLNU );

#ifdef ATL_AXPYMV
   mb1  = N - ( ( N - 1 ) / mb ) * mb;

   for( jan = N - mb; jan > 0; jan -= mb )
   {
      ja = jan - K; ja = Mmax( ja, 0 ); ma = Mmin( mb, K );
      ku = na = jan - ja;  kl = K - ku; kl = Mmax( kl, 0 );
      tbmv0( mb, K, A+jan*lda2, LDA, X+(jan SHIFT) );
      gbmv0( ma, na, kl, ku, one, A+ja*lda2, LDA, X+(ja SHIFT), 1, one,
             X+(jan SHIFT), 1 );
   }
   tbmv0( mb1, K, A, LDA, X );
#else
   nb1  = N - ( n1  = ( ( N - 1 ) / nb ) * nb );

   tbmv0( nb1, K, A + n1*lda2, LDA, X+(n1 SHIFT) );

   for( m = nb1, jan = n1; m < N; m += nb )
   {
      ja = nb - K; ja = jan - nb + Mmax( ja, 0 ); ma = Mmin( m, K );
      ku = na = jan - ja; kl = K - ku; kl = Mmax( kl, 0 );
      gbmv0( ma, na, kl, ku, one, A+ja*lda2, LDA, X+(ja SHIFT), 1, one,
             X+(jan SHIFT), 1 );
      jan -= nb; tbmv0( nb, K, A+jan*lda2, LDA, X+(jan SHIFT) );
   }
#endif
/*
 * End of Mjoin( PATL, tbmvLN )
 */
}
