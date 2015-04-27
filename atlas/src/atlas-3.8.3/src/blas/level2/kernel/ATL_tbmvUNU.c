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

#ifdef TREAL
#define    MinN_tbmvUNU      16
#else
#define    MinN_tbmvUNU       8
#endif

void Mjoin( PATL, tbmvUNU )
(
   const int                  N,
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
 * Mjoin( PATL, tbmvUNU ) performs the following matrix-vector operation
 *
 *    x := A * x,
 *
 * where x is an n-element vector and  A  is an n by n  unit, upper tri-
 * angular band matrix, with ( k + 1 ) diagonals.
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
#ifdef TREAL
#define    one                ATL_rone
#define    lda2               LDA
#else
   const TYPE                 one[2] = { ATL_rone, ATL_rzero };
   const int                  lda2   = ( LDA << 1 );
#endif
   int                        ia, kl, ku, ma, na, n1, n1s, n2;
/* ..
 * .. Executable Statements ..
 *
 */
   if( N <= MinN_tbmvUNU )
   { Mjoin( PATL, reftbmvUNU )( N, K, A, LDA, X, 1 ); }
   else
   {
      n2 = N - ( n1 = ( N >> 1 ) ); n1s = (n1 SHIFT);
      Mjoin( PATL, tbmvUNU )( n1, K, A, LDA, X );
      ia = n1 - K; ia = Mmax( ia, 0 ); ma = n1 - ia;    na = Mmin( n2, K );
      kl = ma - 1; kl = Mmax( kl, 0 ); ku = K - 1 - kl; ku = Mmax( ku, 0 );
      Mjoin( PATL, gbmvN_a1_x1_b1_y1 )( ma, na, kl, ku, one, A+n1*lda2, LDA,
                                        X+n1s, 1, one, X+(ia SHIFT), 1 );
      Mjoin( PATL, tbmvUNU )( n2, K, A+n1*lda2, LDA, X+n1s );
   }
/*
 * End of Mjoin( PATL, tbmvUNU )
 */
}
