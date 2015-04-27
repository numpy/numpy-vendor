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

#ifdef TREAL
#define    MinN_sprL    16
#else
#define    MinN_sprL     8
#endif

void Mjoin( PATL, sprL )
(
   const int                  N,
   const TYPE                 * X,
   const TYPE                 * Y,
   const int                  INCY,
   TYPE                       * A,
   const int                  LDA
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, sprL ) performs the symmetric rank 1 operation
 *
 *    A := alpha * x * x' + A,
 *
 * where  alpha is a scalar, x is an n-element vector and A is an n by n
 * symmetric matrix, supplied in packed form.
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
   register TYPE              t0;
   int                        i, iaij, j, jaj, jy, lda = LDA, n1, n1s, n2;
#else
   const TYPE                 one[2] = { ATL_rone, ATL_rzero };
   register TYPE              t0_i, t0_r;
   const int                  incy2 = 2 * INCY;
   int                        i, iaij, ix, j, jaj, jy, lda = LDA,
                              lda2 = (LDA << 1), n1, n1s, n2;
#endif
#ifdef TREAL
#define   gpr      Mjoin( PATL, gpr1L_a1_x1_yX  )
#else
#define   gpr      Mjoin( PATL, gpr1uL_a1_x1_yX )
#endif
/* ..
 * .. Executable Statements ..
 *
 */
   if( N <= MinN_sprL )
   {
#ifdef TREAL
      for( j = 0, jaj = 0, jy = 0; j < N; j++, jy += INCY )
      {
         t0 = Y[jy];
         for( i = j, iaij = jaj; i < N; i++, iaij++ ) { A[iaij] += X[i] * t0; }
         jaj += lda; lda -= 1;
      }
#else
      for( j = 0, jaj = 0, jy = 0; j < N; j++, jy += incy2, X += 2 )
      {
         t0_r = Y[jy]; t0_i = Y[jy+1];
         for( i = j, iaij = jaj, ix = 0; i < N; i++, iaij += 2, ix += 2 )
         {
            A[iaij  ] += X[ix] * t0_r - X[ix+1] * t0_i;
            A[iaij+1] += X[ix] * t0_i + X[ix+1] * t0_r;
         }
         jaj += lda2; lda2 -= 2;
      }
#endif
   }
   else
   {
      n2 = N - ( n1 = ( N >> 1 ) ); n1s = (n1 SHIFT);
      Mjoin( PATL, sprL )( n1, X, Y, INCY, A, lda ); X += n1s;
      gpr(  n2, n1, one, X, 1, Y, INCY, A + n1s, lda ); Y += n1s * INCY;
      MLpnext( n1, A, lda ); Mjoin( PATL, sprL )( n2, X, Y, INCY, A, lda );
   }
/*
 * End of Mjoin( PATL, sprL )
 */
}
