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
#define    MinN_hprU    16
#else
#define    MinN_hprU     8
#endif

/*
 * Workaround for icc bug on IA64Itan2
 */
#ifdef ATL_IntelIccBugs
#pragma optimize("", off);
#endif
void Mjoin( PATL, hprU )
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
 * Mjoin( PATL, hprU ) performs the Hermitian rank 1 operation
 *
 *    A := alpha* x * conjg( x' ) + A,
 *
 * where  alpha is a real scalar, x is an n-element vector and A is an n
 * by n Hermitian matrix, supplied in packed form.
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
                              lda2 = (LDA SHIFT), n1, n1s, n2;
#endif
#define   gpr      Mjoin( PATL, gpr1cU_a1_x1_yX )
/* ..
 * .. Executable Statements ..
 *
 */
   if( N <= MinN_hprU )
   {
      for( j = 0, jaj = 0, jy = 0; j < N; j++, jy += incy2 )
      {
         t0_r = Y[jy]; t0_i = -Y[jy+1];
         for( i = 0, iaij = jaj, ix = 0; i < j; i++, iaij += 2, ix += 2 )
         {
            A[iaij  ] += X[ix] * t0_r - X[ix+1] * t0_i;
            A[iaij+1] += X[ix] * t0_i + X[ix+1] * t0_r;
         }
         A[iaij  ] += X[ix] * t0_r - X[ix+1] * t0_i;
         A[iaij+1]  = ATL_rzero;
         jaj += lda2; lda2 += 2;
      }
   }
   else
   {
      n2 = N - ( n1 = ( N >> 1 ) ); n1s = (n1 SHIFT);
      Mjoin( PATL, hprU )( n1, X, Y, INCY, A, lda );
      MUpnext( n1, A, lda ); Y += n1s * INCY;
      gpr( n1, n2, one, X, 1, Y, INCY, A - n1s, lda ); X += n1s;
      Mjoin( PATL, hprU )( n2, X, Y, INCY, A, lda );
   }
/*
 * End of Mjoin( PATL, hprU )
 */
}
