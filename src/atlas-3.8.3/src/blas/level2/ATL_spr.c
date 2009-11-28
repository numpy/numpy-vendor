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
#include "atlas_level1.h"
#include "atlas_kernel2.h"
#include "atlas_lvl2.h"
#include "atlas_mv.h"
#include "atlas_r1.h"

#include "atlas_reflvl2.h"          /* temporary for building purposes */
#include "atlas_reflevel2.h"        /* used for gbmv, gpmv and gpr.    */

void Mjoin( PATL, spr )
(
   const enum ATLAS_UPLO      UPLO,
   const int                  N,
   const SCALAR               ALPHA,
   const TYPE                 * X,
   const int                  INCX,
   TYPE                       * A
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, spr ) performs the symmetric rank 1 operation
 *
 *    A := alpha * x * x' + A,
 *
 * where  alpha is a scalar, x is an n-element vector and A is an n by n
 * symmetric matrix, supplied in packed form.
 *
 * This is a blocked version of the algorithm.  For a more detailed des-
 * cription of  the arguments of this function, see the reference imple-
 * mentation in the ATLAS/src/blas/reference directory.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
#ifdef TREAL
#define    one                ATL_rone
#else
   const TYPE                 one[2] = { ATL_rone, ATL_rzero };
#endif
   void                       * vx = NULL;
   TYPE                       * x, * y, * A0, * y0;
   int                        incA0, incX, incY, incy, lda, lda0, mb, mb1,
                              n, nb;
#ifdef TREAL
#define   gpr1L   Mjoin( PATL, gpr1L_a1_x1_yX  )
#define   gpr1U   Mjoin( PATL, gpr1U_a1_x1_yX  )
#else
#define   gpr1L   Mjoin( PATL, gpr1uL_a1_x1_yX )
#define   gpr1U   Mjoin( PATL, gpr1uU_a1_x1_yX )
#endif
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( N == 0 ) || ( SCALAR_IS_ZERO( ALPHA ) ) ) return;

   if( ( INCX == 1 ) && ( SCALAR_IS_ONE( ALPHA ) ) )
   {
      x = y = (TYPE *)(X); incy = INCX;
   }
   else
   {
      vx = (void *)malloc( ATL_Cachelen + ATL_MulBySize( N ) );
      ATL_assert( vx ); x = ATL_AlignPtr( vx );
      Mjoin( PATL, cpsc )( N, ALPHA, X, INCX, x, 1 );
      y = (TYPE *)(X); incy = INCX;
   }

   ATL_GetPartP1( A, N, mb, nb );

   mb1  = N - ( ( N - 1 ) / mb ) * mb;

   if( UPLO == AtlasLower )
   {
      lda = lda0 = N; incA0 = ( incX = (mb SHIFT) );

      Mjoin( PATL, sprL )( mb1, x, y, incy, A, lda );
      A0 = (TYPE *)( A + (mb1 SHIFT) ); MLpnext( mb1, A, lda );
      x += (mb1 SHIFT);

      for( n = mb1; n < N; n += mb, A0 += incA0, x += incX )
      {
         gpr1L( mb, n, one, x, 1, y, incy, A0, lda0 );
         Mjoin( PATL, sprL )( mb, x, y + n*(incy SHIFT), incy, A, lda );
         MLpnext( mb, A, lda );
      }
   }
   else
   {
      incA0 = ( incX = (mb SHIFT) ); incY = incX * incy;
      lda0 = lda = 1; A0 = (TYPE *)(A); MUpnext( mb, A0, lda0 );

      for( n  = N - mb, y0 = y + incY; n > 0;
           n -= mb, x += incX, y += incY, y0 += incY )
      {
         Mjoin( PATL, sprU )( mb, x, y, incy, A, lda );
         gpr1U(  mb, n, one, x, 1, y0, incy, A0 - incA0, lda0 );
         lda = lda0; A = A0; MUpnext( mb, A0, lda0 );
      }
      Mjoin( PATL, sprU )( mb1, x, y, incy, A, lda );
   }

   if( vx ) free( vx );
/*
 * End of Mjoin( PATL, spr )
 */
}
