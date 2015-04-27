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

void Mjoin( PATL, spr2 )
(
   const enum ATLAS_UPLO      UPLO,
   const int                  N,
   const SCALAR               ALPHA,
   const TYPE                 * X,
   const int                  INCX,
   const TYPE                 * Y,
   const int                  INCY,
   TYPE                       * A
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, spr2 ) performs the symmetric rank 2 operation
 *
 *    A := alpha * x * y' + alpha * y * x' + A,
 *
 * where  alpha is a scalar, x and y are n-element vectors and A is an n
 * by n symmetric matrix, supplied in packed form.
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
   TYPE                       * A0, * A1, * x, * x0, * x1, * y, * y0, * y1;
   void                       * vx = NULL, * vy = NULL;
   int                        incA0, incXY, incXY1, j, jb, lda, lda0, lda1,
                              mb, mb1, n, nb;
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

   if( INCX != 1 )
   {
      vx = (TYPE *)malloc( ATL_Cachelen + ATL_MulBySize( N ) );
      ATL_assert( vx ); x = ATL_AlignPtr( vx );
      Mjoin( PATL, cpsc )( N, ALPHA, X, INCX, x, 1 );
   }
   else { x = (TYPE *)(X); }

   if( ( ( vx != NULL ) && ( INCY != 1 ) ) ||
       ( ( vx == NULL ) && ( !( SCALAR_IS_ONE( ALPHA ) ) || ( INCY != 1 ) ) ) )
   {
      vy = (TYPE *)malloc( ATL_Cachelen + ATL_MulBySize( N ) );
      ATL_assert( vy ); y = ATL_AlignPtr( vy );
      if( ( vx != NULL ) && ( INCY != 1 ) )
      {
         Mjoin( PATL, copy )( N,        Y, INCY, y, 1 );
      }
      else
      {
         Mjoin( PATL, cpsc )( N, ALPHA, Y, INCY, y, 1 );
      }
   }
   else { y = (TYPE *)(Y); }

   ATL_GetPartP1( A, LDA, mb, nb );

   mb1    = N - ( ( N - 1 ) / mb ) * mb;
   incXY  = (mb SHIFT);
   incXY1 = (nb SHIFT);

   if( UPLO == AtlasLower )
   {
      lda = N; Mjoin( PATL, spr2L )( mb1, x, y, A, lda );

      lda0 = lda; A0 = (TYPE *)(A); MLpnext( mb1, A, lda );
      x0 = x; x += (mb1 SHIFT); y0 = y; y += (mb1 SHIFT);

      for( n = mb1; n < N; n += mb, x += incXY, y += incXY )
      {
         for( j  = 0, lda1 = lda0, A1 = A0 + (n SHIFT), x1 = x0, y1 = y0;
              j < n; j += nb, x1 += incXY1, y1 += incXY1 )
         {
            jb = n - j; jb = Mmin( jb, nb );
            gpr1L( mb, jb, one, x, 1, y1, 1, A1, lda1 );
            gpr1L( mb, jb, one, y, 1, x1, 1, A1, lda1 );
            MLpnext( jb, A1, lda1 ); A1 -= (jb SHIFT);
         }
         Mjoin( PATL, spr2L )( mb, x, y, A, lda ); MLpnext( mb, A, lda );
      }
   }
   else
   {
      incA0 = incXY; lda0 = lda = 1; A0 = (TYPE *)(A); MUpnext( mb, A0, lda0 );

      for( n  = N - mb,         x0  = x + incXY,         y0  = y + incXY; n > 0;
           n -= mb, x += incXY, x0 += incXY, y += incXY, y0 += incXY )
      {
         Mjoin( PATL, spr2U )( mb, x, y, A, lda );

         for( j = 0, lda1 = lda0, A1 = A0 - incA0, x1 = x0, y1 = y0;
              j < n; j += nb, x1 += incXY1, y1 += incXY1 )
         {
            jb = n - j; jb = Mmin( jb, nb );
            gpr1U(  mb, jb, one, x, 1, y1, 1, A1, lda1 );
            gpr1U(  mb, jb, one, y, 1, x1, 1, A1, lda1 );
            MUpnext( jb, A1, lda1 ); A1 -= (jb SHIFT);
         }
         lda = lda0; A = A0; MUpnext( mb, A0, lda0 );
      }
      Mjoin( PATL, spr2U )( mb1, x, y, A, lda );
   }

   if( vx ) free( vx );
   if( vy ) free( vy );
/*
 * End of Mjoin( PATL, spr2 )
 */
}
