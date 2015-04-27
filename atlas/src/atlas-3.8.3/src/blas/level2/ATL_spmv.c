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

void Mjoin( PATL, spmv )
(
   const enum ATLAS_UPLO      UPLO,
   const int                  N,
   const SCALAR               ALPHA,
   const TYPE                 * A,
   const TYPE                 * X,
   const int                  INCX,
   const SCALAR               BETA,
   TYPE                       * Y,
   const int                  INCY
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, spmv ) performs the matrix-vector operation
 *
 *    y := alpha * A * x + beta * y,
 *
 * where alpha and beta are scalars, x and y are n-element vectors and A
 * is an n by n symmetric matrix, supplied in packed form.
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
   void                       (*gpmv0)( const int, const int, const SCALAR,
                              const TYPE *, const int, const TYPE *, const int,
                              const SCALAR, TYPE *, const int );
   void                       (*gpmv1)( const int, const int, const SCALAR,
                              const TYPE *, const int, const TYPE *, const int,
                              const SCALAR, TYPE *, const int );
   void                       (*gpmvN)( const int, const int, const SCALAR,
                              const TYPE *, const int, const TYPE *, const int,
                              const SCALAR, TYPE *, const int );
#ifdef TREAL
   TYPE                       alphaY, beta0;
#define one                   ATL_rone
#define zero                  ATL_rzero
#else
   const TYPE                 * alphaY, * beta0;
   const TYPE                 one [2] = { ATL_rone,  ATL_rzero },
                              zero[2] = { ATL_rzero, ATL_rzero };
#endif
   void                       * vx = NULL, * vy = NULL;
   TYPE                       * A0, * A1, * x, * x0, * x1, * y, * y00, * y0,
                              * y1;
   int                        incXY, incXY1, j, jb, lda, lda0, lda1, mb, mb1,
                              n, nb;
/* ..
 * .. Executable Statements ..
 *
 */
   if( N == 0 ) return;

   if( SCALAR_IS_ZERO( ALPHA ) )
   {
      if( !( SCALAR_IS_ONE( BETA ) ) ) Mjoin( PATL, scal )( N, BETA, Y, INCY );
      return;
   }

   if( ( INCX != 1 ) || ( ( INCY == 1 ) && !( SCALAR_IS_ONE( ALPHA ) ) ) )
   {
      vx = (void *)malloc( ATL_Cachelen + ATL_MulBySize( N ) );
      ATL_assert( vx ); x = ATL_AlignPtr( vx );
      Mjoin( PATL, cpsc )( N, ALPHA, X, INCX, x, 1 );
      alphaY = one;
   }
   else { x = (TYPE *)(X); alphaY = ALPHA; }

   if( ( INCY != 1 ) || !( SCALAR_IS_ONE( alphaY ) ) )
   {
      vy = malloc( ATL_Cachelen + ATL_MulBySize( N ) );
      ATL_assert( vy ); y00 = y = ATL_AlignPtr( vy );
      beta0 = zero;
   }
   else { y00 = y = (TYPE *)(Y); beta0 = BETA; }

   ATL_GetPartSPMV( A, N, &mb, &nb );

   mb1 = N - ( ( N - 1 ) / mb ) * mb; incXY1 = (nb SHIFT);

   if( UPLO == AtlasUpper )
   {
      if(      SCALAR_IS_ZERO( beta0 ) )
         gpmv0 = Mjoin( PATL, gpmvUT_a1_x1_b0_y1 );
      else if( SCALAR_IS_ONE ( beta0 ) )
         gpmv0 = Mjoin( PATL, gpmvUT_a1_x1_b1_y1 );
      else
         gpmv0 = Mjoin( PATL, gpmvUT_a1_x1_bX_y1 );
      gpmv1 = Mjoin( PATL, gpmvUT_a1_x1_b1_y1 );
      gpmvN = Mjoin( PATL, gpmvUN_a1_x1_b1_y1 );

      lda = 1; lda0 = lda; A0 = (TYPE *)(A); MUpnext( mb, A0, lda0 );
      incXY = (mb SHIFT); x0 = x + incXY; y0 = y + incXY;

      for( n = N - mb; n > 0; n -= mb, x0 += incXY, x += incXY,
           y0 += incXY, y += incXY )
      {
         Mjoin( PATL, spmvU )( mb, A, lda, x, beta0, y );

         for( j = 0, lda1 = lda0, A1 = A0 - (mb SHIFT), x1 = x0, y1 = y0; j < n;
              j += nb, x1 += incXY1, y1 += incXY1 )
         {
            jb = n - j; jb = Mmin( jb, nb );
            gpmv0( jb, mb, one, A1, lda1, x,  1, beta0, y1, 1 );
            gpmvN( mb, jb, one, A1, lda1, x1, 1, one,   y,  1 );
            MUpnext( jb, A1, lda1 ); A1 -= (jb SHIFT);
         }
         beta0 = one; gpmv0 = gpmv1; lda = lda0; A = A0; MUpnext( mb, A0, lda0 );
      }
      Mjoin( PATL, spmvU )( mb1, A, lda, x, beta0, y );
   }
   else
   {
      if(      SCALAR_IS_ZERO( beta0 ) )
         gpmv0 = Mjoin( PATL, gpmvLT_a1_x1_b0_y1 );
      else if( SCALAR_IS_ONE ( beta0 ) )
         gpmv0 = Mjoin( PATL, gpmvLT_a1_x1_b1_y1 );
      else
         gpmv0 = Mjoin( PATL, gpmvLT_a1_x1_bX_y1 );
      gpmv1 = Mjoin( PATL, gpmvLT_a1_x1_b1_y1 );
      gpmvN = Mjoin( PATL, gpmvLN_a1_x1_b1_y1 );

      lda = N; lda0 = lda; A0 = (TYPE *)(A); MLpnext( N, A, lda );
      incXY = (mb SHIFT); x0 = x; y0 = y;

      for( n  = N - mb, x += ((N-mb) SHIFT), y += ((N-mb) SHIFT); n > 0;
           n -= mb, x -= incXY, y -= incXY )
      {
         MLpprev( mb, A, lda );
         Mjoin( PATL, spmvL )( mb, A, lda, x, beta0, y );

         for( j = 0, lda1 = lda0, A1 = A0 + (n SHIFT), x1 = x0, y1 = y0; j < n;
              j += nb, x1 += incXY1, y1 += incXY1 )
         {
            jb = n - j; jb = Mmin( jb, nb );
            gpmv0( jb, mb, one, A1, lda1, x,  1, beta0, y1, 1 );
            gpmvN( mb, jb, one, A1, lda1, x1, 1, one,   y,  1 );
            MLpnext( jb, A1, lda1 ); A1 -= (jb SHIFT);
         }
         beta0 = one; gpmv0 = gpmv1;
      }
      Mjoin( PATL, spmvL )( mb1, A0, lda0, x0, beta0, y0 );
   }

   if( vx ) free( vx );
   if( vy )
   { Mjoin( PATL, axpby )( N, alphaY, y00, 1, BETA, Y, INCY ); free( vy ); }
/*
 * End of Mjoin( PATL, spmv )
 */
}
