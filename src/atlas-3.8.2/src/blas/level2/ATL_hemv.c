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

void Mjoin( PATL, hemv )
(
   const enum ATLAS_UPLO      UPLO,
   const int                  N,
   const SCALAR               ALPHA,
   const TYPE                 * A,
   const int                  LDA,
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
 * Mjoin( PATL, hemv ) performs the matrix-vector operation
 *
 *    y := alpha * A * x + beta * y,
 *
 * where alpha and beta are scalars, x and y are n-element vectors and A
 * is an n by n Hermitian matrix.
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
   void                       (*gemv0)( const int, const int, const SCALAR,
                              const TYPE *, const int, const TYPE *, const int,
                              const SCALAR, TYPE *, const int );
   void                       (*gemv1)( const int, const int, const SCALAR,
                              const TYPE *, const int, const TYPE *, const int,
                              const SCALAR, TYPE *, const int );
   void                       (*gemvN)( const int, const int, const SCALAR,
                              const TYPE *, const int, const TYPE *, const int,
                              const SCALAR, TYPE *, const int );
#ifdef TREAL
   TYPE                       alphaY, beta0;
#define lda2                  LDA
#define one                   ATL_rone
#define zero                  ATL_rzero
   void                       * vx = NULL, * vy = NULL;
   TYPE                       * x, * y;
   TYPE                       * A0, * A1, * x0, * x1, * y00, * y0, * y1;
#else
   const int                  lda2 = ( LDA SHIFT );
   const TYPE                 * alphaY, * beta0;
   const TYPE                 one [2] = { ATL_rone,  ATL_rzero },
                              zero[2] = { ATL_rzero, ATL_rzero };
   void                       * vx = NULL, * vy = NULL;
   TYPE                       * x, * y;
   TYPE                       * A0, * A1, * x0, * x1, * y00, * y0, * y1;
#endif
   int                        incA, incA1, incXY, incXY1, j, jb, mb, mb1, n, nb;
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

   if(      SCALAR_IS_ZERO( beta0 ) ) gemv0 = Mjoin( PATL, gemvC_a1_x1_b0_y1 );
   else if( SCALAR_IS_ONE ( beta0 ) ) gemv0 = Mjoin( PATL, gemvC_a1_x1_b1_y1 );
   else                               gemv0 = Mjoin( PATL, gemvC_a1_x1_bX_y1 );
   gemv1 = Mjoin( PATL, gemvC_a1_x1_b1_y1 );
   gemvN = Mjoin( PATL, gemvS_a1_x1_b1_y1 );

   ATL_GetPartSYMV( A, LDA, &mb, &nb );

   mb1   = N - ( ( N - 1 ) / mb ) * mb;
   incA1 = nb * lda2;  incXY1 = (nb SHIFT);

   if( UPLO == AtlasUpper )
   {
      incA  = ( incXY = (mb SHIFT) ) + mb * lda2;
      A0 = (TYPE *)(A) + mb * lda2; x0 = x + incXY; y0 = y + incXY;

      for( n  = N - mb; n > 0; n -= mb, A0 += incA, A += incA, x0 += incXY,
           x += incXY, y0 += incXY, y += incXY )
      {
         Mjoin( PATL, hemvU )( mb, A, LDA, x, beta0, y );

         for( j  =  0, A1 = A0,     x1  = x0,     y1  = y0;     j < n;
              j += nb, A1 += incA1, x1 += incXY1, y1 += incXY1 )
         {
            jb = n - j; jb = Mmin( jb, nb );
            gemv0( jb, mb, one, A1, LDA, x,  1, beta0, y1, 1 );
            gemvN( mb, jb, one, A1, LDA, x1, 1, one,   y,  1 );
         }
         beta0 = one; gemv0 = gemv1;
      }
      Mjoin( PATL, hemvU )( mb1, A, LDA, x, beta0, y );
   }
   else
   {
      incA = incXY = (mb SHIFT);
      A0 = (TYPE *)(A); x0 = x; y0 = y;

      for( n  = N - mb, A += ((N-mb) SHIFT), x += ((N-mb) SHIFT),
           y += ((N-mb) SHIFT); n > 0; n -= mb, A -= incA, x -= incXY,
           y -= incXY )
      {
         Mjoin( PATL, hemvL )( mb, A+n*lda2, LDA, x, beta0, y );

         for( j  =  0, A1  = (TYPE *)(A), x1  = x0,     y1  = y0;      j < n;
              j += nb, A1 += incA1,       x1 += incXY1, y1 += incXY1 )
         {
            jb = n - j; jb = Mmin( jb, nb );
            gemv0( jb, mb, one, A1, LDA, x,  1, beta0, y1, 1 );
            gemvN( mb, jb, one, A1, LDA, x1, 1, one,   y,  1 );
         }
         beta0 = one; gemv0 = gemv1;
      }
      Mjoin( PATL, hemvL )( mb1, A0, LDA, x0, beta0, y0 );
   }

   if( vx ) free( vx );
   if( vy )
   { Mjoin( PATL, axpby )( N, alphaY, y00, 1, BETA, Y, INCY ); free( vy ); }
/*
 * End of Mjoin( PATL, hemv )
 */
}
