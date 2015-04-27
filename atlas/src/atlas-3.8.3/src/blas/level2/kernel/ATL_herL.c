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
#define    MinN_herL    16
#else
#define    MinN_herL     8
#endif

/*
 * Workaround for icc bug on IA64Itan2
 */
#ifdef ATL_IntelIccBugs
#pragma optimize("", off);
#endif
void Mjoin( PATL, herL )
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
 * Mjoin( PATL, herL ) performs the Hermitian rank 1 operation
 *
 *    A := alpha* x * conjg( x' ) + A,
 *
 * where  alpha is a real scalar, x is an n-element vector and A is an n
 * by n Hermitian matrix.
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
   int                        i, j, jy, n1, n1s, n2;
#else
   const TYPE                 one[2] = { ATL_rone, ATL_rzero };
   register TYPE              t0_i, t0_r;
   const int                  incy2 = 2 * INCY, ldap12 = ( (LDA + 1) SHIFT);
   int                        i, iax, j, jy, n1, n1s, n2;
#endif
#define   ger      Mjoin( PATL, ger1c_a1_x1_yX )
/* ..
 * .. Executable Statements ..
 *
 */
   if( N <= MinN_herL )
   {
      for( j = 0, jy = 0; j < N; j++, jy += incy2, A += ldap12, X += 2 )
      {
         t0_r  = Y[jy]; t0_i = -Y[jy+1];
         (*A) += (*X) * t0_r - X[1] * t0_i;
         A[1]  = ATL_rzero;
         for( i = j+1, iax = 2; i < N; i++, iax += 2 )
         {
            A[iax  ] += X[iax] * t0_r - X[iax+1] * t0_i;
            A[iax+1] += X[iax] * t0_i + X[iax+1] * t0_r;
         }
      }
   }
   else
   {
      n2 = N - ( n1 = ( N >> 1 ) ); n1s = (n1 SHIFT);
      Mjoin( PATL, herL )( n1, X, Y, INCY, A, LDA ); X += n1s;
      ger( n2, n1, one, X, 1, Y, INCY, A + n1s, LDA );  Y += n1s * INCY;
      MLrnext( n1, A, LDA ); Mjoin( PATL, herL )( n2, X, Y, INCY, A, LDA );
   }
/*
 * End of Mjoin( PATL, herL )
 */
}
