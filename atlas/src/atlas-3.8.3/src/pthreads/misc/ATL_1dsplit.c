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
#include "atlas_ptmisc.h"

DIM_1DSPLIT_T ATL_1dsplit
(
   const unsigned int         NT,
   const int                  N,
   const int                  NB,
   unsigned int               * NT1,
   unsigned int               * NT2,
   int                        * N1,
   int                        * N2,
   double                     * RATIO
)
{
/*
 * Purpose
 * =======
 *
 * ATL_1dsplit splits  a  problem  dimension  N  within NT threads. This
 * function returns whether or not the problem size N should be split.
 *
 * Arguments
 * =========
 *
 * NT      (input)                       const unsigned int
 *         On entry, NT specifies the initial total number of threads.
 *
 * N       (input)                       const int
 *         On entry, N specifies the problem dimension to be split.
 *
 * NB      (input)                       const int
 *         On entry,  NB  specifies the blocking factor to be used along
 *         this dimension N.
 *
 * NT1     (output)                      unsigned int *
 *         On exit,  NT1  specifies the number of threads to be used for
 *         the first part of the problem size N1.
 *
 * NT2     (output)                      unsigned int *
 *         On exit,  NT2  specifies the number of threads to be used for
 *         the second part of the problem size N2.
 *
 * N1      (output)                      int *
 *         On exit,  N1  specifies  the length of the problem size to be
 *         run on the NT1 threads.
 *
 * N2      (output)                      int *
 *         On exit,  N2  specifies  the length of the problem size to be
 *         run on the NT2 threads.
 *
 * RATIO   (output)                      double *
 *         On exit, RATIO specifies the relative work distribution ratio
 *         namely | ( N1 / NT1 ) - ( N2 / NT2 ) | / ( N / NT ).
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   double                     n, nt;
   DIM_1DSPLIT_T              split = Atlas1dNoSplit;
   unsigned int               nt1, nt2, ntm1;
   int                        n1, n2, nblks, nbm1;
/* ..
 * .. Executable Statements ..
 *
 */
   nblks = ( N + ( nbm1 = NB - 1 ) ) / NB;
   if( ( nblks < 2 ) || ( NT < 2 ) ) return( split );

   n1 = ( ( ( ( N + ( ntm1 = NT - 1 ) ) / NT ) + nbm1 ) / NB ) * NB;
   n1 = Mmin( n1, N );
   n1 = Mmax( n1, 1 );
   n2 = N - n1;

   n  = (double)(N);
   nt = (double)(NT);

   nt1 = (int)floor( ( (double)(n1) / n ) * nt + 0.5 );
   nt2 = NT - ( nt1 = Mmin( nt1, ntm1 ) );

   *N1    = n1;
   *N2    = n2;
   *NT1   = nt1;
   *NT2   = nt2;
   *RATIO = ( (double)(n1) / (double)(nt1) ) - ( (double)(n2) / (double)(nt2) );
   *RATIO = Mabs( *RATIO ) / ( n / nt );

   return( Atlas1dSplit );
/*
 * End of ATL_1dsplit
 */
}
