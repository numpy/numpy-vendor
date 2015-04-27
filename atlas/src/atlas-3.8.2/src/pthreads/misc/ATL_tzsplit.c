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

DIM_TZSPLIT_T ATL_tzsplit
(
   const enum ATLAS_UPLO      UPLO,
   const unsigned int         NT,
   const int                  M,
   const int                  N,
   const int                  K,
   const int                  NB,
   unsigned int               * NT1,
   unsigned int               * NT2,
   int                        * MNK1,
   int                        * MNK2
)
{
/*
 * Purpose
 * =======
 *
 * ATL_tzsplit splits a trapezoidal matrix into two matrices of the same
 * kind. The initial input matrix sizes are specified by M, N and K:
 *
 *              K       N                           N       K
 *             ______________                 _____________
 *            |   |          |               |            |\
 *          M |   |          |             K |            | \
 *            | __|__________|               |____________|__\
 *             \  |          |               |            |  |
 *          K   \ |          |             M |            |  |
 *               \|__________|               |____________|__|
 *
 *           UPLO = AtlasLower               UPLO = AtlasUpper
 *
 * This function returns which dimension should be split if any, as well
 * as the number of threads to be used for each submtrix.
 *
 * Arguments
 * =========
 *
 * UPLO    (input)                       const enum ATLAS_UPLO
 *         On entry,  UPLO  specifies  whether  the  matrix  is upper or
 *         lower trapezoidal.
 *
 * NT      (input)                       const unsigned int
 *         On entry, NT specifies the initial total number of threads.
 *
 * M       (input)                       const int
 *         On entry, M  specifies  the  number  of  complete rows of the
 *         trapezoidal matrix.
 *
 * N       (input)                       const int
 *         On entry, N  specifies  the number of complete columns of the
 *         trapezoidal matrix.
 *
 * K       (input)                       const int
 *         On entry, K specifies the size of the triangular submatrix.
 *
 * NB      (input)                       const int
 *         On entry,  NB  specifies the blocking factor to be used along
 *         this dimension N.
 *
 * NT1     (output)                      unsigned int *
 *         On exit,  NT1  specifies the number of threads to be used for
 *         the first part of the problem size MNK1.
 *
 * NT2     (output)                      unsigned int *
 *         On exit,  NT2  specifies the number of threads to be used for
 *         the second part of the problem size MNK2.
 *
 * MNK1    (output)                      int *
 *         On exit, MNK1 specifies  the length of the problem size to be
 *         run on the NT1 threads.
 *
 * MNK1    (output)                      int *
 *         On exit, MNK2 specifies  the length of the problem size to be
 *         run on the NT2 threads.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   int                        mnks[8], nts[8];
   double                     gacol=0.0, garow=0.0, k, k2, m, mnk1, mnk2,
                              n, nt, nt1, nt2, p, r, rcol=0.0, rrow=0.0,
                              s, tmnk1, tmnk2, work;
   DIM_TZSPLIT_T              split = AtlasTzNoSplit;
   int                        kblks, mblks, nblks, nbm1, ntm1;
/* ..
 * .. Executable Statements ..
 *
 */
   nbm1  = NB - 1;
   mblks = (M+nbm1) / NB; nblks = (N+nbm1) / NB; kblks = (K+nbm1) / NB;
   r     = (double)(mblks + kblks) * (double)(kblks + nblks);

   if( ( r < 4.0 ) || ( NT < 2 ) ) return( split );

   m     = (double)(M);   n     = (double)(N);   k     = (double)(K);
   k2    = k * k;         nt    = (double)(NT);  ntm1  = NT - 1;
   work  = 2.0 * ( ( k + m ) * n + ( k * m ) ) + k2;

   if( ( M > 0 ) || ( K > 0 ) )                       /* Column split */
   {
      s     = k + m;
      p     = k * m;
      gacol = ( 2.0 * nt * n * s  - work ) / ( nt * s );

      if( ( N > 1 ) && ( gacol >= 0.0 ) )
      {                                                    /* split N */
         mnk1   = ( ( N - (int)(gacol / 2.0) + nbm1 ) / NB ) * NB;
         mnk1   = Mmin( mnk1, N );
         mnk1   = Mmax( mnk1, 1 );
         mnk2   = N - mnk1;

         tmnk1  = (double)(mnk1);
         tmnk2  = (double)(mnk2);

         nt1    = (int)floor( ( ( 2 * s * tmnk1 ) / work ) * nt + 0.5 );
         nt2    = NT - ( nt1 = Mmin( nt1, ntm1 ) );

         rcol   = ( ( 2.0 *   s * tmnk1            ) / (double)(nt1) ) -
                  ( ( 2.0 * ( s * tmnk2 + p ) + k2 ) / (double)(nt2) );

         if( UPLO == AtlasLower )
         {
            mnks[0] = mnk1; mnks[1] = mnk2;
            nts [0] = nt1;  nts [1] = nt2;
         }
         else
         {
            mnks[0] = mnk2; mnks[1] = mnk1;
            nts [0] = nt2;  nts [1] = nt1;
         }
      }
      else
      {                                                    /* split K */
         mnk1   = ( ( (int)(s * ( 1.0 - sqrt( 1.0 + ( gacol / s ) ) )) +
                    nbm1 ) / NB ) * NB;
         mnk1   = Mmin( mnk1, K );
         mnk1   = Mmax( mnk1, 1 );
         mnk2   = K - mnk1;

         tmnk1  = (double)(mnk1);
         tmnk2  = (double)(mnk2);

         r      = work - ( tmnk2 * ( tmnk2 + ( 2.0 * m ) ) );
         nt1    = (int)floor( ( r / work ) * nt + 0.5 );
         nt2    = NT - ( nt1 = Mmin( nt1, ntm1 ) );

         rcol   = 2.0 * ( s * n + tmnk1 * ( m + tmnk2 ) ) +
                  ( tmnk1 * tmnk1 );
         rcol  /= (double)(nt1);
         rcol  -= ( tmnk2 * ( 2.0 * m + tmnk2 ) ) / (double)(nt2);

         if( UPLO == AtlasLower )
         {
            mnks[2] = mnk1; mnks[3] = mnk2;
            nts [2] = nt1;  nts [3] = nt2;
         }
         else
         {
            mnks[2] = mnk2; mnks[3] = mnk1;
            nts [2] = nt2;  nts [3] = nt1;
         }
      }
   }
                                                         /* Row split */
   if( ( N > 0 ) || ( K > 0 ) )
   {
      s     = k + n;
      p     = k * n;
      garow = ( 2.0 * nt * m * s - work ) / ( nt * s );

      if( ( M > 1 ) && ( garow >= 0.0 ) )
      {                                                    /* split M */
         mnk1   = ( ( M - (int)(garow / 2.0) + nbm1 ) / NB ) * NB;
         mnk1   = Mmin( mnk1, M );
         mnk1   = Mmax( mnk1, 1 );
         mnk2   = M - mnk1;

         tmnk1  = (double)(mnk1);
         tmnk2  = (double)(mnk2);

         nt1    = (int)floor( ( ( 2 * s * tmnk1 ) / work ) * nt + 0.5 );
         nt2    = NT - ( nt1 = Mmin( nt1, ntm1 ) );

         rrow   = ( ( 2.0 *   s * tmnk1            ) / (double)(nt1) ) -
                  ( ( 2.0 * ( s * tmnk2 + p ) + k2 ) / (double)(nt2) );

         if( UPLO == AtlasUpper )
         {
            mnks[4] = mnk1; mnks[5] = mnk2;
            nts [4] = nt1;  nts [5] = nt2;
         }
         else
         {
            mnks[4] = mnk2; mnks[5] = mnk1;
            nts [4] = nt2;  nts [5] = nt1;
         }
      }
      else
      {                                                    /* split K */
         mnk1   = ( ( (int)(s * ( 1.0 - sqrt( 1.0 + ( garow / s ) ) )) +
                      nbm1 ) / NB ) * NB;
         mnk1   = Mmin( mnk1, K );
         mnk1   = Mmax( mnk1, 1 );
         mnk2   = K - mnk1;

         tmnk1  = (double)(mnk1);
         tmnk2  = (double)(mnk2);

         r      = work - ( tmnk2 * ( tmnk2 + ( 2.0 * n ) ) );
         nt1    = (int)floor( ( r / work ) * nt + 0.5 );
         nt2    = NT - ( nt1 = Mmin( nt1, ntm1 ) );

         rrow   = 2.0 * ( s * m + tmnk1 * ( n + tmnk2 ) ) +
                  ( tmnk1 * tmnk1 );
         rrow  /= (double)(nt1);
         rrow  -= ( tmnk2 * ( 2.0 * n + tmnk2 ) ) / (double)(nt2);

         if( UPLO == AtlasUpper )
         {
            mnks[6] = mnk1; mnks[7] = mnk2;
            nts [6] = nt1;  nts [7] = nt2;
         }
         else
         {
            mnks[6] = mnk2; mnks[7] = mnk1;
            nts [6] = nt2;  nts [7] = nt1;
         }
      }
   }

   if( ( ( M > 0 ) || ( K > 0 ) ) && ( ( K > 0 ) || ( N > 0 ) ) )
   {
      if( Mabs( rcol ) <= Mabs( rrow ) )              /* Column split */
      {
         if( ( N > 1 ) && ( gacol >= 0.0 ) )
         {
            split = AtlasTzSplitNcol;
            *NT1  = nts [0]; *NT2  = nts [1];
            *MNK1 = mnks[0]; *MNK2 = mnks[1];
         }
         else
         {
            split = AtlasTzSplitKcol;
            *NT1  = nts [2]; *NT2  = nts [3];
            *MNK1 = mnks[2]; *MNK2 = mnks[3];
         }
      }
      else
      {
         if( ( M > 1 ) && ( garow >= 0.0 ) )
         {
            split = AtlasTzSplitMrow;
            *NT1  = nts [4]; *NT2  = nts [5];
            *MNK1 = mnks[4]; *MNK2 = mnks[5];
         }
         else
         {
            split = AtlasTzSplitKrow;
            *NT1  = nts [6]; *NT2  = nts [7];
            *MNK1 = mnks[6]; *MNK2 = mnks[7];
         }
      }
   }
   else if( ( M > 0 ) || ( K > 0 ) )
   {
      if( ( N > 1 ) && ( gacol >= 0.0 ) )
      {
         split = AtlasTzSplitNcol;
         *NT1  = nts [0]; *NT2  = nts [1];
         *MNK1 = mnks[0]; *MNK2 = mnks[1];
      }
      else
      {
         split = AtlasTzSplitKcol;
         *NT1  = nts [2]; *NT2  = nts [3];
         *MNK1 = mnks[2]; *MNK2 = mnks[3];
      }
   }
   else if( ( N > 0 ) || ( K > 0 ) )
   {
      if( ( M > 1 ) && ( garow >= 0.0 ) )
      {
         split = AtlasTzSplitMrow;
         *NT1  = nts [4]; *NT2  = nts [5];
         *MNK1 = mnks[4]; *MNK2 = mnks[5];
      }
      else
      {
         split = AtlasTzSplitKrow;
         *NT1  = nts [6]; *NT2  = nts [7];
         *MNK1 = mnks[6]; *MNK2 = mnks[7];
      }
   }
   else
   {
      split = AtlasTzNoSplit;
   }

   return( split );
/*
 * End of ATL_tzsplit
 */
}
