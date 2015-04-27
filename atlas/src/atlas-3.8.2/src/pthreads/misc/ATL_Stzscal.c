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

PT_TREE_T ATL_Stzscal
(
   const PT_MISC_TYPE_T       * PTYPE,
   const unsigned int         NODE,
   const unsigned int         THREADS,
   pthread_attr_t             * ATTR,
   const int                  NB,
   const enum ATLAS_UPLO      UPLO,
   const int                  M,
   const int                  N,
   const int                  K,
   const void                 * ALPHA,
   void                       * A,
   const int                  LDA
)
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  left, right, tree = NULL;
   PT_TZSCAL_ARGS_T           * a_scal = NULL;
   size_t                     size;
   DIM_TZSPLIT_T              split = AtlasTzNoSplit;
   unsigned int               next = (NODE << 1), nt1, nt2;
   int                        mnk1, mnk2;
/* ..
 * .. Executable Statements ..
 *
 */
   if( THREADS < 1 ) return( tree );

   size  = PTYPE->size;
   split = ATL_tzsplit( UPLO, THREADS, M, N, K, NB, &nt1, &nt2,
                        &mnk1, &mnk2 );

   if(      split == AtlasTzSplitMrow )
   {
      if( UPLO == AtlasLower )
      {
         left  = ATL_Stzscal( PTYPE, next+1, nt1, ATTR, NB, UPLO, mnk1, N, K,
                              ALPHA, A, LDA );
         right = ATL_Sgescal( PTYPE, next+2, nt2, ATTR, NB, mnk2, N+K, ALPHA,
                              Mvptm( A, K + mnk1, size ), LDA );
      }
      else
      {
         left  = ATL_Sgescal( PTYPE, next+1, nt1, ATTR, NB, mnk1, N+K, ALPHA,
                              A, LDA );
         right = ATL_Stzscal( PTYPE, next+2, nt2, ATTR, NB, UPLO, mnk2, N, K,
                              ALPHA, Mvptm( A, mnk1, size ), LDA );
      }
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else if( split == AtlasTzSplitKrow )
   {
      if( UPLO == AtlasLower )
      {
         left  = ATL_Stzscal( PTYPE, next+1, nt1, ATTR, NB, UPLO, 0, N,
                              mnk1, ALPHA, A, LDA );
         right = ATL_Stzscal( PTYPE, next+2, nt2, ATTR, NB, UPLO, M, N+mnk1,
                              mnk2, ALPHA, Mvptm( A, mnk1, size ), LDA );
      }
      else
      {
         left  = ATL_Stzscal( PTYPE, next+1, nt1, ATTR, NB, UPLO, M, N+mnk2,
                              mnk1, ALPHA, A, LDA );
         right = ATL_Stzscal( PTYPE, next+2, nt2, ATTR, NB, UPLO, 0, N,
                              mnk2, ALPHA, Mvptm( A, (M + mnk1) + mnk1 * LDA,
                              size ), LDA );
      }
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else if( split == AtlasTzSplitKcol )
   {
      if( UPLO == AtlasLower )
      {
         left  = ATL_Stzscal( PTYPE, next+1, nt1, ATTR, NB, UPLO, M+mnk2, N,
                              mnk1, ALPHA, A, LDA );
         right = ATL_Stzscal( PTYPE, next+2, nt2, ATTR, NB, UPLO, M,      0,
                              mnk2, ALPHA, Mvptm( A, mnk1 + (N + mnk1) * LDA,
                              size ), LDA );
      }
      else
      {
         left  = ATL_Stzscal( PTYPE, next+1, nt1, ATTR, NB, UPLO, M,      0,
                              mnk1, ALPHA, A, LDA );
         right = ATL_Stzscal( PTYPE, next+2, nt2, ATTR, NB, UPLO, M+mnk1, N,
                              mnk2, ALPHA, Mvptm( A, mnk1 * LDA, size ), LDA );
      }
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else if( split == AtlasTzSplitNcol )
   {
      if( UPLO == AtlasLower )
      {
         left  = ATL_Sgescal( PTYPE, next+1, nt1, ATTR, NB, M+K, mnk1, ALPHA,
                              A, LDA );
         right = ATL_Stzscal( PTYPE, next+2, nt2, ATTR, NB, UPLO, M, mnk2, K,
                              ALPHA, Mvptm( A, mnk1*LDA, size ),
                              LDA );
      }
      else
      {
         left  = ATL_Stzscal( PTYPE, next+1, nt1, ATTR, NB, UPLO, M, mnk1, K,
                              ALPHA, A, LDA );
         right = ATL_Sgescal( PTYPE, next+2, nt2, ATTR, NB, M+K, mnk2, ALPHA,
                              Mvptm( A, (K + mnk1) * LDA, size ),
                              LDA );
      }
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else
   {
      a_scal = (PT_TZSCAL_ARGS_T *)malloc( sizeof( PT_TZSCAL_ARGS_T ) );
      ATL_assert( a_scal != NULL );
      a_scal->up = UPLO;   a_scal->m = M; a_scal->n = N;    a_scal->k = K;
      a_scal->al = ALPHA;  a_scal->a = A; a_scal->la = LDA;
      tree = ATL_init_node( NODE, NULL, NULL, NULL, NULL, PTYPE->fun,
                            (void *)(a_scal) );
   }

   return( tree );
/*
 * End of ATL_Stzscal
 */
}
