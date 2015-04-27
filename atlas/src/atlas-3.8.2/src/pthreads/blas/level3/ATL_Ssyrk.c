
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
#include "atlas_ptlvl3.h"
#include "atlas_ptlevel3.h"

PT_TREE_T ATL_Ssyrk
(
   const PT_LVL3_TYPE_T       * PLVL3,
   const unsigned int         NODE,
   const unsigned int         NTHREADS,
   pthread_attr_t             * ATTR,
   const int                  NB,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const enum ATLAS_TRANS     TGEMM,
   const int                  M,
   const int                  N,
   const int                  K,
   const int                  L,
   const void                 * ALPHA,
   const void                 * A,
   const int                  LDA,
   const void                 * BETA,
   void                       * C,
   const int                  LDC
)
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  left, right, tree = NULL;
   PT_SYRK_ARGS_T             * a_syrk  = NULL;
   size_t                     size;
   DIM_TZSPLIT_T              split = AtlasTzNoSplit;
   unsigned int               next = ( NODE << 1 ), nt1=0, nt2=0;
   int                        incA, incC, mnk1=0, mnk2=0;
/* ..
 * .. Executable Statements ..
 *
 */
   if( NTHREADS < 1 ) return( tree );

   size  = PLVL3->size;
   split = ATL_tzsplit( UPLO, NTHREADS, M, N, K, NB, &nt1, &nt2,
                        &mnk1, &mnk2 );

   if(      split == AtlasTzSplitMrow )
   {
      if( UPLO == AtlasLower )
      {
         if( TRANS == AtlasNoTrans )
         {
            left  = ATL_Ssyrk( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS,
                               TGEMM, mnk1, N, K, L, ALPHA, A, LDA, BETA,
                               C, LDC );
            right = ATL_Sgemm( PLVL3, next+2, nt2, ATTR, NB, AtlasNoTrans,
                               TGEMM, mnk2, N+K, L, ALPHA,
                               Mvpt3( A, K + mnk1, size ), LDA,
                               Mvpt3( A,   - N,    size ), LDA, BETA,
                               Mvpt3( C, K + mnk1, size ), LDC );
         }
         else
         {
            left  = ATL_Ssyrk( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS,
                               TGEMM, mnk1, N, K, L, ALPHA, A, LDA, BETA,
                               C, LDC );
            right = ATL_Sgemm( PLVL3, next+2, nt2, ATTR, NB, TGEMM,
                               AtlasNoTrans, mnk2, N+K, L, ALPHA,
                               Mvpt3( A, ( N + K + mnk1 ) * LDA,
                               size ), LDA, A, LDA, BETA,
                               Mvpt3( C, K + mnk1, size ), LDC );
         }
      }
      else
      {
         if( TRANS == AtlasNoTrans )
         {
            left  = ATL_Sgemm( PLVL3, next+1, nt1, ATTR, NB, AtlasNoTrans,
                               TGEMM, mnk1, N+K, L, ALPHA, A, LDA,
                               Mvpt3( A, M, size ), LDA, BETA,
                               C, LDC );
            right = ATL_Ssyrk( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS,
                               TGEMM, mnk2, N, K, L,  ALPHA,
                               Mvpt3( A, mnk1, size ), LDA, BETA,
                               Mvpt3( C, mnk1, size ), LDC );
         }
         else
         {
            left  = ATL_Sgemm( PLVL3, next+1, nt1, ATTR, NB, TGEMM,
                               AtlasNoTrans, mnk1, N+K, L, ALPHA,
                               Mvpt3( A, - M * LDA, size ), LDA,
                               A, LDA, BETA, C, LDC );
            right = ATL_Ssyrk( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS,
                               TGEMM, mnk2, N, K, L,  ALPHA, A, LDA, BETA,
                               Mvpt3( C, mnk1, size ), LDC );
         }
      }

      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else if( split == AtlasTzSplitKrow )
   {
      if( UPLO == AtlasLower )
      {
         left  = ATL_Ssyrk( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS, TGEMM,
                            0, N, mnk1, L, ALPHA, A, LDA, BETA, C, LDC );
         incA  = ( TRANS == AtlasNoTrans ? mnk1 : 0 );
         incC  = mnk1;
         right = ATL_Ssyrk( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS, TGEMM,
                            M, N + mnk1, mnk2, L, ALPHA, Mvpt3( A, incA, size ),
                            LDA,  BETA, Mvpt3( C, incC, size ), LDC );
      }
      else
      {
         left  = ATL_Ssyrk( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS, TGEMM, M,
                            N + mnk2, mnk1, L, ALPHA, A, LDA, BETA, C, LDC );
         incA  = ( TRANS == AtlasNoTrans ? mnk1 : mnk1 * LDA );
         incC  = mnk1 * ( LDC + 1 );
         right = ATL_Ssyrk( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS, TGEMM, 0,
                            N, mnk2, L, ALPHA, Mvpt3( A, incA, size ),
                            LDA, BETA, Mvpt3( C, incC, size ), LDC );
      }
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else if( split == AtlasTzSplitKcol )
   {
      if( UPLO == AtlasLower )
      {
         left  = ATL_Ssyrk( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS, TGEMM,
                            M + mnk2, N, mnk1, L, ALPHA, A, LDA, BETA,
                            C, LDC );
         incA  = ( TRANS == AtlasNoTrans ? mnk1 : ( N + mnk1 ) * LDA );
         incC  = mnk1 + ( N + mnk1 ) * LDC;
         right = ATL_Ssyrk( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS, TGEMM,
                            M, 0, mnk2, L, ALPHA, Mvpt3( A, incA,
                            size ), LDA, BETA, Mvpt3( C, incC,
                            size ), LDC );
      }
      else
      {
         left  = ATL_Ssyrk( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS, TGEMM,
                            M, 0, mnk1, L, ALPHA, A, LDA, BETA, C, LDC );
         incA  = ( TRANS == AtlasNoTrans ? 0 : mnk1 * LDA );
         incC  = mnk1 * LDC;
         right = ATL_Ssyrk( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS, TGEMM,
                            M + mnk1, N, mnk2, L, ALPHA, Mvpt3( A,
                            incA, size ), LDA, BETA, Mvpt3( C, incC,
                            size ), LDC );
      }
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else if( split == AtlasTzSplitNcol )
   {
      if( UPLO == AtlasLower )
      {
         if( TRANS == AtlasNoTrans )
         {
            incA  = -N;
            left  = ATL_Sgemm( PLVL3, next+1, nt1, ATTR, NB, AtlasNoTrans,
                               TGEMM, M + K, mnk1, L, ALPHA, A, LDA,
                               Mvpt3( A, incA, size ), LDA, BETA,
                               C, LDC );
            incC  = mnk1 * LDC;
            right = ATL_Ssyrk( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS,
                               TGEMM, M, mnk2, K, L, ALPHA, A, LDA, BETA,
                               Mvpt3( C, incC, size ), LDC );
         }
         else
         {
            incA  = ( N + K ) * LDA;
            left  = ATL_Sgemm( PLVL3, next+1, nt1, ATTR, NB, TGEMM,
                               AtlasNoTrans, M + K, mnk1, L, ALPHA,
                               Mvpt3( A, incA, size ), LDA,
                               A, LDA, BETA, C, LDC );
            incA  = mnk1 * LDA;
            incC  = mnk1 * LDC;
            right = ATL_Ssyrk( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS,
                               TGEMM, M, mnk2, K, L, ALPHA,
                               Mvpt3( A, incA, size ), LDA, BETA,
                               Mvpt3( C, incC, size ), LDC );
         }
      }
      else
      {
         left  = ATL_Ssyrk( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS, TGEMM,
                            M, mnk1, K, L, ALPHA, A, LDA, BETA, C, LDC );
         incC  = ( K + mnk1 ) * LDC;

         if( TRANS == AtlasNoTrans )
         {
            incA  = M + K + mnk1;
            right = ATL_Sgemm( PLVL3, next+2, nt2, ATTR, NB, AtlasNoTrans,
                               TGEMM, M + K, mnk2, L, ALPHA, A, LDA,
                               Mvpt3( A, incA, size ), LDA, BETA,
                               Mvpt3( C, incC, size ), LDC );
         }
         else
         {
            incA  = ( K + mnk1 ) * LDA;
            right = ATL_Sgemm( PLVL3, next+2, nt2, ATTR, NB, TGEMM,
                               AtlasNoTrans, M + K, mnk2, L, ALPHA,
                               Mvpt3( A, -M * LDA, size ), LDA,
                               Mvpt3( A,     incA, size ), LDA, BETA,
                               Mvpt3( C,     incC, size ), LDC );
         }
      }
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else
   {
      a_syrk = (PT_SYRK_ARGS_T *)malloc( sizeof( PT_SYRK_ARGS_T ) );
      ATL_assert( a_syrk != NULL );
      a_syrk->up = UPLO;  a_syrk->tr = TRANS; a_syrk->m  = M;
      a_syrk->n  = N;     a_syrk->k  = K;     a_syrk->l  = L;
      a_syrk->al = ALPHA; a_syrk->a  = A;     a_syrk->la = LDA;
      a_syrk->be = BETA;  a_syrk->c  = C;     a_syrk->lc = LDC;

      tree = ATL_init_node( NODE, NULL, NULL, NULL, NULL,
                            ( TGEMM == AtlasTrans ? PLVL3->syrk0 :
                            PLVL3->herk0 ), (void *)(a_syrk) );
   }

   return( tree );
/*
 * End of ATL_Ssyrk
 */
}
