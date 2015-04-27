
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

PT_TREE_T ATL_Ssyr2k
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
   const void                 * ALPHC,
   const void                 * A,
   const int                  LDA,
   const void                 * B,
   const int                  LDB,
   const void                 * BETA,
   void                       * C,
   const int                  LDC
)
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  left, right, tree = NULL, root;
   PT_SYR2K_ARGS_T            * a_syr2k  = NULL;
#if 1
   double                     rn;
   DIM_1DSPLIT_T              splitN;
   size_t                     size;
   unsigned int               next = ( NODE << 1 ), nt1, nt2;
   int                        failed, incA, incB, incC, mnk1, mnk2;
#else
   size_t                     size;
   DIM_TZSPLIT_T              split = AtlasTzNoSplit;
   unsigned int               next = ( NODE << 1 ), nt1, nt2;
   int                        failed, incA, incB, incC, mnk1, mnk2;
#endif
/* ..
 * .. Executable Statements ..
 *
 */
   if( NTHREADS < 1 ) return( tree );

   size  = PLVL3->size;

   if( M < 0 || N < 0 ) return( tree ); /* temporary */

#if 1
   if( TGEMM == AtlasTrans )
   {
      failed = PLVL3->ptsyr2k0( NTHREADS, ATTR, UPLO, TRANS, TGEMM, K, L,
                                ALPHA, A, LDA, B, LDB, BETA, C, LDC );
   }
   else
   {
      failed = PLVL3->pther2k0( NTHREADS, ATTR, UPLO, TRANS, TGEMM, K, L,
                                ALPHA, A, LDA, B, LDB, BETA, C, LDC );
   }

   if( failed )
   {
      splitN = ATL_1dsplit( NTHREADS, K, NB, &nt1, &nt2, &mnk1, &mnk2, &rn );

      if( splitN == Atlas1dSplit )
      {
         left  = ATL_Ssyr2k( PLVL3, next+1, NTHREADS, ATTR, NB, UPLO, TRANS,
                             TGEMM, 0, 0, mnk1, L, ALPHA, ALPHC, A, LDA,
                             B, LDB, BETA, C, LDC );

         if( UPLO == AtlasLower )
         {
            if( TRANS == AtlasNoTrans )
            {
               incA  = mnk1; incB  = mnk1; incC  = mnk1;

               root = PLVL3->ptgemm( NTHREADS, ATTR, AtlasNoTrans, TGEMM, mnk2,
                                     mnk1, L, ALPHA, Mvpt3( A, incA, size ),
                                     LDA, B, LDB, BETA, Mvpt3( C, incC, size ),
                                     LDC );
               ATL_join_tree  ( root );
               ATL_free_tree  ( root );

               root = PLVL3->ptgemm( NTHREADS, ATTR, AtlasNoTrans, TGEMM, mnk2,
                                     mnk1, L, ALPHC, Mvpt3( B, incB, size ),
                                     LDB, A, LDA, PLVL3->one, Mvpt3( C, incC,
                                     size ), LDC );
               ATL_join_tree  ( root );
               ATL_free_tree  ( root );

               incC += mnk1 * LDC;
            }
            else
            {
               incA  = mnk1 * LDA; incB  = mnk1 * LDB; incC  = mnk1;

               root = PLVL3->ptgemm( NTHREADS, ATTR, TGEMM, AtlasNoTrans, mnk2,
                                     mnk1, L, ALPHA, Mvpt3( A, incA, size ),
                                     LDA, B, LDB, BETA, Mvpt3( C, incC, size ),
                                     LDC );
               ATL_join_tree  ( root );
               ATL_free_tree  ( root );

               root = PLVL3->ptgemm( NTHREADS, ATTR, TGEMM, AtlasNoTrans, mnk2,
                                     mnk1, L, ALPHC, Mvpt3( B, incB, size ),
                                     LDB, A, LDA, PLVL3->one, Mvpt3( C, incC,
                                     size ), LDC );
               ATL_join_tree  ( root );
               ATL_free_tree  ( root );

               incC += mnk1 * LDC;
            }

         }
         else
         {
            if( TRANS == AtlasNoTrans )
            {
               incA  = mnk1; incB  = mnk1; incC  = mnk1 * LDC;

               root = PLVL3->ptgemm( NTHREADS, ATTR, AtlasNoTrans, TGEMM, mnk1,
                                     mnk2, L, ALPHA, A, LDA, Mvpt3( B, incB,
                                     size ), LDB, BETA, Mvpt3( C, incC, size ),
                                     LDC );
               ATL_join_tree  ( root );
               ATL_free_tree  ( root );

               root = PLVL3->ptgemm( NTHREADS, ATTR, AtlasNoTrans, TGEMM, mnk1,
                                     mnk2, L, ALPHC, B, LDB, Mvpt3( A, incA,
                                     size ), LDA, PLVL3->one, Mvpt3( C, incC,
                                     size ), LDC );
               ATL_join_tree  ( root );
               ATL_free_tree  ( root );

               incC += mnk1;
            }
            else
            {
               incA  = mnk1 * LDA; incB  = mnk1 * LDB; incC  = mnk1 * LDC;

               root = PLVL3->ptgemm( NTHREADS, ATTR, TGEMM, AtlasNoTrans, mnk1,
                                     mnk2, L, ALPHA, A, LDA, Mvpt3( B, incB,
                                     size ), LDB, BETA, Mvpt3( C, incC, size ),
                                     LDC );
               ATL_join_tree  ( root );
               ATL_free_tree  ( root );

               root = PLVL3->ptgemm( NTHREADS, ATTR, TGEMM, AtlasNoTrans, mnk1,
                                     mnk2, L, ALPHC, B, LDB, Mvpt3( A, incA,
                                     size ), LDA, PLVL3->one, Mvpt3( C, incC,
                                     size ), LDC );
               ATL_join_tree  ( root );
               ATL_free_tree  ( root );

               incC += mnk1;
            }
         }

         right = ATL_Ssyr2k( PLVL3, next+2, NTHREADS, ATTR, NB, UPLO, TRANS,
                             TGEMM, 0, 0, mnk2, L, ALPHA, ALPHC,
                             Mvpt3( A, incA, size ), LDA,
                             Mvpt3( B, incB, size ), LDB, BETA,
                             Mvpt3( C, incC, size ), LDC );
      }
      else
      {
         a_syr2k = (PT_SYR2K_ARGS_T *)malloc( sizeof( PT_SYR2K_ARGS_T ) );
         ATL_assert( a_syr2k != NULL );
         a_syr2k->up = UPLO;  a_syr2k->tr = TRANS; a_syr2k->m  = 0;
         a_syr2k->n  = 0;     a_syr2k->k  = K;     a_syr2k->l  = L;
         a_syr2k->al = ALPHA; a_syr2k->ac = ALPHA; a_syr2k->a  = A;
         a_syr2k->la = LDA;   a_syr2k->b  = B;     a_syr2k->lb = LDB;
         a_syr2k->be = BETA;  a_syr2k->c  = C;     a_syr2k->lc = LDC;

         tree = ATL_init_node( NODE, NULL, NULL, NULL, NULL,
                               ( TGEMM == AtlasTrans ? PLVL3->syr2k0 :
                               PLVL3->her2k0 ), (void *)(a_syr2k) );
      }
   }

   return( tree );
#else
   split = ATL_tzsplit( UPLO, NTHREADS, M, N, K, NB, &nt1, &nt2,
                        &mnk1, &mnk2 );

   if( split == AtlasSplitKrow )
   {
      if( UPLO == AtlasLower )
      {
         left  = ATL_Ssyr2k( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS,
                             TGEMM, 0, N, mnk1, L, ALPHA, ALPHC, A, LDA,
                             B, LDB, BETA, C, LDC );
         incA  = ( TRANS == AtlasNoTrans ? mnk1 : 0 );
         incB  = ( TRANS == AtlasNoTrans ? mnk1 : 0 );
         incC  = mnk1;
         right = ATL_Ssyr2k( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS,
                             TGEMM, M, N + mnk1, mnk2, L, ALPHA, ALPHC,
                             Mvpt3( A, incA, size ), LDA,
                             Mvpt3( B, incB, size ), LDB, BETA,
                             Mvpt3( C, incC, size ), LDC );
      }
      else
      {
         left  = ATL_Ssyr2k( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS,
                             TGEMM, M, N + mnk2, mnk1, L, ALPHA, ALPHC,
                             A, LDA, B, LDB, BETA, C, LDC );
         incA  = ( TRANS == AtlasNoTrans ? mnk1 : mnk1 * LDA );
         incB  = ( TRANS == AtlasNoTrans ? mnk1 : mnk1 * LDB );
         incC  = mnk1 * ( LDC + 1 );
         right = ATL_Ssyr2k( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS,
                             TGEMM, 0, N, mnk2, L, ALPHA, ALPHC,
                             Mvpt3( A, incA, size ), LDA,
                             Mvpt3( B, incB, size ), LDB, BETA,
                             Mvpt3( C, incC, size ), LDC );
      }
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else if( split == AtlasSplitKcol )
   {
      if( UPLO == AtlasLower )
      {
         left  = ATL_Ssyr2k( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS,
                             TGEMM, M + mnk2, N, mnk1, L, ALPHA, ALPHC,
                             A, LDA, B, LDB, BETA, C, LDC );
         incA  = ( TRANS == AtlasNoTrans ? mnk1 : ( N + mnk1 ) * LDA );
         incB  = ( TRANS == AtlasNoTrans ? mnk1 : ( N + mnk1 ) * LDB );
         incC  = mnk1 + ( N + mnk1 ) * LDC;
         right = ATL_Ssyr2k( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS,
                             TGEMM, M, 0, mnk2, L, ALPHA, ALPHC,
                             Mvpt3( A, incA, size ), LDA,
                             Mvpt3( B, incB, size ), LDB, BETA,
                             Mvpt3( C, incC, size ), LDC );
      }
      else
      {
         left  = ATL_Ssyr2k( PLVL3, next+1, nt1, ATTR, NB, UPLO, TRANS,
                             TGEMM, M, 0, mnk1, L, ALPHA, ALPHC, A, LDA,
                             B, LDB, BETA, C, LDC );
         incA  = ( TRANS == AtlasNoTrans ? 0 : mnk1 * LDA );
         incB  = ( TRANS == AtlasNoTrans ? 0 : mnk1 * LDB );
         incC  = mnk1 * LDC;
         right = ATL_Ssyr2k( PLVL3, next+2, nt2, ATTR, NB, UPLO, TRANS,
                             TGEMM, M + mnk1, N, mnk2, L, ALPHA, ALPHC,
                             Mvpt3( A, incA, size ), LDA,
                             Mvpt3( B, incB, size ), LDB, BETA,
                             Mvpt3( C, incC, size ), LDC );
      }
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else
   {
      a_syr2k = (PT_SYR2K_ARGS_T *)malloc( sizeof( PT_SYR2K_ARGS_T ) );
      ATL_assert( a_syr2k != NULL );
      a_syr2k->up = UPLO;  a_syr2k->tr = TRANS; a_syr2k->m  = M;
      a_syr2k->n  = N;     a_syr2k->k  = K;     a_syr2k->l  = L;
      a_syr2k->al = ALPHA; a_syr2k->ac = ALPHA; a_syr2k->a  = A;
      a_syr2k->la = LDA;   a_syr2k->b  = B;     a_syr2k->lb = LDB;
      a_syr2k->be = BETA;  a_syr2k->c  = C;     a_syr2k->lc = LDC;

      tree = ATL_init_node( NODE, NULL, NULL, NULL, NULL,
                            ( TGEMM == AtlasTrans ? PLVL3->syr2k0 :
                            PLVL3->her2k0 ), (void *)(a_syr2k) );
   }

   return( tree );
#endif
/*
 * End of ATL_Ssyr2k
 */
}
