
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

PT_TREE_T ATL_Strmm
(
   const PT_LVL3_TYPE_T       * PLVL3,
   const unsigned int         NODE,
   const unsigned int         NTHREADS,
   pthread_attr_t             * ATTR,
   const int                  NB,
   const enum ATLAS_SIDE      SIDE,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const enum ATLAS_DIAG      DIAG,
   const int                  M,
   const int                  N,
   const void                 * ALPHA,
   const void                 * A,
   const int                  LDA,
   void                       * B,
   const int                  LDB
)
{
/*
 * .. Local Variables ..
 */
#if 0
   PT_TREE_T                  left, right, root = NULL;
   PT_TRMM_ARGS_T             * a_mm = NULL;
   size_t                     size;
   unsigned int               next = ( NODE << 1 ), nt1, nt2;
   int                        incA, incB, m1, m2, n1, n2;
#else
   PT_TREE_T                  left, right, tree = NULL;
   PT_TRMM_ARGS_T             * a_mm = NULL;
   double                     rm, rn;
   DIM_1DSPLIT_T              splitM, splitN;
   size_t                     size;
   unsigned int               next = ( NODE << 1 ), ntm1, ntm2, ntn1, ntn2;
   int                        incB, m1, m2, n1, n2;
#endif
/* ..
 * .. Executable Statements ..
 *
 */
   if( NTHREADS < 1 ) return( tree );

   size = PLVL3->size;

   if( SIDE == AtlasLeft )
   {
#if 0
      if( M > N )
      {
         ATL_gesplit( NTHREADS, M, NB, &nt1, &nt2, &m1, &m2 );
         incA  = m1 + m1 * LDA;
         incB  = m1;
         root = PLVL3->pttrmm( NTHREADS, ATTR, SIDE, UPLO, TRANS, DIAG,
                               m2, N, ALPHA, Mvpt3( A, incA, size ), LDA,
                               Mvpt3( B, incB, size ), LDB );
         ATL_join_tree  ( root );
         ATL_free_tree  ( root );
         root = PLVL3->ptgemm( NTHREADS, ATTR, AtlasNoTrans, AtlasNoTrans,
                               m2, N, m1, ALPHA,
                               Mvpt3( A, m1, size ), LDA,
                               B, LDB, PLVL3->one,
                               Mvpt3( B, m1, size ), LDB );
         ATL_join_tree  ( root );
         ATL_free_tree  ( root );
         root = PLVL3->pttrmm( NTHREADS, ATTR, SIDE, UPLO, TRANS, DIAG,
                               m1, N, ALPHA, A, LDA, B, LDB );
         ATL_join_tree  ( root );
         ATL_free_tree  ( root );

         return( NULL );
      }
      else
      {
         ATL_gesplit( NTHREADS, N, NB, &nt1, &nt2, &n1, &n2 );
         left  = ATL_Strmm( PLVL3, next+1, nt1, ATTR, NB, SIDE, UPLO, TRANS,
                            DIAG, M, n1, ALPHA, A, LDA, B, LDB );
         incB  = n1 * LDB;
         right = ATL_Strmm( PLVL3, next+2, nt2, ATTR, NB, SIDE, UPLO, TRANS,
                            DIAG, M, n2, ALPHA, A, LDA, Mvpt3( B, incB, size ),
                            LDB );
      }
#else
      splitN = ATL_1dsplit( NTHREADS, N, NB, &ntn1, &ntn2, &n1, &n2, &rn );

      if( splitN == Atlas1dSplit )
      {
         left  = ATL_Strmm( PLVL3, next+1, ntn1, ATTR, NB, SIDE, UPLO, TRANS,
                            DIAG, M, n1, ALPHA, A, LDA, B, LDB );
         incB  = n1 * LDB;
         right = ATL_Strmm( PLVL3, next+2, ntn2, ATTR, NB, SIDE, UPLO, TRANS,
                            DIAG, M, n2, ALPHA, A, LDA, Mvpt3( B, incB,
                            size ), LDB );
         tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
      }
      else
      {
         a_mm = (PT_TRMM_ARGS_T *)malloc( sizeof( PT_TRMM_ARGS_T ) );
         ATL_assert( a_mm != NULL );
         a_mm->si = SIDE; a_mm->up = UPLO; a_mm->tr = TRANS; a_mm->di = DIAG;
         a_mm->m  = M;    a_mm->n  = N;    a_mm->al = ALPHA;
         a_mm->a  = A;    a_mm->la = LDA;  a_mm->b  = B;     a_mm->lb = LDB;
         tree = ATL_init_node( NODE, NULL, NULL, NULL, NULL, PLVL3->trmm0,
                               (void *)(a_mm) );
      }
#endif
   }
   else
   {
      splitM = ATL_1dsplit( NTHREADS, M, NB, &ntm1, &ntm2, &m1, &m2, &rm );

      if( splitM == Atlas1dSplit )
      {
         left  = ATL_Strmm( PLVL3, next+1, ntm1, ATTR, NB, SIDE, UPLO, TRANS,
                            DIAG, m1, N, ALPHA, A, LDA, B, LDB );
         incB  = m1;
         right = ATL_Strmm( PLVL3, next+2, ntm2, ATTR, NB, SIDE, UPLO, TRANS,
                            DIAG, m2, N, ALPHA, A, LDA, Mvpt3( B, incB,
                            size ), LDB );
         tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
      }
      else
      {
         a_mm = (PT_TRMM_ARGS_T *)malloc( sizeof( PT_TRMM_ARGS_T ) );
         ATL_assert( a_mm != NULL );
         a_mm->si = SIDE; a_mm->up = UPLO; a_mm->tr = TRANS; a_mm->di = DIAG;
         a_mm->m  = M;    a_mm->n  = N;    a_mm->al = ALPHA;
         a_mm->a  = A;    a_mm->la = LDA;  a_mm->b  = B;     a_mm->lb = LDB;
         tree = ATL_init_node( NODE, NULL, NULL, NULL, NULL, PLVL3->trmm0,
                               (void *)(a_mm) );
      }
   }

   return( tree );
/*
 * End of ATL_Strmm
 */
}
