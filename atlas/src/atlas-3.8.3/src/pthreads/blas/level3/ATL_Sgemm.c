
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

PT_TREE_T ATL_Sgemm
(
   const PT_LVL3_TYPE_T       * PLVL3,
   const unsigned int         NODE,
   const unsigned int         NTHREADS,
   pthread_attr_t             * ATTR,
   const int                  NB,
   const enum ATLAS_TRANS     TRANSA,
   const enum ATLAS_TRANS     TRANSB,
   const int                  M,
   const int                  N,
   const int                  K,
   const void                 * ALPHA,
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
   PT_TREE_T                  left, right, tree = NULL;
   PT_GEMM_ARGS_T             * a_mm  = NULL;
   PT_GEADD_ARGS_T            * a_add = NULL;
   void                       * C0, * vp = NULL;
   double                     rk, rm, rn;
   size_t                     size;
   DIM_1DSPLIT_T              splitM, splitN, splitK;
   unsigned int               next = (NODE << 1), ntk1, ntk2, ntm1, ntm2,
                              ntn1, ntn2;
   int                        incA, incB, incC, k1, k2, m1, m2, n1, n2,
                              split;
/* ..
 * .. Executable Statements ..
 *
 */
   if( NTHREADS < 1 ) return( tree );

   size   = PLVL3->size;
   splitM = ATL_1dsplit( NTHREADS, M, NB, &ntm1, &ntm2, &m1, &m2, &rm );
   splitN = ATL_1dsplit( NTHREADS, N, NB, &ntn1, &ntn2, &n1, &n2, &rn );
   splitK = ATL_1dsplit( NTHREADS, K, NB, &ntk1, &ntk2, &k1, &k2, &rk );

   if( ( splitM == Atlas1dSplit ) && ( splitN == Atlas1dSplit ) )
   {
      split = ( rm < rn ? SPLIT_M : SPLIT_N );
   }
   else if( splitM == Atlas1dSplit )
   {
      split = SPLIT_M;
   }
   else if( splitN == Atlas1dSplit )
   {
      split = SPLIT_N;
   }
   else
   {
      split = NOSPLIT;
      if( ( splitK == Atlas1dSplit ) && ( ( M <= N ) || ( M <= NB ) ) &&
          ( N <= NB ) )
      {
         vp = (void *)malloc( ATL_Cachelen + ( M * N * size ) );
         if( vp != NULL ) split = SPLIT_K;
      }
   }

   if(      split == SPLIT_M )
   {
      left  = ATL_Sgemm( PLVL3, next+1, ntm1, ATTR, NB, TRANSA, TRANSB,
                         m1, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC );
      incA  = ( TRANSA == AtlasNoTrans ? m1 : m1 * LDA ); incC = m1;
      right = ATL_Sgemm( PLVL3, next+2, ntm2, ATTR, NB, TRANSA, TRANSB,
                         m2, N, K, ALPHA, Mvpt3( A, incA, size ), LDA,
                         B, LDB, BETA, Mvpt3( C, incC, size ), LDC );
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else if( split == SPLIT_N )
   {
      left  = ATL_Sgemm( PLVL3, next+1, ntn1, ATTR, NB, TRANSA, TRANSB,
                         M, n1, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC );
      incB  = ( TRANSB == AtlasNoTrans ? n1 * LDB : n1 ); incC  = n1 * LDC;
      right = ATL_Sgemm( PLVL3, next+2, ntn2, ATTR, NB, TRANSA, TRANSB,
                         M, n2, K, ALPHA, A, LDA, Mvpt3( B, incB, size ),
                         LDB, BETA, Mvpt3( C, incC, size ), LDC );
      tree = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else if( split == SPLIT_K )
   {
      left  = ATL_Sgemm( PLVL3, next+1, ntk1, ATTR, NB, TRANSA, TRANSB,
                         M, N, k1, ALPHA, A, LDA, B, LDB, BETA, C, LDC );
      C0    = ATL_AlignPtr( vp );
      incA  = ( TRANSA == AtlasNoTrans ? k1 * LDA : k1 );
      incB  = ( TRANSB == AtlasNoTrans ? k1 : k1 * LDB );
      right = ATL_Sgemm( PLVL3, next+2, ntk2, ATTR, NB, TRANSA, TRANSB,
                         M, N, k2, PLVL3->one, Mvpt3( A, incA, size ),
                         LDA, Mvpt3( B, incB, size ), LDB, PLVL3->zero,
                         C0, M );
      a_add = (PT_GEADD_ARGS_T *)malloc( sizeof( PT_GEADD_ARGS_T ) );
      ATL_assert( a_add != NULL );
      a_add->m = M;  a_add->n  = N;   a_add->al = ALPHA;
      a_add->a = C0; a_add->la = M;   a_add->be = PLVL3->one;
      a_add->c = C;  a_add->lc = LDC;
      tree = ATL_init_node( NODE, left, right, vp, NULL, PLVL3->geadd0,
                            (void *)(a_add) );
   }
   else
   {
      a_mm = (PT_GEMM_ARGS_T *)malloc( sizeof( PT_GEMM_ARGS_T ) );
      ATL_assert( a_mm != NULL );
      a_mm->ta = TRANSA; a_mm->tb = TRANSB;
      a_mm->m  = M;      a_mm->n  = N;   a_mm->k = K;
      a_mm->al = ALPHA;  a_mm->a  = A;   a_mm->la = LDA;
      a_mm->b  = B;      a_mm->lb = LDB; a_mm->be = BETA;
      a_mm->c  = C;      a_mm->lc = LDC;
      tree = ATL_init_node( NODE, NULL, NULL, NULL, NULL, PLVL3->gemm0,
                            (void *)(a_mm) );
   }

   return( tree );
/*
 * End of ATL_Sgemm
 */
}
