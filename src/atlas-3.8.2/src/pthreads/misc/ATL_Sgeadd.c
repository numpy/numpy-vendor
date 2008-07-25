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

PT_TREE_T ATL_Sgeadd
(
   const PT_MISC_TYPE_T       * PTYPE,
   const unsigned int         NODE,
   const unsigned int         THREADS,
   pthread_attr_t             * ATTR,
   const int                  NB,
   const int                  M,
   const int                  N,
   const void                 * ALPHA,
   const void                 * A,
   const int                  LDA,
   const void                 * BETA,
   void                       * C,
   const int                  LDC
)
{
/*
 * Purpose
 * =======
 *
 * ATL_Sgeadd  recursively  builds  the binary tasks tree to be threaded
 * for the matrix add operation.
 *
 * Arguments
 * =========
 *
 * PTYPE   (input)                       const PT_MISC_TYPE_T *
 *         On entry, PTYPE  points  to the data structure containing the
 *         type information.
 *
 * NODE    (input)                       const unsigned int
 *         On entry, NODE specifies the current node number.
 *
 * THREADS (input)                       const unsigned int
 *         On entry, THREADS  specifies the number of threads to be used
 *         for the current operation.
 *
 * ATTR    (input)                       pthread_attr_t *
 *         On entry, ATTR  specifies  the  thread attribute object to be
 *         used for the node functions to be threaded.
 *
 * NB      (input)                       const int
 *         On entry, NB  specifies  the  blocksize  to  be  used for the
 *         problem size partitioning.
 *
 * M       (input)                       const int
 * N       (input)                       const int
 * ALPHA   (input)                       const void *
 * A       (input)                       const void *
 * LDA     (input)                       const int
 * BETA    (input)                       const void *
 * C       (input/output)                void *
 * LDC     (input)                       const int
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  left, right, tree = NULL;
   PT_GEADD_ARGS_T            * a_add = NULL;
   double                     rm, rn;
   DIM_1DSPLIT_T              splitM, splitN;
   size_t                     size;
   int                        split, m1, m2, n1, n2;
   unsigned int               next = (NODE << 1), ntm1, ntm2, ntn1, ntn2;
/* ..
 * .. Executable Statements ..
 *
 */
   if( THREADS < 1 ) return( tree );

   size   = PTYPE->size;
   splitM = ATL_1dsplit( THREADS, M, NB, &ntm1, &ntm2, &m1, &m2, &rm );
   splitN = ATL_1dsplit( THREADS, N, NB, &ntn1, &ntn2, &n1, &n2, &rn );

   if( ( splitM == Atlas1dSplit ) && ( splitN == Atlas1dSplit ) )
   { split = ( rm < rn ? SPLIT_M : SPLIT_N ); }
   else if( splitM == Atlas1dSplit ) { split = SPLIT_M; }
   else if( splitN == Atlas1dSplit ) { split = SPLIT_N; }
   else                              { split = NOSPLIT; }

   if(      split == SPLIT_N )
   {
      left  = ATL_Sgeadd( PTYPE, next+1, ntn1, ATTR, NB, M, n1, ALPHA,
                          A, LDA, BETA, C, LDC );
      right = ATL_Sgeadd( PTYPE, next+2, ntn2, ATTR, NB, M, n2, ALPHA,
                          Mvptm( A, n1*LDA, size ), LDA, BETA,
                          Mvptm( C, n1*LDC, size ), LDC );
      tree  = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else if( split == SPLIT_M )
   {
      left  = ATL_Sgeadd( PTYPE, next+1, ntm1, ATTR, NB, m1, N, ALPHA,
                          A, LDA, BETA, C, LDC );
      right = ATL_Sgeadd( PTYPE, next+2, ntm2, ATTR, NB, m2, N, ALPHA,
                          Mvptm( A, m1, size ), LDA, BETA,
                          Mvptm( C, m1, size ), LDC );
      tree  = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
   }
   else
   {
      a_add = (PT_GEADD_ARGS_T *)malloc( sizeof( PT_GEADD_ARGS_T ) );
      ATL_assert( a_add != NULL );
      a_add->m  = M;     a_add->n = N;
      a_add->al = ALPHA; a_add->a = A; a_add->la = LDA;
      a_add->be = BETA;  a_add->c = C; a_add->lc = LDC;
      tree = ATL_init_node( NODE, NULL, NULL, NULL, NULL, PTYPE->fun,
                            (void *)(a_add) );
   }

   return( tree );
/*
 * End of ATL_Sgeadd
 */
}
