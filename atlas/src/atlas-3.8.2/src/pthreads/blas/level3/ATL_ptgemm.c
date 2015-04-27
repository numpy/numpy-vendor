
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

PT_FUN_ARG_T Mjoin( PATL, ptgemm0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = (PT_TREE_T)(ARGS);
   PT_GEMM_ARGS_T             * arg;
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_wait_tree( root );
   arg = (PT_GEMM_ARGS_T *)(root->arg);

#define  _PT_TRANSA    (arg->ta)
#define  _PT_TRANSB    (arg->tb)
#define  _PT_M         (arg->m)
#define  _PT_N         (arg->n)
#define  _PT_K         (arg->k)

#ifdef TREAL
#define  _PT_ALPHA     (SCALAR)(*((TYPE *)(arg->al)))
#define  _PT_BETA      (SCALAR)(*((TYPE *)(arg->be)))
#else
#define  _PT_ALPHA     (SCALAR)(arg->al)
#define  _PT_BETA      (SCALAR)(arg->be)
#endif

#define  _PT_A         ((TYPE *)(arg->a))
#define  _PT_LDA       (arg->la)
#define  _PT_B         ((TYPE *)(arg->b))
#define  _PT_LDB       (arg->lb)
#define  _PT_C         ((TYPE *)(arg->c))
#define  _PT_LDC       (arg->lc)

   Mjoin( PATL, gemm )( _PT_TRANSA, _PT_TRANSB, _PT_M, _PT_N, _PT_K,
                        _PT_ALPHA, _PT_A, _PT_LDA, _PT_B, _PT_LDB,
                        _PT_BETA,  _PT_C, _PT_LDC );
   ATL_signal_tree( root );

   return( NULL );
/*
 * End of Mjoin( PATL, ptgemm0 )
 */
}

void Mjoin( PATL, ptgemm )
(
   const enum ATLAS_TRANS     TRANSA,
   const enum ATLAS_TRANS     TRANSB,
   const int                  M,
   const int                  N,
   const int                  K,
   const SCALAR               ALPHA,
   const TYPE                 * A,
   const int                  LDA,
   const TYPE                 * B,
   const int                  LDB,
   const SCALAR               BETA,
   TYPE                       * C,
   const int                  LDC
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptgemm )  performs one of the matrix-matrix operations
 *
 *    C := alpha * op( A ) * op( B ) + beta * C,
 *
 * where op( X ) is one of
 *
 *    op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ).
 *
 * Alpha and beta are scalars, and A, B and C are matrices, with op( A )
 * an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
 *
 * For a  more  detailed description of  the arguments of this function,
 * see the reference implementation in the  ATLAS/src/blas/reference di-
 * rectory.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = NULL;
   pthread_attr_t             attr;
#ifdef TREAL
   TYPE                       alpha0 = ALPHA, beta0 = BETA;
#endif
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( M == 0 ) || ( N == 0 ) ||
       ( ( SCALAR_IS_ZERO( ALPHA ) || ( K == 0 ) ) &&
         SCALAR_IS_ONE( BETA ) ) ) return;

   if( SCALAR_IS_ZERO( ALPHA ) )
   {
      if( SCALAR_IS_ZERO( BETA ) ) { Mjoin( PATL, ptgezero )( M, N, C, LDC ); }
      else                   { Mjoin( PATL, ptgescal )( M, N, BETA, C, LDC ); }
      return;
   }

   ATL_thread_init( &attr );
#ifdef TREAL
   root = Mjoin( PATL, ptgemm_nt )( ATL_NTHREADS, &attr, TRANSA, TRANSB,
                                    M, N, K, (void *)(&alpha0), (void *)(A),
                                    LDA, (void *)(B), LDB, (void *)(&beta0),
                                    (void *)(C), LDC );
#else
   root = Mjoin( PATL, ptgemm_nt )( ATL_NTHREADS, &attr, TRANSA, TRANSB,
                                    M, N, K, (void *)(ALPHA), (void *)(A),
                                    LDA, (void *)(B), LDB, (void *)(BETA),
                                    (void *)(C), LDC );
#endif
   ATL_join_tree  ( root );
   ATL_free_tree  ( root );
   ATL_thread_exit( &attr );
/*
 * End of Mjoin( PATL, ptgemm )
 */
}

PT_TREE_T Mjoin( PATL, ptgemm_nt )
(
   const unsigned int         NTHREADS,
   pthread_attr_t             * ATTR,
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
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptgemm_nt )  performs one of the matrix-matrix operations
 *
 *    C := alpha * op( A ) * op( B ) + beta * C,
 *
 * where op( X ) is one of
 *
 *    op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ).
 *
 * Alpha and beta are scalars, and A, B and C are matrices, with op( A )
 * an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
 *
 * For a  more  detailed description of  the arguments of this function,
 * see the reference implementation in the  ATLAS/src/blas/reference di-
 * rectory.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = NULL;
   PT_LVL3_TYPE_T             type;
   double                     tblks;
#ifdef TREAL
   TYPE                       alpha, beta;
#else
   TYPE                       alpha[2], beta[2];
#endif
   unsigned int               nthreads;
   int                        nb, nbm1, NOTHR;
/* ..
 * .. Executable Statements ..
 *
 */
   nbm1  = ( nb = Mjoin( PATL, GetNB )() ) - 1;
   tblks = (double)( (M+nbm1) / nb ) * (double)( (N+nbm1) / nb ) *
           (double)( (K+nbm1) / nb );
/*
 * When two dimensions or more are less than NB, do not thread.
 * Changed by RCW to allow threading if tblks is very large.  Need to replace
 * this with a probe.
 */
   NOTHR = (tblks <= (double)(ATL_XOVER_L3_DEFAULT)) || (NTHREADS <= 1);
   NOTHR = NOTHR & ( ((M<=nb && N<=nb) || (M<=nb && K<=nb) || (N<=nb && K<=nb))
                     && tblks < 256 );
   if (NOTHR)
   {
#ifdef TREAL
      alpha    = *((TYPE *)(ALPHA));   beta     = *((TYPE *)(BETA));
#else
      alpha[0] = ((TYPE *)(ALPHA))[0]; alpha[1] = ((TYPE *)(ALPHA))[1];
      beta [0] = ((TYPE *)(BETA ))[0]; beta [1] = ((TYPE *)(BETA ))[1];
#endif
      Mjoin( PATL, gemm )( TRANSA, TRANSB, M, N, K, alpha, (TYPE *)(A), LDA,
                           (TYPE *)(B), LDB, beta, (TYPE *)(C), LDC );
      return( root );
   }

   if( tblks >= (double)(NTHREADS) ) { nthreads = NTHREADS; }
   else    { nthreads = (unsigned int)floor( tblks + 0.5 ); }

   Mjoin( PATL, ptl3settype )( &type );
   root = ATL_Sgemm( &type, 0, nthreads, ATTR, nb, TRANSA, TRANSB, M, N, K,
                     ALPHA, A, LDA, B, LDB, BETA, C, LDC );
   ATL_thread_tree( root, ATTR );

   return( root );
/*
 * End of Mjoin( PATL, ptgemm_nt )
 */
}
