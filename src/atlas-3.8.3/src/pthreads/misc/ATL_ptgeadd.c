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

PT_FUN_ARG_T Mjoin( PATL, ptgeadd0 )( PT_FUN_ARG_T ARGS )
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptgeadd0 ) adds an m-by-n matrix A to the matrix B.
 *
 * This is a multi-threaded version of the algorithm.
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
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = (PT_TREE_T)(ARGS);
   PT_GEADD_ARGS_T            * arg;
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_wait_tree( root );
   arg = (PT_GEADD_ARGS_T *)(root->arg);

#define  _PT_M         (arg->m)
#define  _PT_N         (arg->n)

#ifdef TREAL
#define  _PT_ALPHA     (SCALAR)(*((TYPE *)(arg->al)))
#define  _PT_BETA      (SCALAR)(*((TYPE *)(arg->be)))
#else
#define  _PT_ALPHA     (SCALAR)(arg->al)
#define  _PT_BETA      (SCALAR)(arg->be)
#endif

#define  _PT_A         ((TYPE *)(arg->a))
#define  _PT_LDA       (arg->la)
#define  _PT_C         ((TYPE *)(arg->c))
#define  _PT_LDC       (arg->lc)

   Mjoin( PATL, geadd )( _PT_M, _PT_N, _PT_ALPHA, _PT_A, _PT_LDA,
                         _PT_BETA, _PT_C, _PT_LDC );
   ATL_signal_tree( root );

   return( NULL );
/*
 * End of Mjoin( PATL, ptgeadd0 )
 */
}

void Mjoin( PATL, ptgeadd )
(
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   const TYPE                 * A,
   const int                  LDA,
   const SCALAR               BETA,
   TYPE                       * C,
   const int                  LDC
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptgeadd ) adds an m-by-n matrix A to the matrix B.
 *
 * This is a multi-threaded version of the algorithm.
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
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   pthread_attr_t             attr;
   PT_TREE_T                  root = NULL;
#ifdef TREAL
   TYPE                       alpha0 = (TYPE)(ALPHA),
                              beta0  = (TYPE)(BETA);
#endif
   void                       * alpha, * beta;
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( M <= 0 ) || ( N <= 0 ) ||
       ( SCALAR_IS_ZERO( ALPHA ) && SCALAR_IS_ONE( BETA ) ) ) return;

#ifdef TREAL
   alpha = (void *)(&alpha0); beta = (void *)(&beta0);
#else
   alpha = (void *)(ALPHA);   beta = (void *)(BETA);
#endif
   ATL_thread_init( &attr );
   root = Mjoin( PATL, ptgeadd_nt )( ATL_NTHREADS, &attr, M, N, alpha,
                                     (void *)(A), LDA, beta, (void *)(C),
                                     LDC );
   ATL_join_tree  ( root );
   ATL_free_tree  ( root );
   ATL_thread_exit( &attr );
/*
 * End of Mjoin( PATL, ptgeadd )
 */
}

PT_TREE_T Mjoin( PATL, ptgeadd_nt )
(
   const unsigned int         THREADS,
   pthread_attr_t             * ATTR,
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
 * Mjoin( PATL, ptgeadd_nt ) adds an m-by-n matrix A to the matrix B.
 *
 * This is a multi-threaded version of the algorithm.
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
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = NULL;
   PT_MISC_TYPE_T             type;
#ifdef TREAL
   TYPE                       alpha,    beta;
#else
   TYPE                       alpha[2], beta[2];
#endif
   double                     tblks;
   unsigned int               nthreads;
   int                        nb, nbm1;
/* ..
 * .. Executable Statements ..
 *
 */
#ifdef TREAL
   alpha    = *((TYPE *)(ALPHA));   beta     = *((TYPE *)(BETA));
#else
   alpha[0] = ((TYPE *)(ALPHA))[0]; alpha[1] = ((TYPE *)(ALPHA))[1];
   beta [0] = ((TYPE *)(BETA ))[0]; beta [1] = ((TYPE *)(BETA ))[1];
#endif
/*
 * Make sure we don't thread this for the time being
 */
   if( THREADS >= 1 )
   {
      Mjoin( PATL, geadd )( M, N, alpha, (TYPE *)(A), LDA, beta,
                            (TYPE *)(C), LDC );
      return( root );
   }

   nbm1  = ( nb = Mjoin( PATL, GetNB )() ) - 1;
   tblks = (double)( (M+nbm1) / nb ) * (double)( (N+nbm1) / nb );

   if( ( tblks <= (double)(ATL_XOVER_MI_DEFAULT) ) || ( THREADS <= 1 ) )
   {
      Mjoin( PATL, geadd )( M, N, alpha, (TYPE *)(A), LDA, beta,
                            (TYPE *)(C), LDC );
      return( root );
   }

   if( tblks >= (double)(THREADS) ) { nthreads = THREADS; }
   else    { nthreads = (unsigned int)floor( tblks + 0.5 ); }
#ifdef TREAL
   type.size = sizeof( TYPE    ); type.fun = Mjoin( PATL, ptgeadd0 );
#else
   type.size = sizeof( TYPE[2] ); type.fun = Mjoin( PATL, ptgeadd0 );
#endif
   root = ATL_Sgeadd( &type, 0, nthreads, ATTR, nb, M, N, ALPHA, A, LDA,
                      BETA, C, LDC );
   ATL_thread_tree( root, ATTR );

   return( root );
/*
 * End of Mjoin( PATL, ptgeadd_nt )
 */
}
