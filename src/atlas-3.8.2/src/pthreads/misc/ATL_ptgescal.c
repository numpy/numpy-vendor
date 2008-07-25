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

PT_FUN_ARG_T Mjoin( PATL, ptgescal0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = (PT_TREE_T)(ARGS);
   PT_GESCAL_ARGS_T           * arg;
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_wait_tree( root );
   arg = (PT_GESCAL_ARGS_T *)(root->arg);

#define  _PT_M         (arg->m)
#define  _PT_N         (arg->n)

#ifdef TREAL
#define  _PT_ALPHA     (SCALAR)(*((TYPE *)(arg->al)))
#else
#define  _PT_ALPHA     (SCALAR)(arg->al)
#endif

#define  _PT_A         ((TYPE *)(arg->a))
#define  _PT_LDA       (arg->la)

   Mjoin( PATL, gescal )( _PT_M, _PT_N, _PT_ALPHA, _PT_A, _PT_LDA );
   ATL_signal_tree( root );

   return( NULL );
/*
 * End of Mjoin( PATL, ptgescal0 )
 */
}

void Mjoin( PATL, ptgescal )
(
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   TYPE                       * A,
   const int                  LDA
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptgescal ) scales an m-by-n matrix A by the scalar alpha.
 *
 * This is a multi-threaded version of the algorithm.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   pthread_attr_t             attr;
   PT_TREE_T                  root = NULL;
#ifdef TREAL
   TYPE                       alpha0 = (TYPE)(ALPHA);
#endif
   void                       * alpha;
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( M == 0 ) || ( N == 0 ) || ( SCALAR_IS_ONE( ALPHA ) ) ) return;

#ifdef TREAL
   alpha = (void *)(&alpha0);
#else
   alpha = (void *)(ALPHA);
#endif
   ATL_thread_init( &attr );
   root = Mjoin( PATL, ptgescal_nt )( ATL_NTHREADS, &attr, M, N, alpha,
                                      (void *)(A), LDA );
   ATL_join_tree  ( root );
   ATL_free_tree  ( root );
   ATL_thread_exit( &attr );
/*
 * End of Mjoin( PATL, ptgescal )
 */
}

PT_TREE_T Mjoin( PATL, ptgescal_nt )
(
   const unsigned int         THREADS,
   pthread_attr_t             * ATTR,
   const int                  M,
   const int                  N,
   const void                 * ALPHA,
   void                       * A,
   const int                  LDA
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptgescal_nt ) scales an m-by-n matrix A by the scalar alpha.
 *
 * This is a multi-threaded version of the algorithm.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = NULL;
   PT_MISC_TYPE_T             type;
#ifdef TREAL
   TYPE                       alpha;
#else
   TYPE                       alpha[2];
#endif
   double                     tblks;
   unsigned int               nthreads;
   int                        nb, nbm1;
/* ..
 * .. Executable Statements ..
 *
 */
#ifdef TREAL
   alpha    = *((TYPE *)(ALPHA));
#else
   alpha[0] = ((TYPE *)(ALPHA))[0]; alpha[1] = ((TYPE *)(ALPHA))[1];
#endif
/*
 * Make sure we don't thread this for the time being
 */
   if( THREADS >= 1 )
   { Mjoin( PATL, gescal )( M, N, alpha, (TYPE *)(A), LDA ); return( root ); }

   nbm1  = ( nb = Mjoin( PATL, GetNB )() ) - 1;
   tblks = (double)( (M+nbm1) / nb ) * (double)( (N+nbm1) / nb );
/*
 * When one or more dimension is less than nb, do not thread.
 */
   if( ( M <= nb ) || ( N <= nb ) || ( THREADS <= 1 ) ||
       ( tblks <= (double)(ATL_XOVER_MI_DEFAULT) ) )
   { Mjoin( PATL, gescal )( M, N, alpha, (TYPE *)(A), LDA ); return( root ); }

   if( tblks >= (double)(THREADS) ) { nthreads = THREADS; }
   else    { nthreads = (unsigned int)floor( tblks + 0.5 ); }
#ifdef TREAL
   type.size = sizeof( TYPE    ); type.fun = Mjoin( PATL, ptgescal0 );
#else
   type.size = sizeof( TYPE[2] ); type.fun = Mjoin( PATL, ptgescal0 );
#endif
   root = ATL_Sgescal( &type, 0, nthreads, ATTR, nb, M, N, ALPHA, A, LDA );
   ATL_thread_tree( root, ATTR );

   return( root );
/*
 * End of Mjoin( PATL, ptgescal_nt )
 */
}
