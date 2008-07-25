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

PT_FUN_ARG_T Mjoin( PATL, ptgezero0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = (PT_TREE_T)(ARGS);
   PT_GEZERO_ARGS_T           * arg;
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_wait_tree( root );
   arg = (PT_GEZERO_ARGS_T *)(root->arg);

#define  _PT_M         (arg->m)
#define  _PT_N         (arg->n)
#define  _PT_A         ((TYPE *)(arg->a))
#define  _PT_LDA       (arg->la)

   Mjoin( PATL, gezero )( _PT_M, _PT_N, _PT_A, _PT_LDA );
   ATL_signal_tree( root );

   return( NULL );
/*
 * End of Mjoin( PATL, ptgezero0 )
 */
}

void Mjoin( PATL, ptgezero )
(
   const int                  M,
   const int                  N,
   TYPE                       * A,
   const int                  LDA
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptgezero ) zeroes an m-by-n matrix A.
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
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( M <= 0 ) || ( N <= 0 ) ) return;

   ATL_thread_init( &attr );
   root = Mjoin( PATL, ptgezero_nt )( ATL_NTHREADS, &attr, M, N, A, LDA );
   ATL_join_tree  ( root );
   ATL_free_tree  ( root );
   ATL_thread_exit( &attr );
/*
 * End of Mjoin( PATL, ptgezero )
 */
}

PT_TREE_T Mjoin( PATL, ptgezero_nt )
(
   const unsigned int         THREADS,
   pthread_attr_t             * ATTR,
   const int                  M,
   const int                  N,
   void                       * A,
   const int                  LDA
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptgezero_nt ) zeroes an m-by-n matrix A.
 *
 * This is a multi-threaded version of the algorithm.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   PT_MISC_TYPE_T             type;
   PT_TREE_T                  root = NULL;
   double                     tblks;
   unsigned int               nthreads;
   int                        nb, nbm1;
/* ..
 * .. Executable Statements ..
 *
 */
/*
 * Make sure we don't thread this for the time being
 */
   if( THREADS >= 1 )
   { Mjoin( PATL, gezero )( M, N, (TYPE *)(A), LDA ); return( root ); }

   nbm1  = ( nb = Mjoin( PATL, GetNB )() ) - 1;
   tblks = (double)( (M+nbm1) / nb ) * (double)( (N+nbm1) / nb );
/*
 * When one or more dimension is less than nb, do not thread.
 */
   if( ( M <= nb ) || ( N <= nb ) || ( THREADS <= 1 ) ||
       ( tblks <= (double)(ATL_XOVER_MI_DEFAULT) ) )
   { Mjoin( PATL, gezero )( M, N, (TYPE *)(A), LDA ); return( root ); }

   if( tblks >= (double)(THREADS) ) { nthreads = THREADS; }
   else    { nthreads = (unsigned int)floor( tblks + 0.5 ); }
#ifdef TREAL
   type.size = sizeof( TYPE    ); type.fun = Mjoin( PATL, ptgezero0 );
#else
   type.size = sizeof( TYPE[2] ); type.fun = Mjoin( PATL, ptgezero0 );
#endif
   root = ATL_Sgezero( &type, 0, nthreads, ATTR, nb, M, N, A, LDA );
   ATL_thread_tree( root, ATTR );

   return( root );
/*
 * End of Mjoin( PATL, ptgezero_nt )
 */
}
