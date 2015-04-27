
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

PT_FUN_ARG_T Mjoin( PATL, ptherk0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = (PT_TREE_T)(ARGS);
   PT_SYRK_ARGS_T             * arg;
   int                        k, l, la, lc, m, n;
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_wait_tree( root );
   arg = (PT_SYRK_ARGS_T *)(root->arg);

#define  _PT_ALPHA     (SCALAR)(arg->al)
#define  _PT_BETA      (SCALAR)(arg->be)
#define  _PT_ALPHA_R   *((TYPE *)(arg->al))
#define  _PT_BETA_R    *((TYPE *)(arg->be))

#define  _PT_GEMM_Nx   Mjoin( PATL, gemmNC )
#define  _PT_GEMM_xN   Mjoin( PATL, gemmCN )

#define  _PT_A         ((TYPE *)(arg->a))
#define  _PT_C         ((TYPE *)(arg->c))

   k = arg->k; l = arg->l; m = arg->m; n = arg->n; la = arg->la; lc = arg->lc;

   if( arg->up == AtlasLower )
   {
      if( arg->tr == AtlasNoTrans )
      {
         if( ( m + k > 0 ) && ( n > 0 ) )
         {
            _PT_GEMM_Nx( m + k, n, l, _PT_ALPHA, _PT_A, la, _PT_A - (n SHIFT),
                         la, _PT_BETA, _PT_C, lc );
         }
         Mjoin( PATL, herk )( AtlasLower, AtlasNoTrans, k, l, _PT_ALPHA_R,
                              _PT_A, la, _PT_BETA_R, _PT_C + ((n * lc) SHIFT),
                              lc );
         if( ( m > 0 ) && ( k > 0 ) )
         {
            _PT_GEMM_Nx( m, k, l, _PT_ALPHA, _PT_A + (k SHIFT), la, _PT_A, la,
                         _PT_BETA, _PT_C + ((k + n * lc) SHIFT), lc );
         }
      }
      else
      {
         if( ( m + k > 0 ) && ( n > 0 ) )
         {
            _PT_GEMM_xN( m + k, n, l, _PT_ALPHA, _PT_A + ((n * la) SHIFT),
                         la, _PT_A, la, _PT_BETA, _PT_C, lc );
         }
         Mjoin( PATL, herk )( AtlasLower, AtlasTrans, k, l, _PT_ALPHA_R,
                              _PT_A + ((n * la) SHIFT), la, _PT_BETA_R,
                              _PT_C + ((n * lc) SHIFT), lc );
         if( ( m > 0 ) && ( k > 0 ) )
         {
            _PT_GEMM_xN( m, k, l, _PT_ALPHA, _PT_A + (( (n + k) * la) SHIFT),
                         la, _PT_A + ((n * la) SHIFT), la, _PT_BETA, _PT_C +
                         ((k + n * lc) SHIFT), lc );
         }
      }
   }
   else
   {
      if( arg->tr == AtlasNoTrans )
      {
         if( ( m > 0 ) && ( k > 0 ) )
         {
            _PT_GEMM_Nx( m, k, l, _PT_ALPHA, _PT_A, la, _PT_A + (m SHIFT),
                         la, _PT_BETA, _PT_C, lc );
         }
         Mjoin( PATL, herk )( AtlasUpper, AtlasNoTrans, k, l, _PT_ALPHA_R,
                              _PT_A + (m SHIFT), la, _PT_BETA_R,
                              _PT_C + (m SHIFT), lc );
         if( ( m + k > 0 ) && ( n > 0 ) )
         {
            _PT_GEMM_Nx( m + k, n, l, _PT_ALPHA, _PT_A, la,
                         _PT_A + ((m + k ) SHIFT), la, _PT_BETA,
                         _PT_C + ((k * lc) SHIFT), lc );
         }
      }
      else
      {
         if( ( m > 0 ) && ( k > 0 ) )
         {
            _PT_GEMM_xN( m, k, l, _PT_ALPHA, _PT_A - ((m * la) SHIFT), la,
                         _PT_A, la, _PT_BETA, _PT_C, lc );
         }
         Mjoin( PATL, herk )( AtlasUpper, AtlasTrans, k, l, _PT_ALPHA_R,
                              _PT_A, la, _PT_BETA_R, _PT_C + (m SHIFT), lc );
         if( ( m + k > 0 ) && ( n > 0 ) )
         {
            _PT_GEMM_xN( m + k, n, l, _PT_ALPHA, _PT_A - ((m * la) SHIFT),
                         la, _PT_A + ((k * la) SHIFT), la, _PT_BETA,
                         _PT_C + ((k * lc) SHIFT), lc );
         }
      }
   }

   ATL_signal_tree( root );

   return( NULL );
/*
 * End of Mjoin( PATL, ptherk0 )
 */
}

void Mjoin( PATL, ptherk )
(
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const int                  N,
   const int                  K,
   const TYPE                 ALPHA,
   const TYPE                 * A,
   const int                  LDA,
   const TYPE                 BETA,
   TYPE                       * C,
   const int                  LDC
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptherk )  performs one of the Hermitian rank k operations
 *
 *    C := alpha * A * conjg( A' ) + beta * C,
 *
 * or
 *
 *    C := alpha * conjg( A' ) * A + beta * C,
 *
 * where alpha and beta are  real  scalars, C is an n by n Hermitian ma-
 * trix and  A is an n by k matrix in the first case and a k by n matrix
 * in the second case.
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
   TYPE                       Calph[2], Cbeta[2];
   void                       * alpha, * beta;
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( N == 0 ) ||
       ( ( ( ALPHA == ATL_rzero ) || ( K == 0 ) ) && ( BETA == ATL_rone ) ) )
      return;

   if( ( ( ALPHA == ATL_rzero ) ) || ( K == 0 ) )
   { Mjoin( PATL, pthescal )( UPLO, N, N, BETA, C, LDC ); return; }

   ATL_thread_init( &attr );
   *Calph = ALPHA; Calph[1] = ATL_rzero; alpha  = (void *)(Calph);
   *Cbeta = BETA;  Cbeta[1] = ATL_rzero; beta   = (void *)(Cbeta);
   root = Mjoin( PATL, ptherk_nt )( ATL_NTHREADS, &attr, UPLO, TRANS,
                                    N, K, alpha, (void *)(A), LDA, beta,
                                    (void *)(C), LDC );
   ATL_join_tree  ( root );
   ATL_free_tree  ( root );
   ATL_thread_exit( &attr );
/*
 * End of Mjoin( PATL, ptherk )
 */
}

PT_TREE_T Mjoin( PATL, ptherk_nt )
(
   const unsigned int         NTHREADS,
   pthread_attr_t             * ATTR,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const int                  N,
   const int                  K,
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
 * Mjoin( PATL, ptherk_nt )  performs one of the Hermitian rank k operations
 *
 *    C := alpha * A * conjg( A' ) + beta * C,
 *
 * or
 *
 *    C := alpha * conjg( A' ) * A + beta * C,
 *
 * where alpha and beta are  real  scalars, C is an n by n Hermitian ma-
 * trix and  A is an n by k matrix in the first case and a k by n matrix
 * in the second case.
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
   nbm1  = ( nb = Mjoin( PATL, GetNB )() ) - 1;
   tblks = (double)( (N+nbm1) / nb ) * (double)( (N+nbm1) / nb ) *
           (double)( (K+nbm1) / nb );

   if( ( tblks <= (double)(ATL_XOVER_L3_DEFAULT) ) || ( NTHREADS <= 1 ) )
   {
#ifdef TREAL
      alpha    = *((TYPE *)(ALPHA));   beta     = *((TYPE *)(BETA));
#else
      alpha[0] = ((TYPE *)(ALPHA))[0]; alpha[1] = ((TYPE *)(ALPHA))[1];
      beta [0] = ((TYPE *)(BETA ))[0]; beta [1] = ((TYPE *)(BETA ))[1];
#endif
      Mjoin( PATL, herk )( UPLO, TRANS, N, K, alpha[0], (TYPE *)(A), LDA,
                           beta[0], (TYPE *)(C), LDC );
      return( root );
   }

   if( tblks >= (double)(NTHREADS) ) { nthreads = NTHREADS; }
   else    { nthreads = (unsigned int)floor( tblks + 0.5 ); }

   Mjoin( PATL, ptl3settype )( &type );
   root = ATL_Ssyrk( &type, 0, nthreads, ATTR, nb, UPLO, TRANS, AtlasConjTrans,
                     0, 0, N, K, ALPHA, A, LDA, BETA, C, LDC );
   ATL_thread_tree( root, ATTR );

   return( root );
/*
 * End of Mjoin( PATL, ptherk_nt )
 */
}
