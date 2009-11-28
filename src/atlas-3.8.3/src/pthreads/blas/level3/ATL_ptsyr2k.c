
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
#include "atlas_lvl3.h"
#include "atlas_kern3.h"

PT_FUN_ARG_T Mjoin( PATL, ptsyr2k0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
#ifdef TREAL
   TYPE                       one    = ATL_rone;
#else
   TYPE                       one[2] = { ATL_rone, ATL_rzero };
#endif
   PT_TREE_T                  root = (PT_TREE_T)(ARGS);
   PT_SYR2K_ARGS_T            * arg;
   int                        k, l, la, lb, lc, m, n;
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_wait_tree( root );
   arg = (PT_SYR2K_ARGS_T *)(root->arg);

#ifdef TREAL
#define  _PT_ALPHA     (SCALAR)(*((TYPE *)(arg->al)))
#define  _PT_ALPHC     (SCALAR)(*((TYPE *)(arg->al)))
#define  _PT_BETA      (SCALAR)(*((TYPE *)(arg->be)))
#define  _PT_BETA_R    (SCALAR)(*((TYPE *)(arg->be)))
#else
#define  _PT_ALPHA     (SCALAR)(arg->al)
#define  _PT_ALPHC     (SCALAR)(arg->al)
#define  _PT_BETA      (SCALAR)(arg->be)
#define  _PT_BETA_R    (SCALAR)(arg->be)
#endif

#define  _PT_GEMM_Nx   Mjoin( PATL, gemmNT )
#define  _PT_GEMM_xN   Mjoin( PATL, gemmTN )


#define  _PT_A         ((TYPE *)(arg->a))
#define  _PT_B         ((TYPE *)(arg->b))
#define  _PT_C         ((TYPE *)(arg->c))

   k  = arg->k;  l  = arg->l;  m  = arg->m;  n  = arg->n;
   la = arg->la; lb = arg->lb; lc = arg->lc;

   if( arg->up == AtlasLower )
   {
      if( arg->tr == AtlasNoTrans )
      {
         if( ( m + k > 0 ) && ( n > 0 ) )
         {
            _PT_GEMM_Nx( m + k, n, l, _PT_ALPHA, _PT_A, la, _PT_B - (n SHIFT),
                         lb, _PT_BETA, _PT_C, lc );
            _PT_GEMM_Nx( m + k, n, l, _PT_ALPHC, _PT_B, lb, _PT_A - (n SHIFT),
                         la, one,      _PT_C, lc );
         }
         Mjoin( PATL, syr2k )( AtlasLower, AtlasNoTrans, k, l, _PT_ALPHA,
                               _PT_A, la, _PT_B, lb, _PT_BETA_R,
                               _PT_C + ((n * lc) SHIFT), lc );
         if( ( m > 0 ) && ( k > 0 ) )
         {
            _PT_GEMM_Nx( m, k, l, _PT_ALPHA, _PT_A + (k SHIFT), la, _PT_B, lb,
                         _PT_BETA, _PT_C + ((k + n * lc) SHIFT), lc );
            _PT_GEMM_Nx( m, k, l, _PT_ALPHC, _PT_B + (k SHIFT), lb, _PT_A, la,
                         one,      _PT_C + ((k + n * lc) SHIFT), lc );
         }
      }
      else
      {
         if( ( m + k > 0 ) && ( n > 0 ) )
         {
            _PT_GEMM_xN( m + k, n, l, _PT_ALPHA, _PT_A + ((n * la) SHIFT),
                         la, _PT_B, lb, _PT_BETA, _PT_C, lc );
            _PT_GEMM_xN( m + k, n, l, _PT_ALPHC, _PT_B + ((n * lb) SHIFT),
                         lb, _PT_A, la, one,      _PT_C, lc );
         }
         Mjoin( PATL, syr2k )( AtlasLower, AtlasTrans, k, l, _PT_ALPHA,
                               _PT_A + ((n * la) SHIFT), la,
                               _PT_B + ((n * lb) SHIFT), lb,
                               _PT_BETA_R, _PT_C + ((n * lc) SHIFT), lc );
         if( ( m > 0 ) && ( k > 0 ) )
         {
            _PT_GEMM_xN( m, k, l, _PT_ALPHA, _PT_A + (( (n + k) * la) SHIFT),
                         la, _PT_B + ((n * lb) SHIFT), lb, _PT_BETA, _PT_C +
                         ((k + n * lc) SHIFT), lc );
            _PT_GEMM_xN( m, k, l, _PT_ALPHC, _PT_B + (( (n + k) * lb) SHIFT),
                         lb, _PT_A + ((n * la) SHIFT), la, one,      _PT_C +
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
            _PT_GEMM_Nx( m, k, l, _PT_ALPHA, _PT_A, la, _PT_B + (m SHIFT),
                         lb, _PT_BETA, _PT_C, lc );
            _PT_GEMM_Nx( m, k, l, _PT_ALPHC, _PT_B, lb, _PT_A + (m SHIFT),
                         la, one,      _PT_C, lc );
         }
         Mjoin( PATL, syr2k )( AtlasUpper, AtlasNoTrans, k, l, _PT_ALPHA,
                               _PT_A + (m SHIFT), la,
                               _PT_B + (m SHIFT), lb, _PT_BETA_R,
                               _PT_C + (m SHIFT), lc );
         if( ( m + k > 0 ) && ( n > 0 ) )
         {
            _PT_GEMM_Nx( m + k, n, l, _PT_ALPHA, _PT_A, la,
                         _PT_B + ((m + k ) SHIFT), lb, _PT_BETA,
                         _PT_C + ((k * lc) SHIFT), lc );
            _PT_GEMM_Nx( m + k, n, l, _PT_ALPHC, _PT_B, lb,
                         _PT_A + ((m + k ) SHIFT), la, one,
                         _PT_C + ((k * lc) SHIFT), lc );
         }
      }
      else
      {
         if( ( m > 0 ) && ( k > 0 ) )
         {
            _PT_GEMM_xN( m, k, l, _PT_ALPHA, _PT_A - ((m * la) SHIFT), la,
                         _PT_B, lb, _PT_BETA, _PT_C, lc );
            _PT_GEMM_xN( m, k, l, _PT_ALPHC, _PT_B - ((m * lb) SHIFT), lb,
                         _PT_A, la, one,      _PT_C, lc );
         }
         Mjoin( PATL, syr2k )( AtlasUpper, AtlasTrans, k, l, _PT_ALPHA,
                               _PT_A, la, _PT_B, lb, _PT_BETA_R,
                               _PT_C + (m SHIFT), lc );
         if( ( m + k > 0 ) && ( n > 0 ) )
         {
            _PT_GEMM_xN( m + k, n, l, _PT_ALPHA, _PT_A - ((m * la) SHIFT),
                         la, _PT_B + ((k * lb) SHIFT), lb, _PT_BETA,
                         _PT_C + ((k * lc) SHIFT), lc );
            _PT_GEMM_xN( m + k, n, l, _PT_ALPHC, _PT_B - ((m * lb) SHIFT),
                         lb, _PT_A + ((k * la) SHIFT), la, one,
                         _PT_C + ((k * lc) SHIFT), lc );
         }
      }
   }

   ATL_signal_tree( root );

   return( NULL );
/*
 * End of Mjoin( PATL, ptsyr2k0 )
 */
}

int Mjoin( PATL, ptsyr2k0_nt )
(
   const unsigned int         NTHREADS,
   pthread_attr_t             * ATTR,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const enum ATLAS_TRANS     TGEMM,
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
#ifdef TREAL
   TYPE                       alpha, beta, zero;
#else
   TYPE                       alpha[2], beta[2], zero[2];
#endif
   PT_TREE_T                  root = NULL;
   void                       * vc = NULL, * zptr;
   TYPE                       * c;
   size_t                     lgth;

   if( ( lgth = ATL_MulBySize(N)*N ) <= ATL_MaxMalloc )
      vc = malloc( ATL_Cachelen + lgth );

   if( vc == NULL ) return( 1 );
#ifdef TREAL
   alpha    = *((TYPE *)(ALPHA)); beta = *((TYPE *)(BETA) );
   zero     = ATL_rzero;          zptr = (void *)(&zero);
#else
   alpha[0] = ((TYPE *)(ALPHA))[0]; alpha[1] = ((TYPE *)(ALPHA))[1];
   beta [0] = ((TYPE *)(BETA ))[0]; beta [1] = ((TYPE *)(BETA ))[1];
   zero [0] = ATL_rzero;            zero [1] = ATL_rzero;
   zptr     = (void *)(zero);
#endif
   c = ATL_AlignPtr( vc );
   if( TRANS == AtlasNoTrans )
   {
      root = Mjoin( PATL, ptgemm_nt )( NTHREADS, ATTR, AtlasNoTrans, TGEMM,
                                       N, N, K, ALPHA, A, LDA, B, LDB, zptr,
                                       (void *)(c), N );
   }
   else
   {
      root = Mjoin( PATL, ptgemm_nt )( NTHREADS, ATTR, TGEMM, AtlasNoTrans,
                                       N, N, K, ALPHA, A, LDA, B, LDB, zptr,
                                       (void *)(c), N );
   }
   ATL_join_tree( root );
   ATL_free_tree( root );

   if( UPLO == AtlasUpper )
   {
      if(      SCALAR_IS_ONE(  beta ) )
      { Mjoin( PATL, syr2k_putU_b1   )( N, c, beta, (TYPE *)(C), LDC ); }
      else if( SCALAR_IS_ZERO( beta ) )
      { Mjoin( PATL, syr2k_putU_b0   )( N, c, beta, (TYPE *)(C), LDC ); }
#ifdef TCPLX
      else if( SCALAR_IS_NONE( beta ) )
      { Mjoin( PATL, syr2k_putU_bn1  )( N, c, beta, (TYPE *)(C), LDC ); }
      else if( beta[1] == ATL_rzero   )
      { Mjoin( PATL, syr2k_putU_bXi0 )( N, c, beta, (TYPE *)(C), LDC ); }
#endif
      else
      { Mjoin( PATL, syr2k_putU_bX   )( N, c, beta, (TYPE *)(C), LDC ); }
   }
   else
   {
      if(      SCALAR_IS_ONE(  beta ) )
      { Mjoin( PATL, syr2k_putL_b1   )( N, c, beta, (TYPE *)(C), LDC ); }
      else if( SCALAR_IS_ZERO( beta ) )
      { Mjoin( PATL, syr2k_putL_b0   )( N, c, beta, (TYPE *)(C), LDC ); }
#ifdef TCPLX
      else if( SCALAR_IS_NONE( beta ) )
      { Mjoin( PATL, syr2k_putL_bn1  )( N, c, beta, (TYPE *)(C), LDC ); }
      else if( beta[1] == ATL_rzero   )
      { Mjoin( PATL, syr2k_putL_bXi0 )( N, c, beta, (TYPE *)(C), LDC ); }
#endif
      else
      { Mjoin( PATL, syr2k_putL_bX   )( N, c, beta, (TYPE *)(C), LDC ); }
   }

   free( vc );

   return( 0 );
/*
 * End of Mjoin( PATL, ptsyr2k0_nt )
 */
}

void Mjoin( PATL, ptsyr2k )
(
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
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
 * Mjoin( PATL, ptsyr2k ) performs one of the symmetric rank 2k operations
 *
 *    C := alpha * A * B' + alpha * B * A' + beta * C,
 *
 * or
 *
 *    C := alpha * A' * B + alpha * B' * A + beta * C,
 *
 * where alpha and beta are scalars, C is an n by n symmetric matrix and
 * A and B are n by k matrices in the first case and k by n  matrices in
 * the second case.
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
   if( ( N == 0 ) || ( ( SCALAR_IS_ZERO( ALPHA ) || ( K == 0 ) ) &&
                       SCALAR_IS_ONE( BETA ) ) ) return;

   if( ( SCALAR_IS_ZERO( ALPHA ) ) || ( K == 0 ) )
   { Mjoin( PATL, pttrscal )( UPLO, N, N, BETA, C, LDC ); return; }

   ATL_thread_init( &attr );
#ifdef TREAL
   root = Mjoin( PATL, ptsyr2k_nt )( ATL_NTHREADS, &attr, UPLO, TRANS, N,
                                     K, (void *)(&alpha0), (void *)(&alpha0),
                                     (void *)(A), LDA, (void *)(B), LDB,
                                     (void *)(&beta0), (void *)(C), LDC );
#else
   root = Mjoin( PATL, ptsyr2k_nt )( ATL_NTHREADS, &attr, UPLO, TRANS, N,
                                     K, (TYPE *)(ALPHA), (TYPE *)(ALPHA),
                                     (void *)(A), LDA, (void *)(B), LDB,
                                     (TYPE *)(BETA), (void *)(C), LDC );
#endif
   ATL_join_tree  ( root );
   ATL_free_tree  ( root );
   ATL_thread_exit( &attr );
/*
 * End of Mjoin( PATL, ptsyr2k )
 */
}

PT_TREE_T Mjoin( PATL, ptsyr2k_nt )
(
   const unsigned int         NTHREADS,
   pthread_attr_t             * ATTR,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const int                  N,
   const int                  K,
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
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptsyr2k_nt ) performs one of the symmetric rank 2k operations
 *
 *    C := alpha * A * B' + alpha * B * A' + beta * C,
 *
 * or
 *
 *    C := alpha * A' * B + alpha * B' * A + beta * C,
 *
 * where alpha and beta are scalars, C is an n by n symmetric matrix and
 * A and B are n by k matrices in the first case and k by n  matrices in
 * the second case.
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
      alpha    = *((TYPE *)(ALPHA));   beta     = *((TYPE *)(BETA ));
#else
      alpha[0] = ((TYPE *)(ALPHA))[0]; alpha[1] = ((TYPE *)(ALPHA))[1];
      beta [0] = ((TYPE *)(BETA ))[0]; beta [1] = ((TYPE *)(BETA ))[1];
#endif
      Mjoin( PATL, syr2k )( UPLO, TRANS, N, K, alpha, (TYPE *)(A), LDA,
                            (TYPE *)(B), LDB, beta, C, LDC );
      return( root );
   }

   if( tblks >= (double)(NTHREADS) ) { nthreads = NTHREADS; }
   else    { nthreads = (unsigned int)floor( tblks + 0.5 ); }

   Mjoin( PATL, ptl3settype )( &type );
   root = ATL_Ssyr2k( &type, 0, nthreads, ATTR, nb, UPLO, TRANS, AtlasTrans,
                      0, 0, N, K, ALPHA, ALPHC, A, LDA, B, LDB, BETA, C, LDC );
   ATL_thread_tree( root, ATTR );

   return( root );
/*
 * End of Mjoin( PATL, ptsyr2k_nt )
 */
}
