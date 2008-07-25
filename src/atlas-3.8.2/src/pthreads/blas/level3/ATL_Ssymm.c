
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

PT_TREE_T ATL_Ssymm
(
   const PT_LVL3_TYPE_T       * PLVL3,
   const unsigned int         NODE,
   const unsigned int         NTHREADS,
   pthread_attr_t             * ATTR,
   const int                  NB,
   const enum ATLAS_TRANS     TGEMM,
   const enum ATLAS_SIDE      SIDE,
   const enum ATLAS_UPLO      UPLO,
   const int                  M,
   const int                  N,
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
   PT_SYMM_ARGS_T             * a_mm = NULL;
   double                     rm, rn;
   DIM_1DSPLIT_T              splitM, splitN;
   size_t                     size;
   unsigned int               next = ( NODE << 1 ), ntm1, ntm2, ntn1, ntn2;
   int                        incB, incC, m1, m2, n1, n2;
/* ..
 * .. Executable Statements ..
 *
 */
   if( NTHREADS < 1 ) return( tree );

   size = PLVL3->size;

   if( SIDE == AtlasLeft )
   {
      splitN = ATL_1dsplit( NTHREADS, N, NB, &ntn1, &ntn2, &n1, &n2, &rn );

      if( splitN == Atlas1dSplit )
      {
         left  = ATL_Ssymm( PLVL3, next+1, ntn1, ATTR, NB, TGEMM, SIDE, UPLO,
                            M, n1, ALPHA, A, LDA, B, LDB, BETA, C, LDC );
         incB  = n1 * LDB;
         incC  = n1 * LDC;
         right = ATL_Ssymm( PLVL3, next+2, ntn2, ATTR, NB, TGEMM, SIDE, UPLO,
                            M, n2, ALPHA, A, LDA, Mvpt3( B, incB, size ), LDB,
                            BETA, Mvpt3( C, incC, size ), LDC );
         tree  = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
      }
      else
      {
         a_mm = (PT_SYMM_ARGS_T *)malloc( sizeof( PT_SYMM_ARGS_T ) );
         ATL_assert( a_mm != NULL );
         a_mm->si = SIDE; a_mm->up = UPLO;
         a_mm->m  = M;    a_mm->n  = N;    a_mm->al = ALPHA;
         a_mm->a  = A;    a_mm->la = LDA;  a_mm->b  = B;     a_mm->lb = LDB;
         a_mm->be = BETA; a_mm->c  = C;    a_mm->lc = LDC;
         tree = ATL_init_node( NODE, NULL, NULL, NULL, NULL,
                               ( TGEMM == AtlasConjTrans ? PLVL3->hemm0 :
                               PLVL3->symm0 ), (void *)(a_mm) );
      }
   }
   else
   {
      splitM = ATL_1dsplit( NTHREADS, M, NB, &ntm1, &ntm2, &m1, &m2, &rm );

      if( splitM == Atlas1dSplit )
      {
         left  = ATL_Ssymm( PLVL3, next+1, ntm1, ATTR, NB, TGEMM, SIDE, UPLO,
                            m1, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC );
         incB  = m1;
         incC  = m1;
         right = ATL_Ssymm( PLVL3, next+2, ntm2, ATTR, NB, TGEMM, SIDE, UPLO,
                            m2, N, ALPHA, A, LDA, Mvpt3( B, incB, size ), LDB,
                            BETA, Mvpt3( C, incC, size ), LDC );
         tree  = ATL_init_node( NODE, left, right, NULL, NULL, NULL, NULL );
      }
      else
      {
         a_mm = (PT_SYMM_ARGS_T *)malloc( sizeof( PT_SYMM_ARGS_T ) );
         ATL_assert( a_mm != NULL );
         a_mm->si = SIDE; a_mm->up = UPLO;
         a_mm->m  = M;    a_mm->n  = N;    a_mm->al = ALPHA;
         a_mm->a  = A;    a_mm->la = LDA;  a_mm->b  = B;     a_mm->lb = LDB;
         a_mm->be = BETA; a_mm->c  = C;    a_mm->lc = LDC;
         tree = ATL_init_node( NODE, NULL, NULL, NULL, NULL,
                               ( TGEMM == AtlasConjTrans ? PLVL3->hemm0 :
                               PLVL3->symm0 ), (void *)(a_mm) );
      }
   }

   return( tree );
/*
 * End of ATL_Ssymm
 */
}
