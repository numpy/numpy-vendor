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

PT_TREE_T ATL_init_node
(
   unsigned int               NODE,
   PT_TREE_T                  LEFT,
   PT_TREE_T                  RIGHT,
   PT_DATA_T                  DATA,
   PT_FUN_VAL_T               * VAL,
   PT_FUN_T                   FUN,
   PT_FUN_ARG_T               ARG
)
{
/*
 * Purpose
 * =======
 *
 * ATL_init_node creates and returns a node of the recursion tree.
 *
 * Arguments
 * =========
 *
 * NODE    (input)                       unsigned int
 *         On entry, NODE specifies the node number.
 *
 * LEFT    (input)                       PT_TREE_T
 *         On entry, LEFT specifies the left child of the current node.
 *
 * RIGHT   (input)                       PT_TREE_T
 *         On entry, RIGHT specifies  the  right  child  of the  current
 *         node.
 *
 * DATA    (input)                       PT_DATA_T
 *         On entry, DATA  is a pointer to the data associated with this
 *         current node.
 *
 * VAL     (input)                       PT_FUN_VAL_T
 *         On entry, VAL  is  a pointer to the memory location where the
 *         the  function  associated  with  this  node should return its
 *         value.
 *
 * FUN     (input)                       PT_FUN_T
 *         On entry, FUN  is  a  pointer to the function associated with
 *         this current node.
 *
 * ARG     (input)                       PT_FUN_ARG_T
 *         On entry, ARG  is  a  pointer to  data  structure  containing
 *         the arguments to be passed to the function FUN.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  t = NULL;
/* ..
 * .. Executable Statements ..
 *
 */
   t = (PT_TREE_T)malloc( sizeof( PT_NODE_T ) ); ATL_assert( t != NULL );

   ATL_assert(!pthread_mutex_init( &(t->mutex), NULL ));
   ATL_assert(!pthread_cond_init ( &(t->cond),  NULL ));
   t->left  = LEFT; t->right = RIGHT; t->node = NODE; t->count = 0;
   t->data  = DATA; t->val   = VAL;   t->fun  = FUN;  t->arg   = ARG;

   return( t );
/*
 * End of ATL_init_node
 */
}
