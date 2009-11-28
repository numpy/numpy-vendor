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

void ATL_thread_init
(
   pthread_attr_t             * ATTR
)
{
/*
 * Purpose
 * =======
 *
 * ATL_thread_init initializes the pthread library.
 *
 * By default, threads are created joinable and scheduled within  system
 * scope. From "Pthreads Programming", B. Nichols & al., O Reilly, 1996:
 *
 * A  thread  detached state can be joinable  (PTHREAD_CREATE_JOINABLE),
 * or it can be detached (PTHREAD_CREATE_DETACHED).  The following  call
 * to the   pthread_attr_setdetachstate  function  adjusts  the detached
 * state of a thread.
 *
 * A thread can use system-scope scheduling  (PTHREAD_SCOPE_SYSTEM),  in
 * which  case  the  operating  system  compares  the  priorities of all
 * runnable  threads of all processes systemwide to select a  thread  to
 * run  on an available CPU.  Alternatively,  it  can  use process-scope
 * scheduling (PTHREAD_SCOPE_PROCESS),  in  which  case only the highest
 * priority runnable  thread  in a  process competes against the highest
 * priority   threads  of  other  processes  in  the  system  scheduling
 * activity. The following function call adjusts the scheduling scope of
 * a thread.
 *
 * Arguments
 * =========
 *
 * ATTR    (output)                      pthread_attr_t *
 *         On exit,  ATTR  specifies  the  thread attribute object to be
 *         initialized.
 *
 * ---------------------------------------------------------------------
 */
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_assert(!pthread_attr_init(ATTR));
#ifdef IBM_PT_ERROR
   ATL_assert(!pthread_attr_setdetachstate( ATTR, PTHREAD_CREATE_UNDETACHED ));
#else
   ATL_assert(!pthread_attr_setdetachstate( ATTR, PTHREAD_CREATE_JOINABLE   ));
#endif
   pthread_attr_setscope(ATTR, PTHREAD_SCOPE_SYSTEM); /* no chk, OK to fail */
/*
 * End of ATL_thread_init
 */
}
