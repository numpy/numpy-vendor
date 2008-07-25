
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

void Mjoin( PATL, ptl3settype )
(
   PT_LVL3_TYPE_T             * PLVL3
)
{
/*
 * .. Local Variables ..
 */
   static TYPE                scalars[6] = { ATL_rnone, ATL_rzero,
                                             ATL_rzero, ATL_rzero,
                                             ATL_rone,  ATL_rzero };
/* ..
 * .. Executable Statements ..
 *
 */
#ifdef TREAL
   PLVL3->size     = sizeof( TYPE    );
#else
   PLVL3->size     = sizeof( TYPE[2] );
#endif
   PLVL3->negone   = (void *)(&(scalars[0]));
   PLVL3->zero     = (void *)(&(scalars[2]));
   PLVL3->one      = (void *)(&(scalars[4]));

   PLVL3->geadd0   = Mjoin( PATL, ptgeadd0    );
   PLVL3->gemm0    = Mjoin( PATL, ptgemm0     );
   PLVL3->symm0    = Mjoin( PATL, ptsymm0     );
#ifdef TREAL
   PLVL3->hemm0    = Mjoin( PATL, ptsymm0     );
#else
   PLVL3->hemm0    = Mjoin( PATL, pthemm0     );
#endif
   PLVL3->syrk0    = Mjoin( PATL, ptsyrk0     );
   PLVL3->syr2k0   = Mjoin( PATL, ptsyr2k0    );
#ifdef TREAL
   PLVL3->herk0    = Mjoin( PATL, ptsyrk0     );
   PLVL3->her2k0   = Mjoin( PATL, ptsyr2k0    );
#else
   PLVL3->herk0    = Mjoin( PATL, ptherk0     );
   PLVL3->her2k0   = Mjoin( PATL, pther2k0    );
#endif
   PLVL3->trmm0    = Mjoin( PATL, pttrmm0     );
   PLVL3->trsm0    = Mjoin( PATL, pttrsm0     );

   PLVL3->ptgemm   = Mjoin( PATL, ptgemm_nt   );
   PLVL3->pttrmm   = Mjoin( PATL, pttrmm_nt   );
#ifdef TREAL
   PLVL3->ptsyr2k0 = Mjoin( PATL, ptsyr2k0_nt );
   PLVL3->pther2k0 = Mjoin( PATL, ptsyr2k0_nt );
#else
   PLVL3->ptsyr2k0 = Mjoin( PATL, ptsyr2k0_nt );
   PLVL3->pther2k0 = Mjoin( PATL, pther2k0_nt );
#endif

/*
 * End of Mjoin( PATL, ptl3settype )
 */
}
