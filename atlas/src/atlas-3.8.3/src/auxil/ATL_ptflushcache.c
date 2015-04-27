/*
 *             Automatically Tuned Linear Algebra Software v3.8.3
 *                    (C) Copyright 2000 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include <pthread.h>
#define DREAL
#include "atlas_misc.h"
#include "atlas_pthreads.h"

typedef struct ATL_FC ATL_FC;
struct ATL_FC
{
   double dret;
   double *dp;
   int N;
};

void *ATL_ptCF(void *vptr)
{
   ATL_FC *tp = vptr;
   int i, N;
   double *cache, dret=0.0;

   N = tp->N;
   if (N > 0)
   {
      cache = tp->dp;
      for (i=0; i != N; i++) dret += cache[i];
   }
   tp->dret = dret;
   return(NULL);
}

double ATL_ptflushcache(int size)
/*
 * flush cache by reading enough mem; note that if the compiler gets
 * really smart, may be necessary to make vp a global variable so it
 * can't figure out it's not being modified other than during setup;
 * the fact that ATL_dzero is external will confuse most compilers
 */
{
  static void *vp=NULL;
  static double *cache=NULL;
  double dret=0.0;
  static int i, N = 0;
  static pthread_attr_t attr;
  static pthread_t tp[ATL_NTHREADS];
  static ATL_FC fct[ATL_NTHREADS];

  if (size < 0) /* flush cache */
  {
     ATL_assert(cache);
     ATL_thread_init(&attr);
     for (i=0; i < ATL_NTHREADS-1; i++)
     {
        fct[i].N = N;
        fct[i].dp = cache+i*N;
        pthread_create(tp+i, &attr, ATL_ptCF, fct+i);
     }
     fct[i].N = N;
     fct[i].dp = cache+i*N;
     ATL_ptCF(fct+i);
     dret = fct[i].dret;
     for (i=0; i < ATL_NTHREADS-1; i++)
     {
        pthread_join(tp[i], NULL);
        dret += fct[i].dret;
     }
     ATL_thread_exit(&attr);
  }
  else if (size > 0) /* initialize */
  {
     vp = malloc(ATL_Cachelen + (size * ATL_NTHREADS));
     ATL_assert(vp);
     cache = ATL_AlignPtr(vp);
     N = size / sizeof(double);
     ATL_dzero(N*ATL_NTHREADS, cache, 1);
  }
  else if (size == 0) /* free cache */
  {
     if (vp) free(vp);
     vp = cache = NULL;
     N = 0;
  }
  return(dret);
}
