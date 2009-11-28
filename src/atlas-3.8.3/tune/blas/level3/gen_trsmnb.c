#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

main(int nargs, char **args)
{
   FILE *fpin, *fpout;
   int stops, stopd;
   fpin = fopen("res/sTRSM_NB", "r");
   assert(fpin);
   assert(fscanf(fpin, "%d", &stops) == 1);
   fclose(fpin);
   fpin = fopen("res/dTRSM_NB", "r");
   assert(fpin);
   assert(fscanf(fpin, "%d", &stopd) == 1);
   fclose(fpin);
   fpout = fopen("res/atlas_trsmNB.h", "w");
   assert(fpout);
   fprintf(fpout, "#ifndef ATLAS_TRSMNB_H\n   #define ATLAS_TRSMNB_H\n\n");
   fprintf(fpout, "   #ifdef SREAL\n      #define TRSM_NB %d\n", stops);
   fprintf(fpout, "   #elif defined(DREAL)\n      #define TRSM_NB %d\n", stopd);
   fprintf(fpout, "   #else\n      #define TRSM_NB 4\n");
   fprintf(fpout, "   #endif\n\n#endif\n");
   fclose(fpout);
   exit(0);
}
