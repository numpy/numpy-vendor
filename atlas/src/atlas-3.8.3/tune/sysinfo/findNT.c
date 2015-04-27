#ifdef ATL_NCPU

#include "atlas_misc.h"
#include "assert.h"

void PrintUsage(char *nam)
{
   fprintf(stderr, "\nUSAGE: %s [-o <outfile>]\n", nam);
   exit(-1);
}

void GetFlags(int nargs, char **args, FILE **fpout)
{
   int i;

   *fpout=stdout;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'o':
         *fpout = fopen(args[++i], "w");
         assert(*fpout);
         break;
      default:
         PrintUsage(args[0]);
      }
   }
}

main(int nargs, char **args)
{
   FILE *fpout;
   GetFlags(nargs, args, &fpout);

   fprintf(fpout, "#ifndef ATLAS_NTHREADS_H\n   #define ATLAS_NTHREADS_H\n\n");
   fprintf(fpout, "   #include \"pthread.h\"\n");
   #if ATL_NCPU != 0
      fprintf(fpout, "   #define ATL_NTHREADS %d\n", ATL_NCPU);
   #else
      fprintf(fpout, "   #define ATL_NTHREADS 4\n");
   #endif
   fprintf(fpout, "\n#endif\n");
   fclose(fpout);
   exit(0);
}
#endif
