#include "atlconf.h"

void GetPmake(int verb, enum OSTYPE OS, int ncpu, char *pmake)
{
   if (ncpu == 0) ncpu = 1;
   if (!OSIsWin(OS)&&CmndResults(NULL, "make DoNothing --version | fgrep GNU"))
   { /* using gnu make */
      if (ncpu == 1) sprintf(pmake, "$(MAKE) -j 2");
      else sprintf(pmake, "$(MAKE) -j %d", ncpu);
   }
   else
   {  /* AIX, HP-UX, SunOS make do not have parallel option */
      switch(OS)
      {
      case OSIRIX:
         strcpy(pmake, "$(MAKE) -P");
         break;
      case OSFreeBSD:
         if (ncpu == 1) sprintf(pmake, "$(MAKE) -j 2");
         else sprintf(pmake, "$(MAKE) -j %d", ncpu);
         break;
      default:
         strcpy(pmake, "$(MAKE)");
      }
   }
}

void PrintUsage(char *name, int iarg, char *arg)
{
   fprintf(stderr, "\nERROR around arg %d (%s).\n", iarg,
           arg ? arg : "unknown");
   fprintf(stderr, "USAGE: %s [flags] where flags are:\n", name);
   fprintf(stderr, "   -v <verb> : verbosity level\n");
   fprintf(stderr, "   -O <enum OSTYPE #>  : set OS type\n");
   fprintf(stderr,
           "   -t <#> : set # of threads (-1: autodect; 0: no threading)\n");
   fprintf(stderr,
      "NOTE: enum #s can be found by : make xprint_enums ; ./xprint_enums\n");
   exit(iarg);
}

void GetFlags(int nargs,                /* nargs as passed into main */
              char **args,              /* args as passed into main */
              int *verb,                /* verbosity setting */
              enum OSTYPE *OS,          /* OS to assume */
              int *nthreads,           /* # of threads */
              char **targ             /* mach to ssh to*/
             )
{
   int i, k, k0, kn, DoInt;
   char *sp, *sp0;

   *verb = 0;
   *targ = NULL;

   *nthreads = 0;
   *OS = 0;
   *verb = 0;
   *nthreads = -1;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 't':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *nthreads = atoi(args[i]);
         break;
      case 'O':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *OS = atoi(args[i]);
         break;
      case 'v':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *verb = atoi(args[i]);
         break;
      case 'T':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *targ = args[i];
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
}

main(int nargs, char **args)
{
   int verb, OS, arch, ncpu;
   char *targ;
   char pmake[256];
   GetFlags(nargs, args, &verb, &OS, &ncpu, &targ);
   GetPmake(verb, OS, ncpu, pmake);
   printf("PMAKE='%s'\n", pmake);
   exit(0);
}
