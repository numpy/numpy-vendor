#include "atlconf.h"

void PrintUsage(char *name, int iarg, char *arg)
{
   fprintf(stderr, "\nERROR around arg %d (%s).\n", iarg,
           arg ? arg : "unknown");
   fprintf(stderr, "USAGE: %s [flags] where flags are:\n", name);
   fprintf(stderr, "   -v <verb> : verbosity level\n");
   fprintf(stderr, "   -O <enum OSTYPE #>  : set OS type\n");
   fprintf(stderr, "   -s <enum ASMDIA #>  : set assembly dialect\n");
   fprintf(stderr, "   -A <enum MACHTYPE #> : set machine/architecture\n");
   fprintf(stderr,
   "   -V #    # = ((1<<vecISA1) | (1<<vecISA2) | ... | (1<<vecISAN))\n");
   fprintf(stderr, "   -b <32/64> : set pointer bitwidth\n");
   fprintf(stderr, "   -o <outfile>\n");
   fprintf(stderr, "   -C [xc,ic,if,sk,dk,sm,dm,al,ac] <compiler>\n");
   fprintf(stderr, "   -F [xc,ic,if,sk,dk,sm,dm,al,ac,gc] '<comp flags>'\n");
   fprintf(stderr,    /* HERE */
           "   -Fa [xc,ic,if,sk,dk,sm,dm,al,ac,gc] '<comp flags to append>'\n");
   fprintf(stderr, "        al: append flags to all compilers\n");
   fprintf(stderr, "        ac: append flags to all C compilers\n");
   fprintf(stderr, "        gc: append flags to gcc compiler used in user-contributed index files.\n");
   fprintf(stderr, "        acg: append to all C compilers & the index gcc\n");
   fprintf(stderr, "        alg: append to all compilers & the index gcc\n");
   fprintf(stderr,
      "   -T <targ> : ssh target for cross-compilation (probably broken)\n");
   fprintf(stderr, "   -D [c,f] -D<mac>=<rep> : cpp #define to add to [CDEFS,F2CDEFS]\n");
   fprintf(stderr,
   "      eg. -D c -DL2SIZE=8388604 -D f -DADD__ -D f -DStringSunStyle\n");
   fprintf(stderr, "   -d [s,b]  : set source/build directory\n");
   fprintf(stderr, "   -f <#> : size (in KB) to flush before timing\n");
   fprintf(stderr,
           "   -t <#> : set # of threads (-1: autodect; 0: no threading)\n");
   fprintf(stderr, "   -m <mhz> : set clock rate\n");
   fprintf(stderr, "   -S[i/s] <handle> <val>  : special int/string arg\n");
   fprintf(stderr,
           "      -Si bozol1 <0/1> : supress/enable bozo L1 defaults\n");
   fprintf(stderr,
           "      -Si archdef <1/0> : enable/supress arch default use\n");
      fprintf(stderr,
        "      -Si nof77 <0/1> : Have/don't have fortran compiler\n");
   fprintf(stderr,
        "      -Ss kern <path to comp> : use comp for all kernel compilers\n");
   fprintf(stderr,
        "      -Ss pmake <parallel make invocation (eg '$(MAKE) -j 4')>\n");
   fprintf(stderr,
"      -Ss f77lib <path to f77 lib needed by C compiler>\n");
   fprintf(stderr,
"      -Ss flapack <path to netlib lapack>: used to build full lapack lib\n");
   fprintf(stderr, "      -Ss [s,d]maflags 'flags'\n");
   fprintf(stderr,
      "NOTE: enum #s can be found by : make xprint_enums ; ./xprint_enums\n");
   exit(iarg);
}

void GetFlags(int nargs,                /* nargs as passed into main */
              char **args,              /* args as passed into main */
              int *verb,                /* verbosity setting */
              enum OSTYPE *OS,          /* OS to assume */
              enum ASMDIA *asmb,        /* assembly dialect to assume */
              int *vec,                 /* Vector ISA extension bitfield */
              enum MACHTYPE *mach,     /* machine/arch to assume */
              int *mhz,                /* Clock rate in Mhz */
              int *ptrbits             /* # of bits in ptr: 32/64 */,
              int *nthreads,           /* # of threads */
              char **comps,
              char **gccflags,        /* append flags for user-contrib gcc */
              char **outfile,
              char **srcdir,          /* path to top of source directory */
              char **bindir,          /* path to top of binary directory */
              int *bozol1,            /* Use untuned L1 */
              int *UseArchDef,        /* Use arch defaults */
              int *NoF77,
              char **f2cdefs,         /* F77-to-C interface defines */
              char **ecdefs,          /* extra cpp defines to add to CDEFS */
              char **pmake,           /* parallel make command */
              char **flapack,         /* netlib F77 LAPACK  */
              char **smaflags,       /* single prec muladd flags */
              char **dmaflags,       /* double prec muladd flags */
              char **f77lib,         /* netlib F77 LAPACK  */
              int *flush,             /* size in KB to flush */
              char **targ             /* mach to ssh to*/
             )
{
   int i, k, k0, kn, DoInt;
   char *sp, *sp0;
   char *gcc3=NULL;
   char *cdefs=NULL, *fdefs=NULL;
   char ln[1024];

   *verb = 0;
   *srcdir = *bindir = NULL;
    *bozol1 = 0;
    *UseArchDef = 1;
    *flapack = NULL;
    *f77lib = NULL;
    *smaflags = *dmaflags = NULL;
    *mhz = 0;
   *outfile = NULL;
   *targ = NULL;
   for (k=0; k < NCOMP*3; k++)
      comps[k] = NULL;
   *gccflags = NULL;

   *flush = 0;
   *nthreads = 0;
   *ptrbits = 0;
   *mhz = 0;
   *mach = 0;
   *vec = 0;
   *asmb = 0;
   *OS = 0;
   *verb = 0;
   *NoF77 = 0;
   *nthreads = -1;
   *pmake = NULL;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'f':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *flush = atoi(args[i]);
         break;
      case 't':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *nthreads = atoi(args[i]);
         break;
      case 'b':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *ptrbits = atoi(args[i]);
         break;
      case 'm':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *mhz = atoi(args[i]);
         break;
      case 'A':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *mach = atoi(args[i]);
         break;
      case 'V':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *vec = atoi(args[i]);
         break;
      case 's':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *asmb = atoi(args[i]);
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
      case 'S':
         if (args[i][2] != 'i' && args[i][2] != 's')
            PrintUsage(args[0], i, "-S needs i or s suffix!");
         DoInt = args[i][2] == 'i';
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         sp0 = args[i];
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         if (DoInt)
            k = atoi(args[i]);
         else
            sp = args[i];
         if (!strcmp(sp0, "archdef"))
            *UseArchDef = k;
         else if (!strcmp(sp0, "bozol1"))
            *bozol1 = k;
         else if (!strcmp(sp0, "nof77"))
            *NoF77 = k;
         else if (!strcmp(sp0, "kern"))
            gcc3 = sp;
         else if (!strcmp(sp0, "pmake"))
            *pmake = sp;
        else if (!strcmp(sp0, "flapack"))
           *flapack = sp;
        else if (!strcmp(sp0, "f77lib"))
           *f77lib = sp;
        else if (!strcmp(sp0, "smaflags"))
           *smaflags = sp;
        else if (!strcmp(sp0, "dmaflags"))
           *dmaflags = sp;
         else
            PrintUsage(args[0], i-1, sp0);
         break;
      case 'o':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *outfile = args[i];
         break;
      case 'D':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         if (args[i-1][0] == 'f')
            fdefs = NewAppendedString(fdefs, args[i]);
         else
            cdefs = NewAppendedString(cdefs, args[i]);
         break;
      case 'd':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         sp = args[i-1];
         if (*sp == 's')
            *srcdir = args[i];
         else if (*sp == 'b')
            *bindir = args[i];
         break;
      case 'C':
      case 'F':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         sp = args[i];
         k = -1;
         if (*sp == 'i' && sp[1] == 'c') k = ICC_;
         else if (*sp == 'i' && sp[1] == 'f') k = F77_;
         else if (*sp == 's' && sp[1] == 'k') k = SKC_;
         else if (*sp == 'd' && sp[1] == 'k') k = DKC_;
         else if (*sp == 's' && sp[1] == 'm') k = SMC_;
         else if (*sp == 'd' && sp[1] == 'm') k = DMC_;
         else if (*sp == 'x' && sp[1] == 'c') k = XCC_;
         if (*sp == 'a' && (sp[1] == 'l' || sp[1] == 'c'))
         {  /* only appended flags can be applied to all compilers */
            if (args[i-1][1] == 'F')
            {
               if (args[i-1][2] == 'a')
               {
                  k0 = NCOMP+NCOMP;
                  kn = k0 + NCOMP;
               }
               else
               {
                  k0 = NCOMP;
                  kn = NCOMP+NCOMP;
               }
            }
            else
            {
               k0 = 0;
               kn = NCOMP;
            }
            if (++i >= nargs)
               PrintUsage(args[0], i, "out of arguments");
            for (k=k0; k < kn; k++)
               if (sp[1] == 'l' || k-2*NCOMP != F77_)
                  comps[k] = args[i];
            if (sp[2] == 'g' && args[i-2][1] == 'F')
               *gccflags = args[i];
         }
         else if (*sp == 'g' && sp[1] == 'c')
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, "out of arguments");
            *gccflags = args[i];
         }
         else
         {
            if (k < 0) PrintUsage(args[0], i, args[i]);
            if (args[i-1][1] == 'F')
            {
               k += NCOMP;
               if (args[i-1][2] == 'a')
                  k += NCOMP;
            }
            if (++i >= nargs)
               PrintUsage(args[0], i, "out of arguments");
            comps[k] = args[i];
         }
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
/*
 * allocate these strings ourselves so we can free them later if necessary
 */
   for (i=0; i < 3*NCOMP; i++)
   {
      if (comps[i])
      {
         if (!strcmp(comps[i], "default"))
            comps[i] = NULL;
         else
         {
            sp = malloc(sizeof(char)*(strlen(comps[i])+1));
            strcpy(sp, comps[i]);
            comps[i] = sp;
         }
      }
   }
/*
 * If the special flag -Ss gcc3 is thrown, force gcc3's use for all kernel
 * compilers (standard gcc assumed to be gcc4)
 */
   if (gcc3)
   {
      for (i=0; i < NCOMP; i++)
      {
         if (!comps[i] && (i == SMC_ || i == DMC_ || i == SKC_ || i == DKC_))
            comps[i] = NewStringCopy(gcc3);
      }
   }
   *f2cdefs = fdefs;
   *ecdefs = cdefs;
   if (*ptrbits != 32 && *ptrbits != 64)
      *ptrbits = 0;
}
char *GetPtrbitsFlag(enum OSTYPE OS, enum MACHTYPE arch, int ptrbits,
                     char *comp)
/*
 * RETURNS: string forcing setting of ptrbits for gcc
 */
{
   char *sp = "";
   int i, j, k;

   if (MachIsIA64(arch))
      return(sp);
   if (MachIsMIPS(arch))
      return((ptrbits == 64) ? "-mabi=64" : "-mabi=n32");
   if (!CompIsGcc(comp))
   {
/*
 *    Add correct 64/32 bit flags for Sun workshop compilers
 */
      if (MachIsUS(arch) && CompIsSunWorkshop(comp))
      {
         if (ptrbits == 64)
            sp = (arch == SunUSI || arch == SunUSII) ?
                 "-xarch=v9" : "-xarch=v9b";
         else
            sp = (arch == SunUSI || arch == SunUSII) ?
                 "-xarch=v8plusa" : "-xarch=v8plusb";
      }
      else if (CompIsIBMXL(comp))  /* IBM xl compilers */
         sp = (ptrbits == 64) ? "-q64" : "-q32";
      return(sp);
   }
   GetGccVers(comp, &k, &j, &k, &k);
   if ( !(j >= 3 && (OS != OSOSX || j > 3 || !CompIsAppleGcc(comp))) )
      return(sp);
   else if (OS == OSAIX)
      sp = (ptrbits == 64) ? "-maix64" : "-maix32";
   else if (ptrbits == 64)
     sp = "-m64";
   else if (ptrbits == 32)
     sp = "-m32";
   return(sp);
}

main(int nargs, char **args)
{
   enum OSTYPE OS;
   enum MACHTYPE mach;
   int h, i, j, k, verb, asmb, nof77, mhz;
   int vecexts, ISAX;
   int ptrbits, l2size;
   int delay=0;  /* change this to come from "special" ints in GetFlags */
   int THREADS=0;
   int Use3DNow=0;  /* this needs to come from getflags */
   int ncpu;
   int USEDEFL1, USEARCHDEF;
   char *targ, *sp, *pmake, *flapack;
   char *comps[3*NCOMP], *comp, *flags, *srcdir, *blddir, *f2cdefs, *cdefs;
   char *outfile, *smaflags, *dmaflags, *f77lib, *gccflags, *goodgcc;
   char targarg[256], ln[1024];
   FILE *fpout;
   char *adnames[NARDEF] = {"sKERNDEF", "dKERNDEF", "sMMDEF", "dMMDEF"};

   GetFlags(nargs, args, &verb, &OS, (enum ASMDIA*) &asmb, &vecexts, &mach,
            &mhz, &ptrbits, &ncpu, comps, &gccflags, &outfile, &srcdir,
            &blddir, &USEDEFL1, &USEARCHDEF, &nof77, &f2cdefs, &cdefs, &pmake,
            &flapack, &smaflags, &dmaflags, &f77lib, &l2size, &targ);
   if (ncpu > 1) THREADS = 1;
   if (!outfile)
      fpout = stdout;
   else
      fpout = fopen(outfile, "w");
   assert(fpout);
   assert(srcdir && blddir);
/*
 * Update l2size, and set f2cdefs/cdefs if they are null
 */
   if (!l2size) l2size = 4*1024*1024;
   else l2size *= 1024;
   if (!f2cdefs) f2cdefs = "";
/*
 * Append any appended flags, and then we have just compilers and flags
 */
   for (i=2*NCOMP; i < 3*NCOMP; i++)
   {
      if (comps[i])
      {
         comps[i-NCOMP] = NewAppendedString(comps[i-NCOMP], comps[i]);
         free(comps[i]);
         comps[i] = NULL;
      }
   }
/*
 * If any C compiler is unspecified, use it to specify the others
 * Use DKC preferentially if it is specified
 */
   if (comps[DKC_] && comps[NCOMP+DKC_])
      k = DKC_;
   else
   {
      k = -1;
      for (i=0; i < F77_; i++)
      {
         if (comps[i] && comps[NCOMP+i])
         {
            k = i;
            break;
         }
      }
      if (k < 0)
      {
         fprintf(stderr, "Need a valid C compiler and flags\n");
         exit(100);
      }
   }
   for (i=0; i < F77_; i++)
   {
      if (!comps[i])
         comps[i] = comps[k];
      if (!comps[NCOMP+i])
         comps[NCOMP+i] = comps[NCOMP+k];
   }
/*
 * If F77 compiler unspecified or nof77 asserted, set it to ICC for linking
 */
   if (!comps[F77_] || nof77)
   {
      comps[F77_] = comps[ICC_];
      comps[NCOMP+F77_] = comps[NCOMP+ICC_];
   }
/*
 * Find dominant ISA extension
 */
   ISAX = 0;
   for (i=1; i < NISA && !ISAX; i++)
      if (vecexts & (1<<i))
         ISAX = i;

   fprintf(fpout, "#  ----------------------------\n");
   fprintf(fpout, "#  Make.inc for ATLAS3.8.3\n");
   fprintf(fpout, "#  ----------------------------\n\n");

   fprintf(fpout, "#  ----------------------------------\n");
   fprintf(fpout, "#  Make sure we get the correct shell\n");
   fprintf(fpout, "#  ----------------------------------\n");
   fprintf(fpout, "   SHELL = /bin/sh\n\n");

   fprintf(fpout, "#  -------------------------------------------------\n");
   fprintf(fpout, "#  Name indicating the platform to configure BLAS to\n");
   fprintf(fpout, "#  -------------------------------------------------\n");
   fprintf(fpout, "   ARCH = %s", machnam[mach]);
   fprintf(fpout, "%d", ptrbits);
   if (ISAX)
      fprintf(fpout, "%s", ISAXNAM[ISAX]);
   fprintf(fpout, "\n\n");

   fprintf(fpout, "#  ----------------------------\n");
   fprintf(fpout, "#  Paths to various directories\n");
   fprintf(fpout, "#  ----------------------------\n");
   fprintf(fpout, "   BLDdir = %s\n", blddir);
   fprintf(fpout, "   SRCdir = %s\n", srcdir);
   fprintf(fpout, "   INCAdir = $(BLDdir)/include\n");
   fprintf(fpout, "   INCSdir = $(SRCdir)/include\n");
   fprintf(fpout, "   BINdir = $(BLDdir)/bin\n");
   fprintf(fpout, "   LIBdir = $(BLDdir)/lib\n\n");
   fprintf(fpout, "   SYSdir = $(BLDdir)/tune/sysinfo\n");
   fprintf(fpout, "   GMMdir = $(BLDdir)/src/blas/gemm\n");
   fprintf(fpout, "   GMVdir = $(BLDdir)/src/blas/gemv\n");
   fprintf(fpout, "   GR1dir = $(BLDdir)/src/blas/ger\n");
   fprintf(fpout, "   L1Bdir = $(BLDdir)/src/blas/level1\n");
   fprintf(fpout, "   L2Bdir = $(BLDdir)/src/blas/level2\n");
   fprintf(fpout, "   L3Bdir = $(BLDdir)/src/blas/level3\n");
   fprintf(fpout, "   TSTdir = $(BLDdir)/src/testing\n");
   fprintf(fpout, "   AUXdir = $(BLDdir)/src/auxil\n");
   fprintf(fpout, "   CBLdir = $(BLDdir)/interfaces/blas/C/src\n");
   fprintf(fpout, "   FBLdir = $(BLDdir)/interfaces/blas/F77/src\n");
   fprintf(fpout, "   PTSdir = $(BLDdir)/src/pthreads\n");
   fprintf(fpout, "   MMTdir = $(BLDdir)/tune/blas/gemm\n");
   fprintf(fpout, "   MVTdir = $(BLDdir)/tune/blas/gemv\n");
   fprintf(fpout, "   R1Tdir = $(BLDdir)/tune/blas/ger\n");
   fprintf(fpout, "   L1Tdir = $(BLDdir)/tune/blas/level1\n");
   fprintf(fpout, "   L3Tdir = $(BLDdir)/tune/blas/level3\n");
   fprintf(fpout, "\n");

   fprintf(fpout,
"#  ---------------------------------------------------------------------\n");
   fprintf(fpout,
"#  Name and location of scripts for running executables during tuning\n");
   fprintf(fpout,
"#  ---------------------------------------------------------------------\n");
   fprintf(fpout, "   ATLRUN = $(BLDdir)/bin/ATLrun.sh\n");
   fprintf(fpout, "   ATLFWAIT = $(BLDdir)/bin/xatlas_waitfile\n\n");

   fprintf(fpout, "#  ---------------------\n");
   fprintf(fpout, "#  Libraries to be built\n");
   fprintf(fpout, "#  ---------------------\n");
   fprintf(fpout, "   ATLASlib = $(LIBdir)/libatlas.a\n");
   fprintf(fpout, "   CBLASlib = $(LIBdir)/libcblas.a\n");
   fprintf(fpout, "   F77BLASlib = $(LIBdir)/libf77blas.a\n");
   if (THREADS)
   {
      fprintf(fpout, "   PTCBLASlib = $(LIBdir)/libptcblas.a\n");
      fprintf(fpout, "   PTF77BLASlib = $(LIBdir)/libptf77blas.a\n");
   }
   fprintf(fpout, "   LAPACKlib = $(LIBdir)/liblapack.a\n");
   fprintf(fpout, "   TESTlib = $(LIBdir)/libtstatlas.a\n\n");

   fprintf(fpout, "#  -------------------------------------------\n");
   fprintf(fpout, "#  Upper bound on largest cache size, in bytes\n");
   fprintf(fpout, "#  -------------------------------------------\n");
   fprintf(fpout, "   L2SIZE = -DL2SIZE=%d\n\n", l2size);

   fprintf(fpout, "#  ---------------------------------------\n");
   fprintf(fpout, "#  Command setting up correct include path\n");
   fprintf(fpout, "#  ---------------------------------------\n");
   fprintf(fpout,
           "   INCLUDES = -I$(INCAdir) -I$(INCSdir) -I$(INCSdir)/contrib\n\n");

   fprintf(fpout, "#  -------------------------------------------\n");
   fprintf(fpout, "#  Defines for setting up F77/C interoperation\n");
   fprintf(fpout, "#  -------------------------------------------\n");
   fprintf(fpout, "   F2CDEFS = %s\n\n", f2cdefs);

   fprintf(fpout, "#  ------------------------------\n");
   fprintf(fpout, "#  Architecture identifying flags\n");
   fprintf(fpout, "#  ------------------------------\n");
   fprintf(fpout, "   ARCHDEFS =");
   if (OS != OSOther)
      fprintf(fpout, " -DATL_OS_%s", osnam[OS]);
   if (mach != MACHOther)
      fprintf(fpout, " -DATL_ARCH_%s", machnam[mach]);
   if (mhz)
      fprintf(fpout, " -DATL_CPUMHZ=%d", mhz);
   if (OS == OSSunOS)
      fprintf(fpout, " -DSUN_HR");
   if (OSIsWin(OS))
      fprintf(fpout, " -DGCCWIN -DUseClock");
   for (i=1; i < NISA; i++)
      if (vecexts & (1<<i))
         fprintf(fpout, " -DATL_%s", ISAXNAM[i]);
   if (Use3DNow) fprintf(fpout, " -DATL_3DNowFLOPS");
   if (ptrbits == 64)
      fprintf(fpout, " -DATL_USE64BITS");
   if (mach == IA64Itan || mach == IA64Itan2 )
      fprintf(fpout, " -DATL_MAXNREG=128");
   if (asmb != ASM_None) fprintf(fpout, " -DATL_%s", ASMNAM[asmb]);
   if (mach == IA64Itan2)
      fprintf(fpout, " -DATL_IntelIccBugs");
/*
 * Need up update handling of apple vs. gnu gcc for altivec
 */
#if 0
   if (ISAX == ISA_AV && strstr(CC, "gcc") && strstr(MCC, "gcc") &&
       !GccIsApple(CC))
      fprintf(fpout, " -DATL_%s -DATL_AVgcc", ISAXNAM[ISAX]);
#endif
   fprintf(fpout, "\n\n");

   fprintf(fpout,
   "#  -------------------------------------------------------------------\n");
   fprintf(fpout,
   "#  NM is the flag required to name a compiled object/executable\n");
   fprintf(fpout,
   "#  OJ is the flag required to compile to object rather than executable\n");
   fprintf(fpout, "#  These flags are used by all compilers.\n");
   fprintf(fpout,
   "#  -------------------------------------------------------------------\n");
   fprintf(fpout, "   NM = -o\n");
   fprintf(fpout, "   OJ = -c\n\n");

   sprintf(ln, "%s/CONFIG/src/CompMake.txt", srcdir);
   DisplayFile(ln, fpout, 0);
   fprintf(fpout, "   CDEFS = $(L2SIZE) $(INCLUDES) $(F2CDEFS) $(ARCHDEFS)");
   if (cdefs) fprintf(fpout, " %s", cdefs);
   if (THREADS)
   {
      fprintf(fpout, " -DATL_NCPU=%d", ncpu);
      if (OS == OSFreeBSD) fprintf(fpout, " -D_THREAD_SAFE -D_REENTRANT");
      if (OS == OSAIX) fprintf(fpout, " -DIBM_PT_ERROR");
      if (OS == OSIRIX) fprintf(fpout, " -D_POSIX_C_SOURCE=199506L");
   }
   if (delay) fprintf(fpout, " -DATL_FOPENDELAY");
   fprintf(fpout, "\n\n");
   for (i=0; i < NCOMP; i++)
   {
      fprintf(fpout, "   %s = %s\n", COMPNAME[i], comps[i]);
      if (i == F77_)
         fprintf(fpout, "   %sFLAGS = %s\n", COMPNAME[i], comps[NCOMP+i]);
      else if (i == ICC_ || i == XCC_)
         fprintf(fpout, "   %sFLAGS = $(CDEFS) %s\n",
                 COMPNAME[i], comps[NCOMP+i]);
      else /* non-interf comps don't include CDEFS by default (added later) */
         fprintf(fpout, "   %sFLAGS = %s\n",
                 COMPNAME[i], comps[NCOMP+i]);
   }
   fprintf(fpout, "   SMAFLAGS =");
   if (smaflags)
      fprintf(fpout, " %s", smaflags);
   fprintf(fpout, "\n   DMAFLAGS =");
   if (dmaflags)
      fprintf(fpout, " %s", dmaflags);
   fprintf(fpout, "\n");
   fprintf(fpout, "   CKC = $(SKC)\n");
   fprintf(fpout, "   ZKC = $(DKC)\n");
   fprintf(fpout, "   sKCFLAGS = $(CDEFS) $(SKCFLAGS)\n");
   fprintf(fpout, "   dKCFLAGS = $(CDEFS) $(DKCFLAGS)\n");
   fprintf(fpout, "   cKCFLAGS = $(CDEFS) $(SKCFLAGS)\n");
   fprintf(fpout, "   zKCFLAGS = $(CDEFS) $(DKCFLAGS)\n");

   for (i=0; i < NCOMP; i++)
   {
      if (i == XCC_) continue;  /* do not accept cross-compiler */
      j = strlen(comps[i]);
      if (j >= 3 && comps[i][j-3] == 'g' &&
          comps[i][j-2] == 'c' && comps[i][j-1] == 'c')
         break;
   }
   goodgcc = (j < NCOMP) ? comps[j] : "gcc";
   fprintf(fpout, "   GOODGCC = %s", goodgcc);
   if (gccflags)
      fprintf(fpout, " %s", gccflags);
   GetGccVers(goodgcc, &i, &j, &k, &k);
   if (OS == OSOSX && j > 3)  /* need bullshit apple annoyance flag */
      fprintf(fpout, " -force_cpusubtype_ALL");
   sp = GetPtrbitsFlag(OS, mach, ptrbits, goodgcc);
   if (strlen(sp) > 0)
       fprintf(fpout, " %s", sp);
   fprintf(fpout, "\n");

   fprintf(fpout, "   LDFLAGS =");
   if (MachIsX86(mach))
   {
      if (ptrbits == 32)
         fprintf(fpout, " -melf_i386");
      else if (ptrbits == 64)
         fprintf(fpout, " -melf_x86_64");
   }
   fprintf(fpout, "\n   F77SYSLIB = %s\n", f77lib ? f77lib : "");
   fprintf(fpout, "   BC = $(ICC)\n");
   fprintf(fpout, "   NCFLAGS = $(ICCFLAGS)\n");
   fprintf(fpout, "\n   CLINKER = $(ICC)\n   CLINKFLAGS = $(ICCFLAGS)\n");
   fprintf(fpout, "   FLINKER = $(F77)\n   FLINKFLAGS = $(F77FLAGS)\n");
   fprintf(fpout, "   FCLINKFLAGS = $(FLINKFLAGS)");
   if (strstr(comps[F77_], "ifort") && !OSIsWin(OS))
      fprintf(fpout, " -nofor_main");
   fprintf(fpout, "\n   ARCHIVER = ar\n");
   fprintf(fpout, "   ARFLAGS  = r\n");
/*
 * JF Mertens says that even x86 OS X still need ranlib for safety
 */
   if (OS == OSOSX /* && MachIsPPC(mach) */)
      fprintf(fpout, "   RANLIB   = ranlib\n\n");
   else
      fprintf(fpout, "   RANLIB   = echo\n\n");

   fprintf(fpout, "#  -------------------------------------\n");
   fprintf(fpout, "#  tar, gzip, gunzip, and parallel make\n");
   fprintf(fpout, "#  -------------------------------------\n");
   fprintf(fpout, "   TAR    = tar\n");
   fprintf(fpout, "   GZIP   = gzip\n");
   fprintf(fpout, "   GUNZIP = gunzip\n");
   fprintf(fpout, "   PMAKE  = %s\n\n", pmake ? pmake : "$(MAKE)");
/*
 * Need to add libs to GetFlags and update GetSysLib to do this right
*/
   fprintf(fpout, "#  ------------------------------------\n");
   fprintf(fpout, "#  Reference and system libraries\n");
   fprintf(fpout, "#  ------------------------------------\n");
   fprintf(fpout, "   BLASlib = \n");
   fprintf(fpout, "   FBLASlib = \n");
   fprintf(fpout, "   FLAPACKlib = ");
   if (flapack) fprintf(fpout, "%s", flapack);
   fprintf(fpout, "\n");
   fprintf(fpout, "   SLAPACKlib = \n");
   if (THREADS)
      fprintf(fpout, "   LIBS = -lpthread");
   else
      fprintf(fpout, "   LIBS =");
   if (!OSIsWin(OS))
      fprintf(fpout, " -lm");
   fprintf(fpout, "\n\n");

   fprintf(fpout,
   "#  --------------------------------------------------------------------\n");
   fprintf(fpout,
   "#  Compiler names for architectural defaults and flags to atlas_install\n");
   fprintf(fpout,
   "#  --------------------------------------------------------------------\n");
   for (i=0; i < NCOMP; i++)
   {
      sp = NewStringCopy(COMPNAME[i]);
      for (j=0; sp[j]; j++)
        sp[j] = tolower(sp[j]);
      fprintf(fpout, "   %sD = ", sp);
      free(sp);
/*
 *    Regardless of actual names, use standard gnu compiler names for defs
 */
      if (CompIsGcc(comps[i]))
      {
         if (i == F77_)
         {
            GetGccVers(comps[i], &k, &j, &k, &k);
            if (j < 4)
               sp = NewStringCopy("g77");
            else
               sp = NewStringCopy("gfortran");
         }
         else sp = NewStringCopy("gcc");
      }
      else
      {
         sp = NameWithoutPath(comps[i]);
         if (!strncmp(sp, "ATLwin_", 7))
            sp = NewStringCopy(comps[i]+7);
      }
      fprintf(fpout, "%s\n", sp);
      free(sp);
   }
   fprintf(fpout, "   INSTFLAGS = -1 %d -a %d\n\n", USEDEFL1, USEARCHDEF);

   fprintf(fpout, "#  ---------------------------------------\n");
   fprintf(fpout, "#  Generic targets needed by all makefiles\n");
   fprintf(fpout, "#  ---------------------------------------\n");
   fprintf(fpout, "do_it: all\n");
   if (delay)
   {
      fprintf(fpout, "   waitfile = wfdefault\n");
      fprintf(fpout, "waitfile:\n\tcd $(BINdir) ; make xatlas_waitfile\n");
      fprintf(fpout, "\t$(ATLFWAIT) -s %d -f $(waitfile)\n", delay);
   }
   else fprintf(fpout, "waitfile:\n");
   if (fpout != stdout && fpout != stderr) fclose(fpout);
   exit(0);
}
