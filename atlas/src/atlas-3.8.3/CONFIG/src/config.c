#include "atlconf.h"

int GetStrProbe(int verb, char *targarg, char *prb, char *id, char *str)
/*
 * Performs probe where output is string delimited by '', leaving string
 * in str.
 * RETURNS: 0 on success, non-zero on error
 */
{
   char *sp;
   char ln[2048], ln2[2048];
   int ierr=0, i;

   ln[0] = ln2[0] = '\0';
   sprintf(ln, "make IRun_%s args=\"-v %d %s\" | fgrep '%s='",
           prb, verb, targarg, id);
   if (verb > 1)
      printf("cmnd=%s\n", ln);
   ierr = CmndOneLine(NULL, ln, ln2);
   if (!ierr)
   {
      sprintf(ln, "%s='", id);
      sp = strstr(ln2, ln);
      if (sp)
      {
         for (sp += strlen(ln); *sp && *sp != '\''; sp++)
            *str++ = *sp;
         if (!(*sp)) ierr = -2;
         *str = '\0';
      }
      else ierr = -4;
   }
   return(ierr);
}

int GetIntProbe(int verb, char *targarg, char *prb, char *id, int N)
{
   char ln[2048], ln2[2048];
   int iret=0, ierr;

   ln[0] = ln2[0] = '\0';
   sprintf(ln, "make IRun_%s args=\"-v %d %s\" | fgrep '%s='",
           prb, verb, targarg, id);
   if (verb > 1)
      printf("cmnd=%s\n", ln);
   ierr = CmndOneLine(NULL, ln, ln2);
   if (!ierr)
      iret = GetLastInt(ln2);
   if (N)
   {
      if (iret > N || iret < 1)
      {
         printf("\nBad %s value=%d, ierr=%d, ln2='%s'\n", id, iret, ierr, ln2);
         iret = 0;
      }
   }
   return(iret);
}

int GetIntProbeSure(int verb, char *targarg, char *prb, char *id,
                    int *sure)
{
   char ln[2048], ln2[2048];
   int iret=0, ierr;

   ln[0] = ln2[0] = '\0';
   sprintf(ln, "make IRun_%s args=\"-v %d %s\" | fgrep '%s='",
           prb, verb, targarg, id);
   if (verb > 1)
      printf("cmnd=%s\n", ln);
   ierr = CmndOneLine(NULL, ln, ln2);
   if (!ierr)
   {
      iret = GetFirstInt(ln2);
      *sure = GetLastInt(ln2);
   }
   return(iret);
}

void TransCompsToFlags(char **comps, char *flags)
{
   int i, j;
/*
 * WARNING: if you change the order of the ICC_, etc, must change
 *          compnames!
 */
   char *compnames[7] = {"ic", "sm", "dm", "sk",  "dk", "xc", "if"};

   *flags = '\0';
/*
 * Pass any override compilers/flags and appends to probe
 */
   for (j=i=0; i < NCOMP; i++)
   {
      if (comps[i])
         j += sprintf(flags+j, "-C %s '%s' ", compnames[i], comps[i]);
      if (comps[NCOMP+i])
         j += sprintf(flags+j, "-F %s '%s' ", compnames[i], comps[NCOMP+i]);
      if (comps[2*NCOMP+i])
         j += sprintf(flags+j, "-Fa %s '%s' ", compnames[i], comps[2*NCOMP+i]);
   }
}

char *ProbeComp(int verb, char *targarg, enum OSTYPE OS, enum MACHTYPE arch,
                 char **comps, int nof77, int nocygwin, int ptrbits)
/*
 * RETURNS: f2c define string
 */
{
   char ln[4096], flags[2048], comp[512];
   int f2cname, f2cint, f2cstr;
   char *f2cdefs;
   int i, if77=0;
   FILE *fpin;

   TransCompsToFlags(comps, flags);
   i = sprintf(ln,
   "make IRun_comp args=\"-v %d -o atlconf.txt -O %d -A %d -Si nof77 %d %s %s",
               verb, OS, arch, nof77, targarg, flags);
   if (ptrbits == 64 || ptrbits == 32)
      i += sprintf(ln+i, " -b %d", ptrbits);
   sprintf(ln+i, "\"");
   if (verb > 1)
      fprintf(stderr, "cmnd='%s'\n", ln);
   assert(!system(ln));
   fpin = fopen("atlconf.txt", "r");
   assert(fpin);
   while (fgets(ln, 4096, fpin))
   {
      if (ln[0] != '#')
      {
         if (isdigit(ln[0]))
         {
            assert(sscanf(ln, "%d '%[^']' '%[^']", &i, comp, flags) == 3);
            assert(i >= 0 && i < NCOMP);
            comps[i] = NewStringCopy(comp);
            comps[NCOMP+i] = NewStringCopy(flags);
         }
         else
         {
            for (i=0; ln[i] && ln[i] != '('; i++);
            assert(ln[i] = '(');
            assert(sscanf(ln+i+1, "%d,%d,%d", &f2cname, &f2cint, &f2cstr) == 3);
            if77 = 1;
         }
      }
   }
   fclose(fpin);
/*
 * For windows only, all non-XCC gnu compilers get -mnocygwin flag if user has
 * asked for it
 */
   if (nocygwin && OSIsWin(OS))
   {
      for (i=0; i < NCOMP; i++)
      {
         if (i != XCC_)
         {
            if (CompIsGcc(comps[i]))
               comps[i+NCOMP] = NewAppendedString(comps[i+NCOMP],"-mno-cygwin");
         }
      }
   }
/*
 * Only f77 and xcc allowed to be NULL
 */
   for (i=0; i < NCOMP; i++)
   {
      if (comps[i] == NULL)
      {
         if (i == XCC_)
         {
            comps[XCC_] = NewStringCopy(comps[ICC_]);
            comps[NCOMP+XCC_] = NewStringCopy(comps[ICC_+NCOMP]);
         }
         else assert(i == F77_);
      }
      if (comps[i+NCOMP] == NULL)
      {
         if (i == XCC_)
            comps[XCC_+NCOMP] = NewStringCopy(comps[ICC_+NCOMP]);
         else assert(i == F77_);
      }
   }
/*
 * Echo compiler info to screen if user has asked for verbose output
 */
   if (verb)
   {
      fprintf(stdout, "Selected compilers:\n");
      for (i=0; i < NCOMP; i++)
      {
         fprintf(stdout, "%s = '%s' '%s'\n", COMPNAME[i],
                 comps[i] ? comps[i]:"none",
                 comps[NCOMP+i]?comps[NCOMP+i]:"none");
      }
      if (if77)
      {
         fprintf(stdout, "\nF77 calling C interface information:\n");
         fprintf(stdout, "   Name decoration = %s\n", f2c_namestr[f2cname]);
         fprintf(stdout, "   Integer style   = %s\n", f2c_intstr[f2cint]);
         fprintf(stdout, "   String style    = %s\n", f2c_strstr[f2cstr]);
      }
      else
         fprintf(stderr, "F77/C interface not defined!");
   }
   if (!if77)
      f2cdefs = NULL;
   else
   {
      sprintf(ln, "-D%s -D%s -DString%s",
              f2c_namestr[f2cname], f2c_intstr[f2cint], f2c_strstr[f2cstr]);
      f2cdefs = malloc(sizeof(char)*(strlen(ln)+1));
      assert(f2cdefs);
      strcpy(f2cdefs, ln);
   }
   return(f2cdefs);
}

int PathLength(char *str)
/*
 * Given string, finds end unix path, allowing \ to mean sticky space, stops at
 * first non-sticky pass.  Skips any leading spaces.
 * RETURNS: index of first non-sticky space
 * NOTE: Assumes str[-1] valid
 */
{
   int i;
   for (i=0; isspace(str[i]); i++);  /* skip leading spaces */
   for (; str[i]; i++)
      if (str[i] == ' ' && str[i-1] != '\\') break;
   return(i);
}

char *ProbeF77LIB(int verb, char *targarg, enum OSTYPE OS,
                  enum MACHTYPE arch, char **comps, char *f2cdefs, int nof77)
/*
 * Tries to find the directory that needs to be included to link in f77
 * routines
 */
{
   char cmnd[1024], res[2048];
   char *f77lib, *F77LIBdir, *sp;
   int i;
   char ch;

   if (nof77 || !f2cdefs) return(NULL);
   if  (strstr(comps[F77_], "g77")) f77lib = "g2c";
   else if (strstr(comps[F77_], "gfortran")) f77lib = "gfortran";
   else
   {
      if (verb > 1)
         fprintf(stderr, "Unknown F77 compiler, leaving F77LIBS blank!\n");
      return(NULL);
   }
   sprintf(cmnd,  "make IRunFlib F77='%s' F77FLAGS='%s'",
           comps[F77_], comps[NCOMP+F77_]);
   if (verb > 1)
      fprintf(stderr, "LIBDIR cmnd = %s\n", cmnd);
   CmndOneLine(NULL, cmnd, res);
   if (verb > 1)
      fprintf(stderr, "LIBDIR res = %s\n", res);
/*
 * Find which -L leads us to f77lib
 */
   sp = res;
   while (sp = strstr(sp, "-L"))
   {
      sp += 2;
      i = PathLength(sp);
      if (i)
      {
         ch = sp[i];
         sp[i] = '\0';
         sprintf(cmnd, "make IRunTestCFLink F77='%s' F77FLAGS='%s' CC='%s' CCFLAGS='%s -L%s' F77LIB='%s' LIBS='-lm'",
                 comps[F77_], comps[NCOMP+F77_],
                 comps[ICC_], comps[NCOMP+ICC_], sp, f77lib);
         if (verb > 1)
            fprintf(stderr, "cmnd = %s\n", cmnd);
         if (verb > 1)
            fprintf(stderr, "Trying F77 link path of %s ... ", sp);
         if (!system(cmnd))
         {
            if (verb > 1) fprintf(stderr, "SUCCESS!\n");
            F77LIBdir = NewStringCopy(sp-2);
            F77LIBdir = NewAppendedString(F77LIBdir, "-l");
            F77LIBdir = NewAppendedString(F77LIBdir, f77lib);
            if (OSIsWin(OS))
               F77LIBdir = NewAppendedString(F77LIBdir, "-lgcc");
            if (verb)
               fprintf(stderr, "F77LIB = %s\n", F77LIBdir);
            return(F77LIBdir);
         }
         if (verb > 1) fprintf(stderr, "REJECTED!!!!\n");
         sp[i] = ch;
         sp += i;
      }
   }
   if (verb) fprintf(stderr, "F77LIB left blank");
   return(NULL);
}

enum OSTYPE ProbeOS(int verb, char *targarg)
{
   enum OSTYPE iret;

   iret = GetIntProbe(verb, targarg, "OS", "OS", NOS);
   printf("\nOS configured as %s (%d)\n", osnam[iret], iret);
   assert(iret);
   return(iret);
}

enum ASMDIA ProbeAsm(int verb, char *targarg, enum OSTYPE OS)
{
   enum ASMDIA asmd=ASM_None;
   char ln[1024];

   sprintf(ln, "%s -O %d", targarg, OS);
   asmd = GetIntProbe(verb, ln, "asm", "ASM", NASMD);
   printf("\nAssembly configured as %s (%d)\n", ASMNAM[asmd], asmd);
   return(asmd);
}

int ProbeVecs(int verb, char *targarg, enum OSTYPE OS, enum ASMDIA asmb)
{
   int i, iret;
   char ln[1024];

   sprintf(ln, "%s -O %d -s %d", targarg, OS, asmb);
   iret = GetIntProbe(verb, ln, "vec", "VECFLAG", (1<<NISA));
   for (i=0; i < NISA && (iret & (1<<i)) == 0; i++);
   if (i == NISA)
      i = 0;
   printf("\nVector ISA Extension configured as  %s (%d,%d)\n",
          ISAXNAM[i], i, iret);
   return(iret);
}

int ProbeArch(int verb, char *targarg, enum OSTYPE OS, enum ASMDIA asmb)
{
   int i, iret;
   char ln[1024];

   sprintf(ln, "%s -O %d -s %d -a", targarg, OS, asmb);
   iret = GetIntProbe(verb, ln, "arch", "MACHTYPE", NMACH);
   printf("\nArchitecture configured as  %s (%d)\n",
          machnam[iret], iret);
   return(iret);
}

int ProbeMhz(int verb, char *targarg, enum OSTYPE OS, enum ASMDIA asmb)
{
   int i, iret;
   char ln[1024];

   sprintf(ln, "%s -O %d -s %d -m", targarg, OS, asmb);
   iret = GetIntProbe(verb, ln, "arch", "CPU MHZ", 16384);
   printf("\nClock rate configured as %dMhz\n", iret);
   return(iret);
}

int ProbeNcpu(int verb, char *targarg, enum OSTYPE OS, enum ASMDIA asmb)
{
   int i, iret;
   char ln[1024];

   sprintf(ln, "%s -O %d -s %d -n", targarg, OS, asmb);
   iret = GetIntProbe(verb, ln, "arch", "NCPU", 2048);
   printf("\nMaximum number of threads configured as  %d\n", iret);
   return(iret);
}

int ProbePtrbits(int verb, char *targarg, enum OSTYPE OS, enum ASMDIA asmb)
{
   int i, iret;
   char ln[1024];

   sprintf(ln, "%s -O %d -s %d -b", targarg, OS, asmb);
   iret = GetIntProbeSure(verb, ln, "arch", "PTR BITS", &i);
/*
 * If it's not 64-bit, make sure it's not just because of flag setting
 */
   if (iret != 64)
   {
   }
   if (iret != 64)
      iret = 32;
   printf("\nPointer width configured as %d\n", iret);
   return(iret);
}

int ProbeCPUThrottle(int verb, char *targarg, enum OSTYPE OS, enum ASMDIA asmb)
{
   int i, iret;
   char ln[1024];
   sprintf(ln, "%s -O %d -s %d -t", targarg, OS, asmb);
   iret = GetIntProbe(verb, ln, "arch", "CPU THROTTLE", 0);
   if (iret) printf("CPU Throttling apparently enabled!\n");
   else printf("Cannot detect CPU throttling.\n");
   return(iret);
}


void Comps2Flags(char **comps, char *ln)
/*
 * Takes the comps array (1st NCOMP entries are compilers, next NCOMP entries
 * flags for those compilers) and translates them into the corresponding flags
 * for xspew (or indeed this config.c)
 * NOTE: assumes appended flags (2*NCOMP+i) already appended to flags (NCOMP+i)
 */
{
   char *cname[NCOMP] = {"ic", "sm", "dm", "sk", "dk", "xc", "if"};
   int i, j=0;
   ln[0] = '\0';
   for (i=0; i < NCOMP; i++)
   {
      if (comps[i])
         j += sprintf(ln+j, "-C %s '%s' ", cname[i], comps[i]);
      if (comps[NCOMP+i])
         j += sprintf(ln+j, "-F %s '%s' ", cname[i], comps[NCOMP+i]);
   }
}

char *ProbePmake(int verb, enum OSTYPE OS, int ncpu)
/*
 * WARNING: if cross-comp really worked, this would be ncpu of front-end,
 *          not backend!
 */
{
   char args[256], res[1024];
   char *sp;

   sprintf(args, "-O %d -t %d", OS, ncpu);
   if (!GetStrProbe(verb, args, "pmake", "PMAKE", res))
   {
      sp = malloc((strlen(res)+1)*sizeof(char));
      strcpy(sp, res);
   }
   else sp = NULL;
   if (sp) printf("Parallel make command configured as '%s'\n", sp);
   else printf("Parallel make not configured.\n");
   return(sp);
}

void SpewItForth(int verb, enum OSTYPE OS, enum MACHTYPE arch, int mhz,
                 enum ASMDIA asmb, int vecexts, int ptrbits, int ncpu,
                 int l2size, char *srcdir, char *bindir, int bozol1,
                 int archdef, int nof77, char **comps, char *gccflags,
                 char *f2cdefs, char *cdefs, char *pmake, char *flapack,
                 char *smaflags, char *dmaflags, char *f77libs)
/*
 * Calls xspew with correct arguments to build required Make.inc
 */
{
   char ln[4096], compsflags[1024], archflags[1024];
   int i;
   assert(CmndResults(NULL, "make xspew"));
/*
 * Translate compiler/flag array to xspew flags
 */
   Comps2Flags(comps, compsflags);
   i = sprintf(ln, "./xspew -v %d -O %d -A %d -m %d -s %d -V %d -b %d -t %d -f %d -d s '%s' -d b '%s' -D c '%s' -D f '%s' %s -Si archdef %d -Si bozol1 %d -Si nof77 %d -o Make.inc",
               verb, OS, arch, mhz, asmb, vecexts, ptrbits, ncpu, l2size,
               srcdir, bindir, cdefs ? cdefs:"", f2cdefs? f2cdefs : "",
               compsflags, archdef, bozol1, nof77);
   if (pmake)
      i += sprintf(ln+i, " -Ss pmake '%s'", pmake);
   if (flapack)
      i += sprintf(ln+i, " -Ss flapack '%s'", flapack);
   if (smaflags)
      i += sprintf(ln+i, " -Ss smaflags '%s'", smaflags);
   if (dmaflags)
      i += sprintf(ln+i, " -Ss dmaflags '%s'", dmaflags);
   if (f77libs)
      i += sprintf(ln+i, " -Ss f77lib '%s'", f77libs);
   if (gccflags)
      i += sprintf(ln+i, " -Fa gc '%s'", gccflags);
   if (verb > 1)
      fprintf(stderr, "cmnd='%s'\n", ln);
   assert(!system(ln));
}

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
        "      -Si nocygwin <0/1> : Do/don't depend on GPL cygwin library\n");
      fprintf(stderr,
        "                           (Windows compiler/cygwin install only)\n");
      fprintf(stderr,
        "      -Si cputhrchk <0/1> : Ignore/heed CPU throttle probe\n");
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
              int *NoCygwin,
              int *ThrChk,
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
   *NoCygwin = 0;
   *NoF77 = 0;
   *ThrChk = 1;
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
         else if (!strcmp(sp0, "nocygwin"))
            *NoCygwin = k;
         else if (!strcmp(sp0, "cputhrchk"))
            *ThrChk = k;
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

main(int nargs, char **args)
{
   enum OSTYPE OS;
   enum MACHTYPE mach;
   int i, verb, asmb, f2cname, f2cint, f2cstr, ncpu, nof77, nocygwin;
   int thrchk, mhz;
   int j, k, h, vecexts;
   int ptrbits, l2size, bozol1, archdef;
   char *targ, *f2cdefs, *cdefs, *srcdir, *bindir, *outfile, *sp;
   char targarg[256];
   char *comps[3*NCOMP], *gccflags;
   char *pmake, *flapack, *smaflags, *dmaflags, *f77libs;

   GetFlags(nargs, args, &verb, &OS, (enum ASMDIA*) &asmb, &vecexts, &mach,
            &mhz, &ptrbits, &ncpu, comps, &gccflags, &outfile, &srcdir, &bindir,
            &bozol1, &archdef, &nof77, &nocygwin, &thrchk, &f2cdefs, &cdefs,
            &pmake, &flapack, &smaflags, &dmaflags, &f77libs, &l2size, &targ);
   if (targ)
      sprintf(targarg, "-T %s", targ);
   else
      targarg[0] = '\0';
   if (OS == OSOther)
      OS = ProbeOS(verb, targarg);
   if (asmb == ASM_None)
      asmb = ProbeAsm(verb, targarg, OS);
   else if (asmb < 0)
      asmb = 0;
   if (!vecexts)
      vecexts = ProbeVecs(verb, targarg, OS, asmb);
   else if (vecexts < 0)
      vecexts = 0;
   if (mach == MACHOther)
      mach = ProbeArch(verb, targarg, OS, asmb);
   if (!mhz)
      mhz = ProbeMhz(verb, targarg, OS, asmb);
   if (ncpu < 0)
      ncpu = ProbeNcpu(verb, targarg, OS, asmb);
   if (!pmake && ncpu > 1)
      pmake = ProbePmake(verb, OS, ncpu);
   if (ptrbits == 0)
   {
      if (asmb == gas_x86_64)
         ptrbits = 64;
      else
         ptrbits = ProbePtrbits(verb, targarg, OS, asmb);
   }
   if (ProbeCPUThrottle(verb, targarg, OS, asmb))
   {
      fprintf(stderr,
         "It appears you have cpu throttling enabled, which makes timings\n");
      fprintf(stderr,
              "unreliable and an ATLAS install nonsensical.  Aborting.\n");
      fprintf(stderr,
              "See ATLAS/INSTALL.txt for further information\n");
      if (thrchk) exit(1);
      else fprintf(stderr, "Ignoring CPU throttling by user override!\n\n");
   }
/*
 * Override 32/64 bit assembler if asked
 */
   if (asmb == gas_x86_64 && ptrbits == 32)
      asmb = gas_x86_32;
   else if (asmb == gas_x86_32 && ptrbits == 64)
      asmb = gas_x86_64;

   sp = ProbeComp(verb, targarg, OS, mach, comps, nof77, nocygwin, ptrbits);
   if (nof77)
      f2cdefs = "-DATL_NoF77";
   else if (!f2cdefs) f2cdefs = sp;
   if (!f77libs)
      f77libs = ProbeF77LIB(verb, targarg, OS, mach, comps, f2cdefs, nof77);
/*
 * If user has not specified muladd flags (which are suffixed to kernel flags),
 * add flags to keep gcc 4 from hanging, if necessary
 */
   if (!smaflags)
   {
      if (CompIsGcc(comps[SKC_]))
      {
         GetGccVers(comps[SKC_], &h, &i, &j, &k);
         if (i >= 4)
            smaflags = "-fno-tree-loop-optimize";
         else if (i == 3) /* reduce opt to avoid compiler hang */
            smaflags = "-O0";
      }
      else if (OS == OSIRIX && CompIsMIPSpro(comps[SKC_]))
         smaflags = " -O2";
   }
   if (!dmaflags)
   {
      if (CompIsGcc(comps[DKC_]))
      {
         GetGccVers(comps[DKC_], &h, &i, &j, &k);
         if (i >= 4)
            dmaflags = "-fno-tree-loop-optimize";
         else if (i == 3) /* reduce opt to avoid compiler hang */
            dmaflags = "-O0";
      }
      else if (OS == OSIRIX && CompIsMIPSpro(comps[DKC_]))
         dmaflags = " -O2";
   }
   SpewItForth(verb, OS, mach, mhz, asmb, vecexts, ptrbits, ncpu, l2size,
               srcdir, bindir, bozol1, archdef, nof77, comps, gccflags,
               f2cdefs, cdefs, pmake, flapack, smaflags, dmaflags, f77libs);
/*
 * Cleanup directory, and exit
 */
   system("make confclean");
   exit(0);
}
