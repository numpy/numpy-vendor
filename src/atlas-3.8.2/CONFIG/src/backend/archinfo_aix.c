#include "atlconf.h"

void PrintUsage(char *name, int i)
{
   fprintf(stderr, "USAGE: %s -v (verb) -b (@ bits) -a (arch) -n (ncpu) -c <ncache> -C <lvl> (cache size) -m (Mhz) -t (cpu throttling)\n", name);
   exit(i);
}

int GetFlags(int nargs, char **args, int *CacheLevel)
{
   int i, flag = 0;

   *CacheLevel = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0], i);
      switch(args[i][1])
      {
      case 'n':
         flag |= Pncpu;
         break;
      case 'c':
         flag |= Pncache;
         break;
      case 'C':
         if (++i > nargs)
            PrintUsage(args[0], i);
         *CacheLevel = atoi(args[i]);
         break;
      case 'v':
         flag |= Pverb;
         break;
      case 'm':
         flag |= PMhz;
         break;
      case 'a':
         flag |= Parch;
         break;
      case 'b':
         flag |= P64;
         break;
      case 't':
         flag |= Pthrottle;
         break;
      default:
         PrintUsage(args[0], i);
      }
   }
   if (!flag)
     flag = Parch | P64;
   return(flag);
}

enum MACHTYPE ProbeArch()
{
   enum ARCHFAM fam;
   enum MACHTYPE mach=MACHOther;
   int ierr, i;
   char res[1024];

   fam = ProbeArchFam(NULL);
   switch(fam)
   {
   case AFPPC: /* don't know */
      if (!CmndOneLine(NULL, "/usr/sbin/prtconf | fgrep 'Processor Type'", res))
      {
         if (strstr(res, "PowerPC_POWER5"))
            mach = IbmPwr5;
         else if (strstr(res, "PowerPC_POWER4"))
            mach = IbmPwr4;
      }
      break;
   }
   return(mach);
}

int ProbeNCPU()
{
   int ncpu = 0;
   char *reslns, res[1024];

   if (!CmndOneLine(NULL, "/usr/sbin/prtconf | fgrep 'Number Of Processors'",
                    res))
      ncpu = GetLastInt(res);
   return(ncpu);
}

int ProbePointerBits(int *sure)
{
   int iret = 32;
   char *uname;
   char cmnd[1024], res[1024];

/*
 * Note this is a weak probe, archinfo_x86 much better . . .
 */
   *sure = 0;
   uname = FindUname(NULL);
   sprintf(cmnd, "%s -a", uname);
/*
 * This probe should be running on backend; if its ptr length is 8, we've
 * definitely got a 64 bit machine
 * NOTE: getting 4 could be a result of compiler flags on a 64-bit arch,
 * so reverse is not dispositive
 */
   if (sizeof(void*) == 8)
   {
      iret = 64;
      *sure = 1;
   }
   else if (!CmndOneLine(NULL, "/usr/sbin/prtconf | fgrep 'Kernel Type'", res))
   {
      if (GetLastInt(res) == 64)
      {
         iret = 64;
         *sure = 1;
      }
   }
   else if (!CmndOneLine(NULL, cmnd, res))
   {
/*
 *    If uname is a known 64-bit platform, we're sure we've got OS support
 *    for 64bits (may not have compiler support, but that's not our fault)
 */
      if (strstr(res, "x86_64") || strstr(res, "ppc64") || strstr(res, "ia64"))
      {
         iret = 64;
         *sure = 1;
      }
   }
   return(iret);
}

int ProbeMhz()
{
   int mhz=0;
   char res[1024];
   if (!CmndOneLine(NULL, "/usr/sbin/prtconf | fgrep 'Processor Clock Speed'",
                    res))
      mhz = GetLastInt(res);  /* assumes clock speed always given in MHz */
   return(mhz);
}

int ProbeThrottle()
/*
 * RETURNS: 1 if cpu throttling is detected, 0 otherwise
 */
{
   int iret=0;
/*
 * I have no idea how to do this for SunOS/Irix/AIX
 */
   return(iret);
}

main(int nargs, char **args)
{
   int flags, CacheLevel, ncpu, mhz, bits, sure;
   enum MACHTYPE arch=MACHOther;

   flags = GetFlags(nargs, args, &CacheLevel);
   if (flags & Parch)
   {
      arch = ProbeArch();
      if (flags & Pverb)
         printf("Architecture detected as %s.\n", machnam[arch]);
      printf("MACHTYPE=%d\n", arch);
   }
   if (flags & Pncpu)
      printf("NCPU=%d\n", ProbeNCPU());
   if (flags & PMhz)
      printf("CPU MHZ=%d\n", ProbeMhz());
   if (flags & Pthrottle)
      printf("CPU THROTTLE=%d\n", ProbeThrottle());
   if (flags & P64)
   {
      bits = ProbePointerBits(&sure);
      printf("PTR BITS=%d, SURE=%d\n", bits, sure);
   }

/*
 * Here for future, presently unsupported
 */
   if (flags & Pncache)
      printf("NCACHES=0\n");
   if (flags & PCacheSize)
      printf("%d Cache size (kb) = 0\n", CacheLevel);
   exit(0);
}
