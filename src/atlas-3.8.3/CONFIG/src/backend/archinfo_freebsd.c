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
      if (!CmndOneLine(NULL, "sysctl hw.model", res))
      {
         if (strstr(res, "PowerMac"))
         {
            if (strstr(res,"c1,2")||strstr(res,"c3,1")||strstr(res,"c3,2")||
                strstr(res,"c3,3")||strstr(res,"c3,4")||strstr(res,"c3,5")||
                strstr(res,"c3,6")||strstr(res,"c4,2")||strstr(res,"c4,5")||
                strstr(res,"c5,1"))
               mach = PPCG4;
            else if (strstr(res,"c11,2")|| strstr(res,"c12,1")||
                     strstr(res,"c7,2") || strstr(res,"c7,3") ||
                     strstr(res,"c8,1") || strstr(res,"c8,1") ||
                     strstr(res,"c8,2") || strstr(res,"c9,1"))
               mach = PPCG5;
         }
         else if (strstr(res, "PowerBook"))
         {
            if (strstr(res,"k3,2") || strstr(res,"k3,3") || strstr(res,"k3,4")||
                strstr(res,"k3,5") || strstr(res,"k5,1") || strstr(res,"k5,2")||
                strstr(res,"k5,3") || strstr(res,"k5,4") || strstr(res,"k5,5")||
                strstr(res,"k5,6") || strstr(res,"k5,7") || strstr(res,"k5,8")||
                strstr(res,"k5,9") || strstr(res,"k6,1") || strstr(res,"k6,2")||
                strstr(res,"k6,3") || strstr(res,"k6,4") || strstr(res,"k6,5")||
                strstr(res,"k6,7") || strstr(res,"k6,8"))
               mach = PPCG4;
         }
         else if (strstr(res, "RackMac"))
         {
            if (strstr(res, "c1,1") || strstr(res, "c1,2"))
               mach = PPCG4;
            else if (strstr(res, "c3,1"))
               mach = PPCG5;
         }
      }
      break;
   case AFSPARC: /* don't know */
      break;
   case AFALPHA:
      #if 0
      if (!CmndOneLine(NULL, "sysctl hw.model", res))
      {
         if (strstr(res, "433au")) mach = Dec21164;
         else if (strstr(res, "XP1000")) mach = Dec21264;
      }
      #endif
      break;
   case AFIA64: /* don't know */
      break;
   case AFX86:
      if (!CmndOneLine(NULL, "sysctl hw.model", res))
      {
         if (strstr(res, "Pentium Pro")) mach = IntPPRO;
         else if (strstr(res, "Pentium III")) mach = IntPIII;
         else if (strstr(res, "Pentium II ")) mach = IntPII;
         else if (strstr(res, "Athlon")) mach = AmdAthlon;
         else if (strstr(res, "AMD-K7")) mach = AmdAthlon;
         else if (strstr(res, "32 bit Hammer")) mach = AmdHammer;
         else if (strstr(res, "64 bit Hammer")) mach = AmdHammer;
         else if (strstr(res, "Pentium/P55C")) mach = IntP5MMX; /* sent by */
         else if (strstr(res, "Pentium")) mach=IntP5;       /* Nakata Maho */
         else if (strstr(res, "iMac4,1")) mach=IntCoreDuo;
      }
      break;
   default:;
   }
   return(mach);
}

int ProbeNCPU()
{
   int ncpu = 0;
   char *reslns, res[1024];

   if (!CmndOneLine(NULL, "sysctl hw.ncpu", res)) ncpu = GetLastInt(res);
   return(ncpu);
}

int ProbePointerBits(int *sure)
{
   int iret = 32;
   char *uname;
   char cmnd[1024], res[1024];

   *sure = 0;
/*
 * Note this is a weak probe, archinfo_x86 much better . . .
 */
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
   if (!CmndOneLine(NULL, "sysctl hw.cpufrequency", res) )
      mhz = GetFirstDouble(res) / 1000000;
   return(mhz);
}

int ProbeThrottle()
/*
 * RETURNS: 1 if cpu throttling is detected, 0 otherwise
 */
{
   int iret=0;
   int imax=0, imin=0, icur=0;
   char res[1024];

   if (!CmndOneLine(NULL, "sysctl hw.cpufrequency_max", res) )
      imax = GetFirstInt(res);
   if (!CmndOneLine(NULL, "sysctl hw.cpufrequency_min", res) )
      imin = GetFirstInt(res);
   if (imax)
   {
      if (imax != imin)
         iret = 1;
   }
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
