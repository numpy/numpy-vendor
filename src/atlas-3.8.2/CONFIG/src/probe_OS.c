#include "atlconf.h"

enum OSTYPE ProbeOS(int verb, char *targ)
{
   int ierr;
   char ln[1024], ln2[1024];
   enum OSTYPE i, OS;
   char *unam;

   if (verb) printf("Probing to make operating system determination:\n");
   unam = FindUname(targ);

   sprintf(ln2, "%s -s", unam);
   ierr = CmndOneLine(targ, ln2, ln);
   if (ierr == 0)
   {
      if(strstr(ln, "Linux")) OS = OSLinux;
      else if(strstr(ln, "FreeBSD")) OS = OSFreeBSD;
      else if (strstr(ln, "Darwin")) OS = OSOSX;
      else if(strstr(ln, "SunOS"))
      {
         sprintf(ln2, "%s -r", unam);
         CmndOneLine(targ, ln2, ln);
         if (ln[0] == '4') OS = OSSunOS4;
         else OS = OSSunOS;
      }
      else if(strstr(ln, "OSF1")) OS = OSOSF1;
      else if(strstr(ln, "IRIX")) OS = OSIRIX;
      else if(strstr(ln, "AIX")) OS = OSAIX;
      else if(strstr(ln, "WIN"))
      {
         if (strstr(ln, "95") || strstr(ln, "98") || strstr(ln, "_ME"))
            OS = OSWin9x;
/*
 *       Need to confirm what is returned under cygwin for XP, etc.
 */
         else if (strstr(ln, "NT")) OS = OSWinNT;  /* check this */
         else ierr = 1;
      }
      else if (strstr(ln, "HP-UX")) OS = OSHPUX;
      else ierr = 1;
   }
   if (ierr)
   {
      printf("\n\nUnable to determine OS, quitting\n\n");
      exit(-1);
   }
   if (verb)
      printf("Operating system configured as %s\n\n", osnam[OS]);

   return(OS);
}

void PrintUsage(char *name, int i)
{
   fprintf(stderr, "Error around argument %d\n", i);
   fprintf(stderr, "USAGE: %s -v <verbose #> -T <targ machine>\n", name);
   exit(-1);
}

int GetFlags(int nargs, char **args, char *targ)
{
   int verb=1, i;
   *targ = '\0';
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i);
      switch(args[i][1])
      {
      case 'v':
         if (++i >= nargs)
            PrintUsage(args[0], i);
         verb = atoi(args[i]);
         break;
      case 'T':                    /* target machine for spawn */
         if (++i >= nargs)
            PrintUsage(args[0], i);
         strcpy(targ, args[i]);
         break;
      default:
         PrintUsage(args[0], i);
      }
   }
   return(verb);
}

main(int nargs, char **args)
{
   int verb;
   char targ[1024];
   enum OSTYPE OS;
   verb = GetFlags(nargs, args, targ);
   OS = ProbeOS(verb, *targ == '\0' ? NULL : targ);
   printf("OS=%d\n", OS);
   exit(OS == OSOther);
}
