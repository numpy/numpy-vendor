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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )
#define Mmin(x, y) ( (x) > (y) ? (y) : (x) )
#define Mmax(x, y) ( (x) > (y) ? (x) : (y) )
int GetPower2(int n)
{
   int pwr2, i;

   if (n == 1) return(0);
   for (pwr2=0, i=1; i < n; i <<= 1, pwr2++);
   if (i != n) pwr2 = 0;
   return(pwr2);
}

#define ShiftThresh 2
char *GetDiv(int N, char *inc)
{
   static char ln[256];
   int pwr2 = GetPower2(N);
   if (N == 1) sprintf(ln, "%s", inc);
   else if (pwr2) sprintf(ln, "((%s) >> %d)", inc, pwr2);
   else sprintf(ln, "((%s) / %d)", inc, N);
   return(ln);
}

char *GetInc(int N, char *inc)
{
   static char ln0[256];
   char ln[256];
   char *p=ln;
   int i, n=N, iPLUS=0;

   if (n == 0)
   {
      ln[0] = '0';
      ln[1] = '\0';
   }
   while(n > 1)
   {
      for (i=0; n >= (1<<i); i++);
      if ( (1 << i) > n) i--;
      if (iPLUS++) *p++ = '+';
      sprintf(p, "((%s) << %d)", inc, i);
      p += strlen(p);
      n -= (1 << i);
   }
   if (n == 1)
   {
      if (iPLUS++) *p++ = '+';
      sprintf(p, "%s", inc);
   }
   if (iPLUS > ShiftThresh) sprintf(ln0, "(%d*(%s))", N, inc);
   else if (iPLUS) sprintf(ln0, "(%s)", ln);
   else sprintf(ln0, "%s", ln);
   return(ln0);
}

static void GetAij(int nY, int ngap, int gap, int gapmul, int ia,
                   int *icol, int *igap, int *gapi, int *imul)
{
   ia = ia % (nY * ngap * gap * gapmul);
   *imul = ia / (nY * ngap * gap);
   ia -= *imul * nY * ngap * gap;
   *gapi = ia / (nY * ngap);
   ia -= *gapi * nY * ngap;
   *igap = ia / nY;
   ia -= *igap * nY;
   *icol = ia;
}

static int GetAi(int ngap, int gap, int igap, int gapi, int imul)
{
   return(imul*gap*ngap+igap*gap+gapi);
}

static void mvTXbody(FILE *fpout, char *spc, int lat, int nY, int ngap, int gap,
                     int gapmul, int pfA, int pfX, int Yregs, int CLEANUP)
{
   int i, j, k, nops, xelts, op=0;
   int FETCHA=1, FETCHX=1;
   int ia=pfA, ix=pfX, iy=0, myregs;
   int icol, iacc, igap, imul, Acol, Aacc, agap, Amul, xcol, xacc, xgap, xmul;

   nops = gapmul * gap * ngap * nY;
   xelts = gapmul * gap * ngap;
   if (lat)
   {
      op = lat;
      ix += lat/nY;
      if (ix > xelts) ix = ix % xelts;
   }

   for (i=0; i < nops; i++, ia++, iy++)
   {
      if (ia == nops)
      {
         ia = 0;
         for (j=0; j < nY; j++)
         {
            if (xelts > 1) fprintf(fpout, "%s   pA%d += %d;\n", spc, j, xelts);
            else fprintf(fpout, "%s   pA%d++;\n", spc, j);
         }
         if (CLEANUP) FETCHA = 0;
      }
      if (ix == xelts)
      {
         ix = 0;
         if (xelts > 1) fprintf(fpout, "%s   X += %d;\n", spc, xelts);
         else fprintf(fpout, "%s   X++;\n", spc);
         if (CLEANUP) FETCHX = 0;
      }
      GetAij(nY, ngap, gap, gapmul, i, &icol, &iacc, &igap, &imul);
      GetAij(nY, ngap, gap, gapmul, ia, &Acol, &Aacc, &agap, &Amul);
      myregs = Yregs / nY;
      if (icol < Yregs % nY) myregs++;
      if (lat)  /* seperate multiply & add code */
      {
         fprintf(fpout, "%s   rY%d += rA%d;\n", spc,
                 icol+((iy/nY)%myregs)*nY, i %pfA);
         if (!CLEANUP || op < nops)
            fprintf(fpout, "%s   rA%d *= rX%d;\n", spc, op % pfA,
                    ((i+lat)/nY)%pfX);
      }
      else fprintf(fpout, "%s   rY%d += rA%d * rX%d;\n",
                   spc, icol+((iy/nY)%myregs)*nY, i%pfA, (i/nY)%pfX);
      op++;
      if (FETCHA)
      {
         k = GetAi(ngap, gap, Aacc, agap, Amul);
         if (k) fprintf(fpout, "%s   rA%d = pA%d[%d];\n", spc, i%pfA, Acol, k);
         else fprintf(fpout, "%s   rA%d = *pA%d;\n", spc, i%pfA, Acol);
      }
      if (op%nY==0 && FETCHX)
      {
         GetAij(1, ngap, gap, gapmul, ix, &xcol, &xacc, &xgap, &xmul);
         k = GetAi(ngap, gap, xacc, xgap, xmul);
         if (k) fprintf(fpout, "%s   rX%d = X[%d];\n",spc, ((i+lat)/nY)%pfX, k);
         else fprintf(fpout, "%s   rX%d = *X;\n", spc, ((i+lat)/nY)%pfX);
         ix++;
      }
   }
}

static void FetchY(FILE *fpout, char *spc, int pre, int beta, int nY, int Yregs,
                   char *breg)
{
   int i, j;

   if (beta != 0)
   {
      fprintf(fpout, "%s   rY0 = *Y;\n", spc);
      for (j=1; j < nY; j++) fprintf(fpout, "%s   rY%d = Y[%d];\n", spc, j, j);
   }
   if (beta != 0 && beta != 1)
   {
      fprintf(fpout, "%s   %s = beta;\n", spc, breg);
      for (j=0; j < nY; j++) fprintf(fpout, "%s   rY%d *= %s;\n", spc,j, breg);
   }
   if (Yregs > nY || beta == 0)
   {
      fprintf(fpout, "%s   ", spc);
      for (i=(beta!=0)*nY; i < Yregs; i++) fprintf(fpout, "rY%d = ", i);
      if (pre == 's' || pre == 'c') fprintf(fpout, "0.0f;\n");
      else fprintf(fpout, "0.0;\n");
   }
}

static void FetchAX(FILE *fpout, char *spc, int nY, int ngap, int gap,
                    int gapmul, int pfA, int pfX)
{
   int ix, ia, icol, iacc, igap, imul, k;

   for (ix=0; ix < pfX; ix++)
   {
      GetAij(1, ngap, gap, gapmul, ix, &icol, &iacc, &igap, &imul);
      k = GetAi(ngap, gap, iacc, igap, imul);
      if (k) fprintf(fpout, "%s   rX%d = X[%d];\n", spc, ix, k);
      else fprintf(fpout, "%s   rX%d = *X;\n", spc, ix);
   }
   for (ia = 0; ia < pfA; ia++)
   {
      GetAij(nY, ngap, gap, gapmul, ia, &icol, &iacc, &igap, &imul);
      k = GetAi(ngap, gap, iacc, igap, imul);
      if (k) fprintf(fpout, "%s   rA%d = pA%d[%d];\n", spc, ia, icol, k);
      else fprintf(fpout, "%s   rA%d = *pA%d;\n", spc, ia, icol);
   }
}

static void StartPipe(FILE *fpout, char *spc, int lat, int nY,
                      int ngap, int gap, int gapmul, int pfX)
{
   int i, k, icol, wgap, gapi, imul;

   for (i=0; i < lat; i++)
   {
      fprintf(fpout, "%s   rA%d *= rX%d;\n", spc, i, (i/nY)%pfX);
      if ((i+1) % nY == 0)
      {
         GetAij(1, ngap, gap, gapmul, pfX+i/nY, &icol, &wgap, &gapi, &imul);
         k = GetAi(ngap, gap, wgap, gapi, imul);
         if (k) fprintf(fpout, "%s   rX%d = X[%d];\n", spc, (i/nY)%pfX, k);
         else fprintf(fpout, "%s   rX%d = *X;\n", spc, (i/nY)%pfX);
      }
   }
}

static void CombY(FILE *fpout, char *spc, int nY, int Yregs)
/*
 * use binary tree for adding up multiple accumulators
 */
{
   int i, j, d, n=(Yregs+nY-1)/nY;
   for (d=1; d < n; d <<= 1)
      for (i=0; i < n; i += (d<<1))
         for (j=0; j < nY; j++)
            if (j+(i+d)*nY < Yregs)
               fprintf(fpout, "%s   rY%d += rY%d;\n", spc, j+i*nY, j+(i+d)*nY);
}

static void XCleanup(FILE *fpout, char *spc, char pre, int lat, int nY,
                     int ngap, int gap, int gapmul, int pfA, int pfX, int Yregs)
{
#if 1
   int mingap, j, ia, ix, iy, il, FirstTime=1;

   fprintf(fpout, "%s   if (X != stXN)\n%s   {\n", spc, spc);
   spc -= 3;
   if (pre == 's') mingap = 8;
   else mingap = 4;
   for (j=1; j < ngap*gap*gapmul; j <<= 1);
   if (j > ngap*gap*gapmul) /* not power of two */
   {
      if (GetPower2(gap))
      {
         if (!GetPower2(ngap)) ngap = 2;
         gapmul = j / (ngap*gap);
      }
      else
      {
         ngap = gapmul = 1;
         gap = j;
      }
   }
   if (j <= 2) FirstTime=0;
   for (; j; j >>= 1)
   {
      if (FirstTime)
         fprintf(fpout, "%s   if ( (ptrdiff_t)(stXN-X) < %d ) goto cu%d;\n",
                 spc, j, j>>1);
      else if (j == 1)
         fprintf(fpout, "%s   if ( (ptrdiff_t)(stXN-X) == 1 )\n", spc);
      else fprintf(fpout, "%s   if ( (ptrdiff_t)(stXN-X) >= %d )\n", spc, j);
      for (ia=pfA; (nY*ngap*gap*gapmul) % ia; ia--);
      for (ix=pfX; (ngap * gap * gapmul)% ix; ix--);
      if (lat)
      {
         for (il=Mmin(ia,lat); (nY*ngap*gap*gapmul) % il; il--);
         if (lat - il > 2 || il < 2 || il < pfA) il = 0;
      }
      else il = 0;
      iy = Mmin(Yregs, nY*ngap*gap*gapmul);
      fprintf(fpout,
"%s   { /* lat=%d, ngap=%d, gap=%d, gapmul=%d, pfA=%d, pfX=%d, yregs=%d */\n",
              spc, il, ngap, gap, gapmul, ia, ix, iy);
      spc -= 3;
      FetchAX(fpout, spc, nY, ngap, gap, gapmul, ia, ix);
      if (il) StartPipe(fpout, spc, il, nY, ngap, gap, gapmul, pfX);
      mvTXbody(fpout, spc, il, nY, ngap, gap, gapmul, ia, ix, iy, 1);
      spc += 3;
      fprintf(fpout, "%s   }\n", spc);
      if (FirstTime) fprintf(fpout, "cu%d:\n", j>>1);
      if (gapmul > 1) gapmul >>= 1;
      else if (ngap > 2) ngap >>= 1;
      else if (gap > mingap) gap >>= 1;
      else if (ngap > 1) ngap >>= 1;
      else if (gap > 1) gap >>= 1;
      FirstTime = 0;
   }
   spc += 3;
   fprintf(fpout, "%s   } /* done X cleanup */\n\n", spc);
#else
   int i;

   fprintf(fpout, "%s   if (X != stXN)\n%s   {\n", spc, spc);
   spc -= 3;
   if (nY == 1) i = 1;
   else for (i=pfA; nY % i; i--);
   FetchAX(fpout, spc, nY, 1, 1, 1, i, 1);

   fprintf(fpout, "%s   if (X != stXN_1)\n%s   {\n", spc, spc);
   spc -= 3;
   fprintf(fpout, "%s   do /* while (X != stXN_1) */\n%s   {\n", spc, spc);
   spc -= 3;
   mvTXbody(fpout, spc, lat, nY, 1, 1, 1, i, 1, nY, 0);
   spc += 3;
   fprintf(fpout, "%s   }\n%s   while(X != stXN_1);\n", spc, spc);
   spc += 3;
   fprintf(fpout, "%s   }\n", spc);

   mvTXbody(fpout, spc, lat, nY, 1, 1, 1, i, 1, nY, 1);
   spc += 3;
   fprintf(fpout, "%s   } /* finish cleanup */\n\n", spc);
#endif
}

static void Xloop(FILE *fpout, char *spc, char pre, int lat, int nY,
                  int ngap, int gap, int gapmul, int pfA, int pfX, int Yregs)
{
   int i, nu = ngap * gap * gapmul;

   fprintf(fpout, "%s   if (N >= %d)\n%s   {\n", spc, 2*nu, spc);
   spc -= 3;
   FetchAX(fpout, spc, nY, ngap, gap, gapmul, pfA, pfX);
   if (lat) StartPipe(fpout, spc, lat, nY, ngap, gap, gapmul, pfX);
   fprintf(fpout, "%s   do /* while (X != stX) */\n%s   {\n", spc, spc);
   spc -= 3;
   mvTXbody(fpout, spc, lat, nY, ngap, gap, gapmul, pfA, pfX, Yregs, 0);
   spc += 3;
   fprintf(fpout, "%s   }\n%s   while(X != stX);\n", spc, spc);
   mvTXbody(fpout, spc, lat, nY, ngap, gap, gapmul, pfA, pfX, Yregs, 1);
   spc += 3;
   fprintf(fpout, "%s   }\n\n", spc);

   XCleanup(fpout, spc, pre, lat, nY, ngap, gap, gapmul, pfA, pfX, Yregs);
   CombY(fpout, spc, nY, Yregs);
}

static void emit_mvT
   (FILE *fpout, char *rout, char pre, char *styp, char *typ, int lat,
    int beta, int nY, int ngap, int gap, int gapmul, int pfA, int pfX,
    int Yregs)
{
   char *spcs = "                                        ";
   char *spc = spcs+40;
   char ln[128];
   char *bnam[3] = {"b0", "b1", "bX"};
   int i, j, nu = ngap * gap * gapmul;
   if (beta != 0 && beta != 1) beta = 2;

   if (nY > 1) fprintf(fpout, "#include \"atlas_level1.h\"\n");
   fprintf(fpout, "#include <stddef.h>\n");
   if (rout == NULL)
   {
      rout = ln;
      sprintf(ln, "ATL_%cgemvT_a1_x1_%s_y1\n", pre, bnam[beta]);
   }
   fprintf(fpout, "void %s", rout);
   fprintf(fpout, "/*\n * lat=%d, beta=%s, nY=%d, ngap=%d, gap=%d, gapmul=%d, pfA=%d, pfX=%d\n */\n",
           lat, bnam[beta], nY, ngap, gap, gapmul, pfA, pfX);
   fprintf(fpout, "   (const int M, const int N, const %s alpha,\n", styp);
   fprintf(fpout,
           "   const %s *A, const int lda, const %s *X, const int incX,\n",
           typ, typ);
   fprintf(fpout, "   const %s beta, %s *Y, const int incY)\n", styp, typ);
   fprintf(fpout, "{\n");
   fprintf(fpout, "   register %s rX0", typ);
   for (i=1; i < pfX; i++) fprintf(fpout, ", rX%d", i);
   for (i=0; i < pfA; i++) fprintf(fpout, ", rA%d", i);
   for (i=0; i < Yregs; i++) fprintf(fpout, ", rY%d", i);
   fprintf(fpout, ";\n");
   fprintf(fpout, "   const %s *pA0=A", typ);
   for (i=1; i < nY; i++) fprintf(fpout, ", *pA%d=pA%d+lda", i, i-1);
   fprintf(fpout, ";\n");
   fprintf(fpout, "   const int n = %s;\n", GetInc(nu,GetDiv(nu, "N")));
   fprintf(fpout, "   const int m = %s;\n", GetInc(nY,GetDiv(nY, "M")));
   fprintf(fpout,
           "   const %s *stX = X + n - %d, *stXN = X + N, *stXN_1 = stXN-1;\n",
           typ, nu);
   fprintf(fpout, "   %s *stY = Y + m;\n", typ);
   fprintf(fpout, "%s   const int incA = %s - N;\n", spc, GetInc(nY, "lda"));
   fprintf(fpout, "\n");

   if (nY > 1)
   {
      fprintf(fpout, "%s   if (m)\n%s   {\n", spc, spc);
      spc -= 3;
   }
   fprintf(fpout, "%s   do /* while (Y != stY) */\n%s   {\n", spc, spc);
   spc -= 3;

   FetchY(fpout, spc, pre, beta, nY, Yregs, "rX0");

   Xloop(fpout, spc, pre, lat, nY, ngap, gap, gapmul, pfA, pfX, Yregs);

   fprintf(fpout, "%s   X -= N;\n", spc);
   for (j=0; j < nY; j++)
   {
      fprintf(fpout, "%s   pA%d += incA;\n", spc, j);
      fprintf(fpout, "%s   Y[%d] = rY%d;\n", spc, j, j);
   }
   if (nY > 1) fprintf(fpout, "%s   Y += %d;\n", spc, nY);
   spc += 3;
   if (nY == 1) fprintf(fpout, "%s   }\n%s   while(++Y != stY);\n", spc, spc);
   else fprintf(fpout, "%s   }\n%s   while(Y != stY);\n", spc, spc);
   if (nY > 1)
   {
      spc += 3;
      fprintf(fpout, "%s   }\n", spc);
      fprintf(fpout, "%s   if (m != M)\n%s   {\n", spc, spc);
      spc -= 3;
      if (beta != 0 && beta != 1) fprintf(fpout, "%s   rX0 = beta;\n", spc);
      if (nY > 2)
      {
         fprintf(fpout, "%s   stY += M-m;\n", spc);
         fprintf(fpout, "%s   do /* while (Y != stY) */\n%s   {\n", spc, spc);
         spc -= 3;
      }
      if (beta != 0) fprintf(fpout, "%s   rY0 = *Y;\n", spc);
      if (beta != 0 && beta != 1) fprintf(fpout, "%s   rY0 *= rX0;\n", spc);
      if (beta == 0) fprintf(fpout, "%s   rY0 = ", spc);
      else fprintf(fpout, "%s   rY0 += ", spc);
      fprintf(fpout, "ATL_%cdot(N, pA0, 1, X, 1);\n", pre);
      fprintf(fpout, "%s   *Y = rY0;\n", spc);
      if (nY > 2)
      {
         fprintf(fpout, "%s   pA0 += lda;\n", spc);
         spc += 3;
         fprintf(fpout, "%s   }\n%s   while (++Y != stY);\n", spc, spc);
      }
      spc += 3;
      fprintf(fpout, "%s   } /* end Y cleanup */;\n", spc);
   }
   fprintf(fpout, "}\n");
}

void PrintUsage(char *nam)
{
   fprintf(stderr, "USAGE: %s -l <latency> -Y <nY> -G <ngap> -g <gap> -M <gapmul> -A <# of regs for A> -X <# of regs for X> -y <# regs for Y> -f <file> -R <rout> -b <beta>\n", nam);
   exit(-1);
}

void GetFlags(int nargs, char **args, int *lat, int *beta, int *nY,
              int *ngap, int *gap, int *gapmul, int *pfA, int *pfX,
              int *Yregs, char *pre, char *styp, char *typ,
              FILE **fpout, char **rout)
{
   int i;
   *beta = -3;
   *lat = 0;
   *nY = *Yregs = *ngap = *gap = *gapmul = *pfA = *pfX = 1;
   *pre = 'd';
   *rout = NULL;
   *fpout = stdout;

   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'R':
         *rout = args[++i];
         break;
      case 'f':
         *fpout = fopen(args[++i], "w");
         assert(*fpout);
         break;
      case 'b':
         *beta = atoi(args[++i]);
         break;
      case 'l':
         *lat = atoi(args[++i]);
         break;
      case 'y':
         *nY = atoi(args[++i]);
         break;
      case 'G':
         *ngap = atoi(args[++i]);
         break;
      case 'g':
         *gap = atoi(args[++i]);
         break;
      case 'M':
         *gapmul = atoi(args[++i]);
         break;
      case 'A':
         *pfA = atoi(args[++i]);
         break;
      case 'X':
         *pfX = atoi(args[++i]);
         break;
      case 'Y':
         *Yregs = atoi(args[++i]);
         break;
      case 'p':
         i++;
         *pre = Mlowcase(args[i][0]);
         break;
      default:
         PrintUsage(args[0]);
      }
   }
   i = *nY * *ngap * *gap * *gapmul;
   assert(i % *pfA == 0);
   if (*lat)
   {
      assert(i % *lat == 0);
      assert(*pfX + *lat / *nY <= i);
      assert(*pfA > *lat);
   }
   assert((*ngap * *gap * *gapmul)% (*pfX)  == 0);
   assert(*Yregs >= *nY);

   switch (*pre)
   {
   case 'z':
      sprintf(styp, "double*");
      sprintf(typ, "double");
      break;
   case 'd':
      sprintf(styp, "double");
      sprintf(typ, "double");
      break;
   case 'c':
      sprintf(styp, "float*");
      sprintf(typ, "float");
      break;
   case 's':
      sprintf(styp, "float");
      sprintf(typ, "float");
      break;
   default:
      PrintUsage(args[0]);
   }
}

main(int nargs, char **args)
{
   char pre, styp[32], typ[32], *rout;
   int lat;     /* 0: combined muladd inst; X: seperate multiply & add, lat=X */
   int nY;      /* number of columns of A to operate on (# of dot products) */
   int Yregs;   /* number of registers to use for Y */
   int ngap;    /* number of accumulaters to use within each dot prod */
   int gap;     /* gap size in column */
   int gapmul;  /* number of unrollings of whole mess to do */
   int pfA;     /* number of registers to use in prefetching A */
   int pfX;     /* number of registers to use in prefetching X */
   int i, beta;
   FILE *fpout;

   GetFlags(nargs, args, &lat, &beta, &nY, &ngap, &gap, &gapmul, &pfA, &pfX,
            &Yregs, &pre, styp, typ, &fpout, &rout);

   if (beta == -3) /* generate all beta cases */
   {
      for (i=0; i < 3; i++)
      {
         if (i == 0 || i == 1) fprintf(fpout, "#ifdef BETA%d\n\n", i);
         else fprintf(fpout, "#ifdef BETAX\n\n");
         emit_mvT(fpout, rout, pre, styp, typ, lat, i, nY, ngap, gap, gapmul,
                  pfA, pfX, Yregs);
         fprintf(fpout, "\n\n#endif\n");
      }
   }
   else emit_mvT(fpout, rout, pre, styp, typ, lat, beta, nY, ngap, gap,
                 gapmul, pfA, pfX, Yregs);
   exit(0);
}
