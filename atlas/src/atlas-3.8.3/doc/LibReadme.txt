***************************** FILE DESCRIPTIONS *******************************
The files in this archive are:

Make.inc      : The Make include file used to build these libs
SUMMARY.LOG   : The SUMMARY.LOG created by atlas_install.  Provides further
                information about the installation machine.
cblas.h       : The C header file for the C interface to the BLAS.
clapack.h     : The C header file for the C interface to LAPACK.
liblapack.a   : The LAPACK routines provided by ATLAS.  If you want a full
                lapack library, the .o in this lib can be archived into the
                f77 lapack lib without error.
libcblas.a    : The ANSI C interface to the BLAS.
libf77blas.a  : The Fortran77 interface to the BLAS.
libatlas.a    : The main ATLAS library, providing low-level routines for all
                interface libs.

Your archive may also contain additional libraries, if it has posix thread
support.  These optional libs are:

libptcblas.a    : The ANSI C interface to the threaded (SMP) BLAS.
libptf77blas.a  : The Fortran77 interface to the threaded (SMP) BLAS.

********************************* LINKING *************************************
When linking, remember that order is important.  So, if you want uniprocessor
libs, your link line would contain IN THIS ORDER:
   -LLIBDIR -llapack -lcblas -lf77blas -latlas
And if you want to utilize an SMP version, it would be:
   -LLIBDIR -llapack -lptcblas -lptf77blas -latlas

NOTE: On Apple's OS X, the above won't work, since system directories are
searched before the -L directories.  OS X includes ATLAS internally, so
the above link line will always get you OS X's libs instead of the ones
that you have built.  The easiest solution is to explicitly link to each
library using the full path, rather than using the -L/-l combo.

************************** GETTING A FULL LAPACK LIB **************************
ATLAS does not provide a full lapack library.  However, there is a simple way
to get ATLAS to provide its faster LAPACK routines to a full LAPACK library.
ATLAS's internal routines are distinct from LAPACK's, so it is safe to compile
ATLAS's LAPACK routines directly into a netlib-style LAPACK library.
If you install LAPACK before ATLAS, this can be automated as described in
ATLAS/INSTALL.txt.  The following describes adding netlib LAPACK to an
already built ATLAS library:

First, obtain the LAPACK src from netlib and build the LAPACK library as
normal.  Then, in the lib/ subdirectory of your ATLAS build directory
(where you should have a liblapack.a), issue the following commands:
  mkdir tmp
  cd tmp
  ar x ../liblapack.a
  cp <your LAPACK path & lib> ../liblapack.a
  ar d ../liblapack.a lsame.o xerbla.o csrot.o zdrot.o
  ar r ../liblapack.a *.o
  cd ..
  rm -rf tmp

Just linking in ATLAS's liblapack.a first will not get you the best LAPACK
performance, mainly because LAPACK's untuned ILAENV will be used instead
of ATLAS's slightly tuned one.
