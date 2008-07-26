This directory contains the sources and tools necessary for building mandatory
dependencies of numpy/lapack, with a special focus on Windows.

The directory is organized as followed:
        * src: sources of the dependencies
        * tools: various scripts to build the binaries of the dependencies

In Particular, the tools directory contains:
        * a build.py script to build blas/lapack on Win32 from a configuration
        file.
        * basic.nsi: a NSIS script to build an installer for blas/lapack under
        different architectures.

The build.py tool
=================

dependencies:
-------------

cygwin with make, gcc, g77 and python

How to use:
-----------

python tools/build.py tools/sse3.cfg

Will build a blas/lapack using configuration in tools/sse3.cfg. In particular,
it contains informations on whether to build atlas or not, which fortran
compiler to use, etc...

The basic.nsi script
====================

You need NSIS to use it: http://nsis.sourceforge.net/Main_Page You also need a
small plugin cpu_caps, to tell NSIS wether the running CPU supports SSE2,
SSE3, etc...

It will look for binaries in the binaries directory: one directory per arch. A
subdir mingw32 is used for differentiatin with the not-yet supported Win64.

Note
----

Although the binaries themselves have to be built on windows (ATLAS is not
cross-compilable AFAIK), the installer itself can be built under linux. Debian
contains nsis ported on Linux:

makensis tools/basic.nsis
