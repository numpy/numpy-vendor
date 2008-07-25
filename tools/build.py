#! /usr/bin/python

from subprocess import Popen
import shutil
import os.path
from os.path import join as pjoin
import glob
import tarfile

# Configuration (this should be put in a config file at some point)
LAPACK_SRC = pjoin('src', 'lapack-lite-3.1.1')
LAPACK_LIB = pjoin(LAPACK_SRC, 'lapack_MINGW32.a')
ATLAS_SRC = pjoin('src', 'atlas-3.8.2')
ATLAS_BUILDDIR = os.path.join(ATLAS_SRC, "MyObj")
ARCH = "SSE2-P4"
ATLAS_TARBALL = 'atlas-3.8.2-%s.tbz2' % ARCH
FC = 'g77'
# Use INT_ETIME for lapack if building with gfortran
#TIMER = 'INT_ETIME'
LAPACK_F77_FLAGS = ''

def build_atlas_tarball():
    print "====== Building ATLAS tarbal ======"
    files = glob.glob(ATLAS_BUILDDIR + '/lib/*.a')
    fid = tarfile.open(ATLAS_TARBALL, 'w:bz2')
    try:
        for f in files:
            fid.add(f)
    finally:
        fid.close()

def build_atlas():
    print "====== Building ATLAS ======"
    p = Popen(['make'], cwd = ATLAS_BUILDDIR)
    os.waitpid(p.pid, 0)

def configure_atlas():
    print "====== Configuring ATLAS ======"
    if os.path.exists(ATLAS_BUILDDIR):
        shutil.rmtree(ATLAS_BUILDDIR)
    os.makedirs(ATLAS_BUILDDIR)
    p = Popen(['../configure',  '--with-netlib-lapack=%s' % LAPACK_LIB, '-C', 'if', FC], cwd = ATLAS_BUILDDIR)
    os.waitpid(p.pid, 0)

def build_lapack():
    print "====== Build LAPACK ======"
    p = Popen(["make", "lapacklib", "FORTRAN=%s" % FC, "LOADER=%s" % FC, "OPTS=%s" % LAPACK_F77_FLAGS], cwd = LAPACK_SRC)
    os.waitpid(p.pid, 0)

def clean_lapack():
    print "====== Clean LAPACK ======"
    p = Popen(['make', 'cleanlib'],  cwd = LAPACK_SRC)
    os.waitpid(p.pid, 0)

def clean_atlas():
    print "====== Clean ATLAS ======"
    if os.path.exists(ATLAS_BUILDDIR):
        shutil.rmtree(ATLAS_BUILDDIR)

def clean():
    clean_atlas()
    clean_lapack()

if __name__ == '__main__':
    clean()
    build_lapack()
    configure_atlas()
    build_atlas()
    build_atlas_tarball()
