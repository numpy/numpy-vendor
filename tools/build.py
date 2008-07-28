#! /usr/bin/python

from subprocess import Popen
import shutil
import os.path
from os.path import join as pjoin, dirname, abspath, basename
import glob
import tarfile
from ConfigParser import ConfigParser
import sys

if os.name == 'nt' and not sys.platform == 'cygwin':
    raise ValueError("You should use cygwin python on windows for now !")

# Configuration (this should be put in a config file at some point)
# All path are relative to top directory (vendor)
LAPACK_SRC = pjoin('src', 'lapack-lite-3.1.1')
LAPACK_LIB = abspath(pjoin(LAPACK_SRC, 'lapack_MINGW32.a'))
BLAS_LIB = abspath(pjoin(LAPACK_SRC, 'blas_MINGW32.a'))
LAPACK_TARBALL = 'lapack-3.1.1.tbz2'
ATLAS_SRC = pjoin('src', 'atlas-3.8.2')
ATLAS_BUILDDIR = pjoin(ATLAS_SRC, "MyObj")
# Use INT_ETIME for lapack if building with gfortran
#TIMER = 'INT_ETIME'

def build_atlas_tarball():
    print "====== Building ATLAS tarbal ======"
    files = glob.glob(ATLAS_BUILDDIR + '/lib/*.a')
    fid = tarfile.open(ATLAS_TARBALL, 'w:bz2')
    try:
        for f in files:
            af = pjoin(ARCH.lower(), 'mingw32', basename(f))
            fid.add(f)
    finally:
        fid.close()

def build_lapack_tarball():
    print "====== Building BLAS/LAPACK tarbal ======"
    fid = tarfile.open(LAPACK_TARBALL, 'w:bz2')
    try:
        af = pjoin(ARCH.lower(), 'mingw32', 'libblas.a')
        fid.add(BLAS_LIB, af)
        af = pjoin(ARCH.lower(), 'mingw32', 'liblapack.a')
        fid.add(LAPACK_LIB, af)
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
    p = Popen(['../configure',  '-C', 'if', FC, "-b", str(ATLAS_PW), "-m", str(ATLAS_MHZ), '--with-netlib-lapack=%s' % LAPACK_LIB], cwd = ATLAS_BUILDDIR)
    os.waitpid(p.pid, 0)

def build_lapack():
    print "====== Build LAPACK ======"
    p = Popen(["make", "lapacklib", "FORTRAN=%s" % FC, "LOADER=%s" % FC, "OPTS=%s" % LAPACK_F77_FLAGS], cwd = LAPACK_SRC)
    os.waitpid(p.pid, 0)

def build_blas():
    print "====== Build BLAS ======"
    p = Popen(["make", "blaslib", "FORTRAN=%s" % FC, "LOADER=%s" % FC, "OPTS=%s" % LAPACK_F77_FLAGS], cwd = LAPACK_SRC)
    os.waitpid(p.pid, 0)

def clean_lapack():
    print "====== Clean LAPACK ======"
    p = Popen(['make', 'cleanlib'],  cwd = LAPACK_SRC)
    os.waitpid(p.pid, 0)

def clean_blas():
    print "====== Clean BLAS ======"
    p = Popen(['make', 'clean'],  cwd = LAPACK_SRC)
    os.waitpid(p.pid, 0)

def clean_atlas():
    print "====== Clean ATLAS ======"
    if os.path.exists(ATLAS_BUILDDIR):
        shutil.rmtree(ATLAS_BUILDDIR)

def clean():
    clean_atlas()
    clean_lapack()

TARGETS = {'atlas' : [configure_atlas, build_atlas, build_atlas_tarball],
        'lapack' : [build_lapack],
        'blas-lapack' : [build_blas, build_lapack, build_lapack_tarball]}

class Config(object):
    def __init__(self):
        self.arch = 'NOSSE'
        self.cpuclass = 'i386'
        self.freq = 0
        self.pw = 32
        self.targets = ['blas', 'lapack']
        self.lapack_flags = ""
        self.f77 = 'g77'

    def __repr__(self):
        r = ["Cpu Configurations: "]
        r += ['\tArch: %s' % self.arch]
        r += ['\tCpu Class: %s' % self.cpuclass]
        r += ['\tFreq: %d MHz' % self.freq]
        r += ['\tPointer width: %d bits' % self.pw]
        r += ["Targets to build: %s" % str(self.targets)]
        return "\n".join(r)

def read_config(file):
    cfgp = ConfigParser()
    f = cfgp.read(file)
    if len(f) < 1:
        raise IOError("file %s not found" % file)

    cfg = Config()
    if cfgp.has_section('CPU'):
        if cfgp.has_option('CPU', 'ARCH'):
            cfg.arch = cfgp.get('CPU', 'ARCH')
        if cfgp.has_option('CPU', 'CLASS'):
            cfg.cpuclass = cfgp.get('CPU', 'CLASS')
        if cfgp.has_option('CPU', 'MHZ'):
            cfg.freq = cfgp.getint('CPU', 'MHZ')
    if cfgp.has_section('BUILD_OPTIONS'):
        if cfgp.has_option('BUILD_OPTIONS', 'TARGETS'):
            cfg.targets = cfgp.get('BUILD_OPTIONS', 'TARGETS').split(',')
        if cfgp.has_option('BUILD_OPTIONS', 'F77'):
            cfg.f77 = cfgp.get('BUILD_OPTIONS', 'F77')
        if cfgp.has_option('BUILD_OPTIONS', 'LAPACK_F77FLAGS'):
            cfg.lapack_flags = " ".join(cfgp.get('BUILD_OPTIONS', 'LAPACK_F77FLAGS').split(','))

    return cfg

if __name__ == '__main__':
    argc = len(sys.argv)
    if argc < 2:
        cfg = Config()
    else:
        cfg = read_config(sys.argv[1])

    ARCH = cfg.arch
    FC = cfg.f77
    LAPACK_F77_FLAGS = cfg.lapack_flags
    ATLAS_TARBALL = 'atlas-3.8.2-%s.tbz2' % ARCH
    ATLAS_PW = cfg.pw
    ATLAS_MHZ = cfg.freq

    clean()
    for t in cfg.targets:
        target = TARGETS[t]
        for action in target:
            action()
