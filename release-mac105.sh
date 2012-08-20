#! /bin/sh
# Script to build OS X installers on OS X
# Copy this script to a Mac 10.5 box and run it.

set -e

# Note that we build the corresponding set of OS X binaries to the python.org
# downloads, i.e. two versions for Python 2.7. The Intel 32/64-bit version is
# for OS X 10.6+, the other dmg installers are for 10.3+ and are built on 10.5

#---------------
# Build tarballs
#---------------
paver sdist


#--------------------
# Build documentation
#--------------------
# Check we're using the correct g++/c++ for the 32-bit 2.6 version we build for
# the docs and the 64-bit 2.7 dmg installer.
# We do this because for Python 2.6 we use a symlink on the PATH to select
# /usr/bin/g++-4.0, while for Python 2.7 we need the default 4.2 version.
export PATH=~/Code/tmp/gpp40temp/:$PATH
gpp="$(g++ --version | grep "4.0")"
if [ -z "$gpp" ]; then
    echo "Wrong g++ version, we need 4.0 to compile scipy with Python 2.6"
    exit 1
fi

# bootstrap needed to ensure we build the docs from the right scipy version
paver bootstrap
source bootstrap/bin/activate
python setup.py install

paver dmg -p 2.7   # 32/64-bit version

paver write_release_and_log
