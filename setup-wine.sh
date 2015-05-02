#! /bin/bash

set -ex

# This script sets up the Wine environment so that one can start compiling
# NumPy for Windows.
#
# Setup:
#
# apt-get install wine
# mkdir $HOME/repos
# cd $HOME/repos
# git clone https://github.com/certik/numpy-vendor
#
# The $tarballs directory should be equal to the numpy-vendor git repo.

tarballs=$HOME/repos/numpy-vendor
TMPDIR=/tmp/numpy-windows
rm -rf $TMPDIR
mkdir $TMPDIR

# The script below should just work without modifications:

rm -rf $HOME/.wine
# Use bogus DISPLAY, so that wine doesn't pop up any windows
export DISPLAY=:100

# The ALLUSERS=1 is essential --- the following line:
msiexec /i $tarballs/python-3.4.0.msi /qn ALLUSERS=1
msiexec /i $tarballs/python-3.3.5.msi /qn ALLUSERS=1
msiexec /i $tarballs/python-3.2.5.msi /qn ALLUSERS=1
msiexec /i $tarballs/python-2.7.6.msi /qn ALLUSERS=1
msiexec /i $tarballs/python-2.6.6.msi /qn ALLUSERS=1
# is then exactly equivalent to the following line:
#msiexec /i $tarballs/python-2.7.3.msi
# But doesn't require any user interaction. Without ALLUSERS=1, the installed
# Python will segfault on basic operations.

cd $HOME/.wine/drive_c/

# Install Atlas binaries:
mkdir -p local/lib/
cd local/lib
cp -r $tarballs/atlas/binaries/ atlas
ln -s atlas/binaries yop  # This is needed for older versions of NumPy
cd ../..

mkdir MinGW
cd MinGW
tar xzf $tarballs/binutils-2.20-1-mingw32-bin.tar.gz
tar xzf $tarballs/gcc-g77-3.4.5-20051220-1.tar.gz
tar xzf $tarballs/gcc-core-3.4.5-20051220-1.tar.gz
tar xzf $tarballs/mingw-runtime-3.10.tar.gz
tar xzf $tarballs/gcc-g++-3.4.5-20051220-1.tar.gz
tar xzf $tarballs/w32api-3.7.tar.gz

cp $tarballs/msvcr90.dll lib/
cp $tarballs/msvcr100.dll lib/
# This library is installed by Wine, but it doesn't work with our "objdump", so
# we remove it:
rm -f $HOME/.wine/drive_c/windows/winsxs/x86_microsoft.vc90.crt_1fc8b3b9a1e18e3b_9.0.30729.4148_none_deadbeef/msvcr90.dll

if [ ! -d "$HOME/.wine/drive_c/Program\ Files\ \(x86\)/" ]; then
    # On 32 bits, this directory is not created, but we expect it to be there:
    ln -s $HOME/.wine/drive_c/Program\ Files $HOME/.wine/drive_c/Program\ Files\ \(x86\)
fi

# Install NSIS
cd $HOME/.wine/drive_c/Program\ Files\ \(x86\)/
tar xzf $tarballs/NSIS.tar.gz


# Add MinGW and NSIS into the PATH:
cd $TMPDIR
cat > regtmp <<EOF
REGEDIT4

[HKEY_CURRENT_USER\Environment]
"PATH"="C:\\\\MinGW\\\\bin;C:\\\\Program Files (x86)\\\\NSIS"
EOF
wine regedit regtmp


# Install Nose:

tar xzf $tarballs/distribute-0.6.49.tar.gz
cd distribute-0.6.49
# The "clean" target here is essential, otherwise distribute fails to install
# properly in Python 3 (silently!) and then the Nose install fails at install
# time.
wine "C:\Python26\python" setup.py install clean
wine "C:\Python27\python" setup.py install clean
wine "C:\Python32\python" setup.py install clean
wine "C:\Python33\python" setup.py install clean
wine "C:\Python34\python" setup.py install clean
cd ..

tar xzf $tarballs/nose-1.1.2.tar.gz
cd nose-1.1.2
# Apply a patch to make Nose work with Python 3.2
# https://github.com/nose-devs/nose/issues/539
cat > nose.patch <<EOF
diff --git a/setup.py b/setup.py
index d8615c3..2b423f0 100644
--- a/setup.py
+++ b/setup.py
@@ -43,7 +43,7 @@ try:
     # This is required by multiprocess plugin; on Windows, if
     # the launch script is not import-safe, spawned processes
     # will re-run it, resulting in an infinite loop.
-    if sys.platform == 'win32':
+    if sys.platform == '_disabled_win32':
         import re
         from setuptools.command.easy_install import easy_install
EOF
patch -p1 < nose.patch


# Again, the "clean --all" target is essential (especially the "--all").
# Without it it (silently!) installs and then fails at import time.
wine "C:\Python26\python" setup.py install clean --all
wine "C:\Python27\python" setup.py install clean --all
wine "C:\Python32\python" setup.py install clean --all
wine "C:\Python33\python" setup.py install clean --all
wine "C:\Python34\python" setup.py install clean --all
cd ..

# Install Paver only in Python 2.7:

tar xzf $tarballs/Paver-1.2.2.tar.gz
cd Paver-1.2.2
wine "C:\Python27\python" setup.py install
cd ..

# Install SCons only in Python 2.7:

tar xzf $tarballs/scons-2.3.0.tar.gz
cd scons-2.3.0
wine "C:\Python27\python" setup.py install
cd ..
