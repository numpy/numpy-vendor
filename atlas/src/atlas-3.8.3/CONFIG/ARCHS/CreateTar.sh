#!/bin/sh
# Creates the appropriate .tgz from all subdirs
for arch in `ls -F | fgrep "/" | sed -e 's/\///'`
do
tar cvf $arch.tar $arch/*
gzip --best $arch.tar
mv $arch.tar.gz $arch.tgz
done
