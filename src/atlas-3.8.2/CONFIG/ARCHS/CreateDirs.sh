#!/bin/sh
# Creates dirs from .tgz
for arch in $*
do
gunzip -c $arch | tar xvf -
done
