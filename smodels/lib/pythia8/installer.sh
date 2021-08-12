#!/bin/sh

VER=`cat pythiaversion`
# https://pythia.org/download/pythia83/pythia8306.tgz
[ -z "$1" ] && { echo "pythia8 version is not given, use pythia$VER"; } || 
		{ VER="$1"; echo "set to version pythia$VER"; }

NCORES="4"
command -v nproc && NCORES=`nproc`

TARBALL="pythia$VER.tgz"
# URL=http://home.thep.lu.se/~torbjorn/pythia8/$TARBALL
URL=http://pythia.org/download/pythia83/$TARBALL
#TARBALL="pythia${VER}_fixed.tgz"

test -e .$TARBALL && { echo "[installer] tarball $TARBALL exists"; } || { echo "[installer] getting $TARBALL"; wget $URL 2>/dev/null || curl -O $URL; mv $TARBALL .$TARBALL; };

test -e pythia$VER && { echo "[installer] pythia$VER directory exists"; } || { echo "[installer] $TARBALL"; tar xvf .$TARBALL; };

test -e pythia$VER/lib/libpythia8.a && { echo "[installer] libpythia8.a exists"; } || { echo "[installer] building pythia8"; cd pythia$VER; ./configure; make -j $NCORES; };
