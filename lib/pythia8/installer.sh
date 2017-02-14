#!/bin/sh

VER='8223'
TARBALL="pythia$VER.tgz"
URL=http://home.thep.lu.se/~torbjorn/pythia8/$TARBALL

test -e $TARBALL && { echo "[installer] tarball $TARBALL exists"; } || { echo "[installer] getting $TARBALL"; wget $URL 2>/dev/null || curl -O $URL; };

test -e pythia$VER && { echo "[installer] pythia$VER directory exists"; } || { echo "[installer] $TARBALL"; tar xvf $TARBALL; };

test -e pythia$VER/lib/libpythia8.a && { echo "[installer] libpythia8.a exists"; } || { echo "[installer] building pythia8"; cd pythia$VER; ./configure; make; };
