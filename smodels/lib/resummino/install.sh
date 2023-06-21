#!/bin/bash


install_dir=$PWD/smodels/lib/resummino

mkdir -p $install_dir

cd $install_dir

wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.3.0.tar.gz -O LHAPDF-6.3.0.tar.gz
tar xf LHAPDF-6.3.0.tar.gz
cd LHAPDF-6.3.0
./configure --prefix=$install_dir/lhapdf --disable-python
make
make install
cd ..
cd ..

wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/PDF4LHC21_40.tar.gz -O- | tar xz -C $install_dir/lhapdf/share/LHAPDF
cd $install_dir

wget https://resummino.hepforge.org/downloads/?f=resummino-3.1.2.zip -O resummino-3.1.2.zip

sudo apt-get install libboost-dev libgsl-dev

unzip resummino-3.1.2.zip
mkdir -p resummino_install

cd resummino-3.1.2
cmake . -DLHAPDF=$install_dir/lhapdf -DCMAKE_INSTALL_PREFIX=$install_dir/resummino_install
make
make install
cd ..
