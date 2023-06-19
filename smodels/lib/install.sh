#!/bin/bash

smodels_dir=/home/supersymmetry/Documents/lpsc/smodels

install_dir=$PWD/lib

mkdir -p $install_dir

cd $install_dir

wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.3.0.tar.gz -O LHAPDF-6.3.0.tar.gz
tar xf LHAPDF-6.3.0.tar.gz
cd LHAPDF-6.3.0
./configure --prefix=$install_dir/lhapdf
make
make install
cd ..
cd ..


cd $install_dir

wget https://resummino.hepforge.org/downloads/?f=resummino-3.1.2.zip -O resummino-3.1.2.zip

sudo apt-get install libboost-dev libgsl-dev

unzip resummino-3.1.2.zip
mkdir -p resummino_install

cd resummino-3.1.2
cmake . -DLHAPDF=$install_dir/LHAPDF-6.3.0 -DCMAKE_INSTALL_PREFIX=$install_dir/resummino_install
make
make install
cd ..
