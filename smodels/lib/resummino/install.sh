#!/bin/bash


install_dir=$PWD
LHAPDF_VERSION="6.5.4"
RESUMMINO_VERSION="3.1.2"


mkdir -p $install_dir

cd $install_dir


get_cpu_cores() {
  local num_cores=$(nproc)
  if [ "$num_cores" -ge 3 ]; then
    echo "$((num_cores / 2 - 1))"
  else
    echo "1"
  fi
}

num_cores_to_use=$(get_cpu_cores)



wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-$LHAPDF_VERSION.tar.gz -O "LHAPDF-$LHAPDF_VERSION.tar.gz
tar xf LHAPDF-6.3.0.tar.gz
cd LHAPDF-6.3.0
./configure --prefix=$install_dir/lhapdf --disable-python
make -j"$num_cores_to_use"
make install
cd ..
cd ..

wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/PDF4LHC21_40.tar.gz -O- | tar xz -C $install_dir/lhapdf/share/LHAPDF
cd $install_dir

wget https://resummino.hepforge.org/downloads/?f=resummino-$RESUMMINO_VERSION.zip -O resummino-$RESUMMINO_VERSION.zip


unzip resummino-$RESUMMINO_VERSION.zip
mkdir -p resummino_install

cd resummino-$RESUMMINO_VERSION
cmake . -DLHAPDF=$install_dir/lhapdf -DCMAKE_INSTALL_PREFIX=$install_dir/resummino_install
make -j"$num_cores_to_use"
make install
cd ..
