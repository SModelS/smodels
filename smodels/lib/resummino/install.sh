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

cat <<EOF > boost_check.cpp
#include <iostream>
#ifdef BOOST_VERSION
#include <boost/version.hpp>
#endif

int main() {
    #ifdef BOOST_VERSION
    std::cout << "installed" << std::endl;
    #else
    std::cout << "not_installed" << std::endl;
    #endif
    return 0;
}
EOF

g++ boost_check.cpp -o boost_check

output=$(./boost_check)

if [ "$output" = "installed" ]; then
    echo "Boost is installed. Continuing script execution."
else
    echo "Boost is not installed. Stopping script execution."
    rm boost_check.cpp boost_check
    exit 1
fi

# Nettoyer les fichiers temporaires si Boost est installé
rm boost_check.cpp boost_check

if command -v gsl-config &>/dev/null; then
    echo "GSL is installed."
    GSL_VERSION=$(gsl-config --version)
    echo "GSL version: $GSL_VERSION"
else
    echo "GSL is not installed. Stopping script execution."
    exit 1
fi

num_cores_to_use=$(get_cpu_cores)

# Vérification de l'existence de LHAPDF
if [ ! -d "$install_dir/lhapdf" ]; then
    wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-$LHAPDF_VERSION.tar.gz -O "LHAPDF-$LHAPDF_VERSION.tar.gz"
    tar xf LHAPDF-$LHAPDF_VERSION.tar.gz
    cd LHAPDF-$LHAPDF_VERSION
    ./configure --prefix=$install_dir/lhapdf --disable-python
    make -j"$num_cores_to_use"
    make install
    cd ..
    cd ..
    wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/PDF4LHC21_40.tar.gz -O- | tar xz -C $install_dir/lhapdf/share/LHAPDF
    cd $install_dir
fi

# Vérification de l'existence de RESUMMINO
if [ ! -d "$install_dir/resummino_install" ]; then
    wget https://resummino.hepforge.org/downloads/?f=resummino-$RESUMMINO_VERSION.zip -O resummino-$RESUMMINO_VERSION.zip
    unzip resummino-$RESUMMINO_VERSION.zip
    mkdir -p resummino_install
    cd resummino-$RESUMMINO_VERSION
    cmake . -DLHAPDF=$install_dir/lhapdf -DCMAKE_INSTALL_PREFIX=$install_dir/resummino_install
    make -j"$num_cores_to_use"
    make install
    cd ..
fi

echo "LHAPDF_version = $LHAPDF_VERSION" > versions.txt
echo "resummino_version = $RESUMMINO_VERSION" >> versions.txt
