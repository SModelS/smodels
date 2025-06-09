#!/bin/bash

install_dir=$PWD
LHAPDF_VERSION="6.5.5"
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
#include <boost/version.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

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

# Clean temporary files if Boost is installed
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

download_file() {
    url="$1"
    output="$2"

    if command -v wget > /dev/null; then
        echo "Using wget to download $url"
        wget "$url" -O "$output"
    elif command -v curl > /dev/null; then
        echo "Using curl to download $url"
        curl -L "$url" -o "$output"
    else
        echo "Error: Neither curl nor wget is available for downloading files."
        exit 1
    fi
}

download_and_install_lhapdf() {
    cd $install_dir
    download_file "$1" "LHAPDF-$LHAPDF_VERSION.tar.gz"
    tar xzf "LHAPDF-$LHAPDF_VERSION.tar.gz"
    if [ ! -d "LHAPDF-$LHAPDF_VERSION" ]; then
        echo "LHAPDF-$LHAPDF_VERSION directory missing, something went wrong. exiting."
        exit -1
    fi
    cd "LHAPDF-$LHAPDF_VERSION"
    ./configure --prefix=$install_dir/lhapdf --disable-python
    make -j"$num_cores_to_use"
    make install
    cd $install_dir
    download_file "http://lhapdfsets.web.cern.ch/lhapdfsets/current/PDF4LHC21_40.tar.gz" "PDF4LHC21_40.tar.gz"
    tar xz -C $install_dir/lhapdf/share/LHAPDF -f PDF4LHC21_40.tar.gz
    cd $install_dir
}

if [ ! -d "$install_dir/lhapdf" ]; then
    download_and_install_lhapdf "https://smodels.github.io/resummino/LHAPDF-$LHAPDF_VERSION.tar.gz"
#    download_file "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-$LHAPDF_VERSION.tar.gz" "LHAPDF-$LHAPDF_VERSION.tar.gz"
fi

if [ ! -d "$install_dir/lhapdf" ]; then
    echo "Failed to download from smodels.github.io, trying lhapdf.hepforge.org ..."
    download_and_install_lhapdf "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-$LHAPDF_VERSION.tar.gz"
    # we could also try to download from gitlab, but then we need to use the autotools, and we dont want that:
    # download_and_install_lhapdf "https://gitlab.com/hepcedar/lhapdf/-/archive/lhapdf-$LHAPDF_VERSION/lhapdf-lhapdf-$LHAPDF_VERSION.tar.gz"
fi

# Checking for the existence of RESUMMINO
download_and_install_resummino() {
    download_file "$1" "resummino-$RESUMMINO_VERSION.zip"
    if [ -f "resummino-$RESUMMINO_VERSION.zip" ]; then
        unzip "resummino-$RESUMMINO_VERSION.zip"
        mkdir -p resummino_install
        cd "resummino-$RESUMMINO_VERSION" || { echo "Error: Unable to change directory to resummino-$RESUMMINO_VERSION"; return 1; }
        cmake . -DLHAPDF=$install_dir/lhapdf -DCMAKE_INSTALL_PREFIX=$install_dir/resummino_install
        make -j"$num_cores_to_use"
        make install
        cd ..
		    cp $install_dir/resummino-$RESUMMINO_VERSION/bin/resummino $install_dir/resummino_install/bin/
        return 0
    else
        return 1
    fi
}

if [ ! -d "$install_dir/resummino_install" ]; then
    if ! download_and_install_resummino "https://smodels.github.io/resummino/resummino-$RESUMMINO_VERSION.zip"; then
        echo "Failed to download from smodels.github.io, trying hepforge.org..."
        download_and_install_resummino "https://resummino.hepforge.org/downloads/?f=resummino-$RESUMMINO_VERSION.zip"
    fi
fi

echo "# Note that this file gets overwritten when calling install.sh" > versions.txt 
echo "# so do not define the versions here!" >> versions.txt
echo "# Instead, look at install.sh." >> versions.txt
echo "LHAPDF_version = $LHAPDF_VERSION" >> versions.txt
echo "resummino_version = $RESUMMINO_VERSION" >> versions.txt
