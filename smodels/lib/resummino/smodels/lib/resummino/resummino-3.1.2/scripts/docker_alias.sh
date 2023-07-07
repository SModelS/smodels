#!/usr/bin/env bash
VER="3.1.2-latest"
if ! command -v lhapdf-config &> /dev/null
then
    echo "<lhapdf-config> could not be found. Using limited preinstalled PDFs."
# Preinstalled pdfs are listed here https://github.com/APN-Pucky/Dockerfiles/tree/master/debian/lhapdf
# Create a pull request/issue for different pdfs
    alias resummino='docker run -i --rm -u `id -u $USER`:`id -g` -v $PWD:$PWD -w $PWD apnpucky/resummino:'$VER' resummino'
else
    alias resummino='docker run -i --rm -u `id -u $USER`:`id -g` -v $PWD:$PWD -v `lhapdf-config --datadir`:/usr/share/LHAPDF -w $PWD apnpucky/resummino:'$VER' resummino'
fi


