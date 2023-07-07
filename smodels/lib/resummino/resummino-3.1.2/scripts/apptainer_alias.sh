#!/usr/bin/env bash
#VER="latest"
VER="3.1.2-latest"
PROG="apptainer"
if ! command -v apptainer &> /dev/null
then
  PROG="singularity"
  if ! command -v singularity &> /dev/null
  then  
    echo "Could not find singularity or apptainer."
    exit 1
  fi
fi
if ! command -v lhapdf-config &> /dev/null
then
    echo "<lhapdf-config> could not be found. Using limited preinstalled PDFs."
# Preinstalled pdfs are listed here https://github.com/APN-Pucky/Dockerfiles/tree/master/debian/lhapdf
# Create a pull request/issue for different pdfs
    alias resummino=$PROG' exec docker://apnpucky/resummino:'$VER' resummino'
else
    alias resummino=$PROG' exec -B `lhapdf-config --datadir`:/usr/share/LHAPDF:ro docker://apnpucky/resummino:'$VER' resummino'
fi


