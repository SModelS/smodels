#!/bin/bash


install_dir=$PWD/smodels/lib/resummino

cd $install_dir
rm -r -v !("install.sh" | "clean.sh")
