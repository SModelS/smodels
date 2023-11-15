#!/bin/bash


install_dir=$PWD
cd $install_dir

LINES=$(ls | grep -v '\.sh$' | grep -v Makefile | grep -v versions.txt )
# echo "lines x${LINES}x"
[ -z "$LINES" ] || {
    echo $LINES | xargs rm -r;
}
