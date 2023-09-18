#!/bin/bash


install_dir=$PWD

cd $install_dir

ls | grep -b '\.sh$' | grep -v Makefile | xargs rm -r
