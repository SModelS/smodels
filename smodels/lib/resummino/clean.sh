#!/bin/bash


install_dir=$PWD
cd $install_dir

ls | grep -v '\.sh$' | grep -v Makefile | xargs rm -r
