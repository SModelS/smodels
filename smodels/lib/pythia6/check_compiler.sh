#!/bin/sh

command -v f77 && { exit 0; }
command -v gfortran && { exit 0; }
command -v ifort && { exit 0; }
