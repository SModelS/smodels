#!/bin/sh

command -v gfortran && { exit 0; }
command -v ifort && { exit 0; }
command -v f77 && { exit 0; }
