#! /usr/bin/env python

import os
from distutils.core import setup
from glob import glob
from distutils.extension import Extension

incdir_src = os.path.abspath("../../include")
incdir_build = os.path.abspath("../../include")
libdir = os.path.abspath("../../src/.libs")


## Configure the C++ extension and LHAPDF package
ext = Extension("lhapdf",
                ["lhapdf.cpp"],
                include_dirs = [incdir_src, incdir_build],
                extra_compile_args=["-I/home/supersymmetry/Documents/lpsc/test/smodels/smodels/lib/resummino/lhapdf/include"],
                library_dirs = [libdir],
                language = "C++",
                libraries = ["stdc++", "LHAPDF"])
setup(name = "LHAPDF",
      version = "6.3.0",
      ext_modules = [ext])
