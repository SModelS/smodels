# Copyright 2014 David R. Lamprea.
# Licensed under the European Union Public Licence (EUPL) 1.1.
#
# Exports the following variables:
# * LHAPDF_FOUND
# * LHAPDF_INCLUDE_DIRS - Set of paths to all required headers
# * LHAPDF_LIBRARIES - Set of all required libraries

find_library(LHAPDF_LIBRARIES LHAPDF HINTS ${LHAPDF}/lib64 ${LHAPDF}/lib $ENV{LHAPDF}/lib64 $ENV{LHAPDF}/lib ${LHAPDF_LIBDIR})

find_path(LHAPDF_INCLUDE_DIRS LHAPDF/LHAPDF.h HINTS ${LHAPDF}/include $ENV{LHAPDF}/include ${LHAPDF_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LHAPDF DEFAULT_MSG LHAPDF_LIBRARIES LHAPDF_INCLUDE_DIRS)
