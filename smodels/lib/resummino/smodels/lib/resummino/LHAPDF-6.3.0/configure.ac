## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([LHAPDF],[6.3.0],[lhapdf-support@cern.ch],[LHAPDF])
if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported: use a different --prefix argument])
fi
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE()
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])


## Library version flags (at least update the version comment with each new release version)
AC_DEFINE_UNQUOTED(LHAPDF_VERSION, "$PACKAGE_VERSION", "LHAPDF version string")
# TODO: improve this sed'ing so that it'll work with version digits > 9
#PACKAGE_VERSION_CODE=[`echo $PACKAGE_VERSION | sed -e 's/\./0/g' -e 's/\([0-9]\+\).*/\1/'`] # sed \+ doesn't work on OS X!
PACKAGE_VERSION_CODE=[`echo $PACKAGE_VERSION | sed -e 's/\./0/g' -e 's/\([0-9]*\).*/\1/'`]
PACKAGE_MAJOR_VERSION=[`echo $PACKAGE_VERSION | sed -e 's/^\(.\..\).*/\1/'`]
AC_DEFINE_UNQUOTED(LHAPDF_VERSION_CODE, $PACKAGE_VERSION_CODE, "LHAPDF version as an int")


## Mac OS X compatibility
AC_CHECK_TOOL(SWVERS, sw_vers)
if test x$SWVERS != x; then
  PROD_NAME=$($SWVERS -productName | cut -f 2 -d:)
fi
AM_CONDITIONAL(WITH_OSX, [test "$PROD_NAME" == "Mac OS X"])
if test "$PROD_NAME" == "Mac OS X"; then
  MACOSX_DEPLOYMENT_TARGET=$($SWVERS -productVersion | cut -f 1,2 -d.)
  AC_MSG_NOTICE([MACOSX_DEPLOYMENT_TARGET = $MACOSX_DEPLOYMENT_TARGET])
  AM_CXXFLAGS="$AM_CXXFLAGS -Dunix"
  AM_CXXFLAGS="$AM_CXXFLAGS -flat_namespace"
fi


## Checks and config headers
AC_LANG(C++)
AC_CONFIG_SRCDIR([src/GridPDF.cc])
AC_CONFIG_HEADERS([config/config.h include/LHAPDF/Version.h])


## Set default C++ optimisation level to -O3 and drop the -g debug flag
if test -z $CXXFLAGS; then
    CXXFLAGS='-O3'
fi

## MPI support wanted?
AC_ARG_ENABLE([mpi],
  AC_HELP_STRING([--enable-mpi], [enable MPI-safe file reading]),
  [enable_mpi=yes], [])

## Test for C++ compiler, possibly with MPI support
AX_PROG_CXX_MPI([test x$enable_mpi = xyes], [have_mpi=yes], [have_mpi=no])
if test x$enable_mpi = xyes; then
  if test x$have_mpi = xyes; then
    AC_DEFINE([HAVE_MPI], "1", [using MPI])
    using_mpi=yes
  else
    AC_MSG_FAILURE([MPI compiler not found.])
    using_mpi=no
  fi
fi
AM_CONDITIONAL(USING_MPI, [test x$using_mpi = xyes])

## Require C++11
AX_CXX_COMPILE_STDCXX([11], [noext], [mandatory])

## Substitute compiler flags for the lhapdf-config script
LHAPDF_CXX="$CXX"
AC_SUBST(LHAPDF_CXX)
LHAPDF_CXXFLAGS="$CXXFLAGS"
AC_SUBST(LHAPDF_CXXFLAGS)


## More standard build tool checks
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL


## Enable LHAGLUE compatibility fns for Fortran and the old C++ interface
AC_ARG_ENABLE([lhaglue],
  [AC_HELP_STRING(--disable-lhaglue, [build without LHAPDF5 compatibility routines])],
  [], [enable_lhaglue=yes])
AC_ARG_ENABLE([lhaglue-cxx],
  [AC_HELP_STRING(--disable-lhaglue-cxx, [build without LHAPDF5 C++ compatibility routines])],
  [], [enable_lhaglue_cxx=yes])
if test x$enable_lhaglue == xyes; then
  AC_MSG_NOTICE(Building LHAGLUE Fortran wrappers)
  if test x$enable_lhaglue_cxx == xyes; then
    AC_MSG_NOTICE(Building LHAGLUE C++ wrappers)
  else
    AC_MSG_NOTICE(Not building LHAGLUE C++ wrappers)
  fi
else
   AC_MSG_NOTICE(Not building LHAGLUE wrappers)
fi
AM_CONDITIONAL(ENABLE_LHAGLUE, [test x$enable_lhaglue == xyes])
if test x$enable_lhaglue == xyes; then flag=1; else flag=0; fi
AC_DEFINE_UNQUOTED(LHAPDF_PDFLIB, $flag, "Whether PDFLIB Fortran compatibility is available")
AM_CONDITIONAL(ENABLE_LHAGLUE_CXX, [test x$enable_lhaglue_cxx == xyes])
if test x$enable_lhaglue_cxx == xyes; then flag=1; else flag=0; fi
AC_DEFINE_UNQUOTED(LHAPDF_LHA5CXX, $flag, "Whether LHAPDF5 C++ compatibility is available")


## Find pdflatex if possible
AC_PATH_PROG(PDFLATEX, pdflatex)
AM_CONDITIONAL(WITH_PDFLATEX, [test x$PDFLATEX != x])


## Find Doxygen if possible
AC_PATH_PROG(DOXYGEN, doxygen)
AM_CONDITIONAL(WITH_DOXYGEN, [test x$DOXYGEN != x])


## Find Python and Cython if possible
AC_ARG_ENABLE([python],
  [AC_HELP_STRING(--disable-python, [disable the build of the Python interface (default=enabled)])],
  [], [enable_python=yes])
if test x$enable_python == xyes; then
  AC_PATH_PROG(PYTHON, python)
  if test x$PYTHON != x; then
    PYTHON_VERSION=`$PYTHON -c "import sys; print('.'.join(map(str, sys.version_info@<:@:2@:>@)))"`
    AC_SUBST(PYTHON_VERSION)
    PYTHON_PATH=`$PYTHON -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib(prefix='$prefix', plat_specific=True))"`
    AC_SUBST(PYTHON_PATH)
  fi
  AM_CHECK_CYTHON([0.19], [:], [:])
  if test x$CYTHON_FOUND = xyes; then
    AC_MSG_NOTICE([Cython >= 0.19 found: Python extension source can be rebuilt (for developers)])
  fi
fi
AM_CONDITIONAL(WITH_PYTHON, [test x$PYTHON != x])
AM_CONDITIONAL(WITH_CYTHON, [test x$CYTHON_FOUND = xyes])


dnl ## Build Octave interface if possible
dnl AC_ARG_ENABLE([octave],
dnl   [AC_HELP_STRING(--disable-octave, [don't try to make Octave interface])],
dnl   [], [enable_octave=yes])
dnl if test x$enable_octave == xyes; then
dnl    AC_PATH_PROG(OCTAVE, octave)
dnl    AC_PATH_PROG(MKOCTFILE, mkoctfile)
dnl fi
dnl AM_CONDITIONAL(WITH_OCTAVE, [test x$OCTAVE != x && test x$MKOCTFILE != x])


## Set base compiler flags
AC_CEDAR_CHECKCXXFLAG([-pedantic], [AM_CXXFLAGS="$AM_CXXFLAGS -pedantic "])
AC_CEDAR_CHECKCXXFLAG([-Wall], [AM_CXXFLAGS="$AM_CXXFLAGS -Wall "])
AC_CEDAR_CHECKCXXFLAG([-Wno-long-long], [AM_CXXFLAGS="$AM_CXXFLAGS -Wno-long-long "])
AC_CEDAR_CHECKCXXFLAG([-Qunused-arguments], [AM_CPPFLAGS="$AM_CPPFLAGS -Qunused-arguments "])


## Include $prefix in the compiler flags for the rest of the configure run
if test x$prefix != xNONE; then
  CPPFLAGS="$CPPFLAGS -I$prefix/include"
  LDFLAGS="$LDFLAGS -L$prefix/lib"
fi


## Export base compiler flags
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(AM_LDFLAGS)


## Output
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/yamlcpp/Makefile])
AC_CONFIG_FILES([include/Makefile include/LHAPDF/Makefile])
AC_CONFIG_FILES([bin/Makefile])
AC_CONFIG_FILES([bin/lhapdf-config])
AC_CONFIG_FILES([examples/Makefile])
AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([doc/Makefile doc/Doxyfile])
AC_CONFIG_FILES([lhapdf.pc])
AC_CONFIG_FILES([wrappers/Makefile])
AC_CONFIG_FILES([wrappers/python/Makefile wrappers/python/setup.py])
dnl AC_CONFIG_FILES([wrappers/octave/Makefile])
AC_OUTPUT

## Message about getting PDF sets
pkgdatadir=`eval "echo $datarootdir/$PACKAGE_TARNAME"`
echo
echo "****************************************************************"
echo "IMPORTANT INFORMATION ABOUT PDF SETS"
echo
echo "LHAPDF no longer bundles PDF set data in the package tarball."
echo "The sets are instead all stored online at"
echo "  http://lhapdfsets.web.cern.ch/lhapdfsets/current/"
echo "and you should install those that you wish to use into"
echo "  $pkgdatadir"
echo
echo "The downloadable PDF sets are packaged as tarballs, which"
echo "must be expanded to be used. The simplest way to do this is with"
echo "the 'lhapdf' script, e.g. to install the CT10nlo PDF set:"
echo "  lhapdf install CT10nlo"
echo "The same effect can be achieved manually with, e.g.:"
echo "  wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/CT10nlo.tar.gz -O- | tar xz -C $pkgdatadir"
echo
echo "****************************************************************"
echo
