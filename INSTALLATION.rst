
                               ==============
                                SModelS v1.1
                               ==============


Installation
============

SModelS is a Python library that requires Python version 2.6 or later with the Python packages setuptools, unum, numpy, argparse, docutils (>=0.3), scipy (>=0.9.0),
and pyslha (>=3.1.).
if SModelS is to compute cross sections, it uses the following tools internally:

 * Pythia 8.226 (requires a C++ compiler)
 * Pythia 6.4.27 (obsolete)
 * NLL-fast 1.2 (7 TeV), 2.1 (8 TeV), and 3.1 (13 TeV) (requires a gfortran compiler)

These tools are built into SModelS and require gfortran and C++.
They need not be installed separately, as the SModelS build system takes care of that. As per default, NLL-fast and both Pythia6 and Pythia8 are installed. However, the user can easily adapt the Makefile, to fit his or her needs.

(The database browser interface provided by smodelsTools.py also requires IPython.
However, all the other SModelS functionalities are independent of IPython.)

If Python's *setuptools* is installed in your machine, SModelS and its dependencies
can be installed with::
  python setup.py install

If the python libraries are installed in a system folder (as is the default behavior),
it will be necessary to run the install command with superuser privilege.
Alternatively, one can run setup.py with the "--user" flag::
  python setup.py install --user

If *setuptools* is not installed, you can try to install the external libraries
manually and then rerun setup.py.
For Ubuntu, SL6 machines and other platforms, a recipe is given below.

There is also a diagnostic tool available: ::

   python smodels/tools/toolBox.py

should list and check all internal tools (Pythia and NLL-fast) and external
(numpy, scipy, unum, ... ) dependencies.

In case everything fails, please contact smodels-users@lists.oeaw.ac.at


Installation on Ubuntu >=16.04
------------------------------

Installation on Ubuntu machines should be straightforward with superuser privileges
(if you do not have superuser privileges see instructions below):

 * sudo apt install gfortran python-setuptools python-scipy python-numpy python-docutils python-argparse
 * python setup.py install
Note that the last command you either run as superuser, or with the
"--user" flag.

Installation on SL7
-------------------

Installation on an SL7 or CentOS7 is straightforward:

 * yum install gcc-c++ scipy numpy

 * pip install unum pyslha argparse

Installation on SL6
-------------------

Installation on an SL6 (Scientific Linux 6 or Scientific Linux CERN 6) machine
is tricky, because SModelS requires a more recent version of *scipy* (>=0.9.0)
than is provided by SL6 (0.7.2).
SModelS can be installed on SL6 by doing:

 * yum install gcc-c++ libstdc++-devel libevent-devel python-devel lapack \
               lapack-devel blas blas-devel libgfortran python-distutils-extra

followed by:

 * pip install nose unum argparse numpy pyslha scipy

Note, that these steps can safely be done within a Python ``virtualenv``.
Pip can also be called with the "--user" flag.


Installation on SL5 and similar distributions
---------------------------------------------

In some distributions like SL5, the python default version may be smaller than
2.6.  In these cases, ``virtualenv`` has to be set up for a python version >=
2.6.  E.g. for python 2.6, do ``virtualenv --python=python2.6 <envname>``, and
modify by hand the first line in the executable from ``#!/usr/bin/env python``
to ``#!/usr/bin/env python2.6``.
Then perform the steps listed under ``Installation on SL6``.


Installation on other platforms or without superuser privileges using Anaconda
----------------------------------------------------------------------------------

Another easy and platform independent way of installing SModelS
without superuser priviledges is via Anaconda (https://www.continuum.io/downloads).
Anaconda provides a local installation of pip as well as several additional python packages.
Here we assume a version of gfortran is already installed in your system.

 * download and install Anaconda for Python 2.7 (https://www.continuum.io/downloads)
 * make sure Anaconda's bin and lib folders are added to your system and python paths: ::

    PATH="<anaconda-folder>/bin:$PATH"
    PYTHONPATH=$PYTHONPATH:"<anaconda-folder>/lib/python2.7/site-packages"

and then install SModelS as a user:

 * python setup.py install --user

In order to make sure all libraries have been correctly installed, you can run

 * python smodels/tools/toolBox.py


Installation of C++ interface
-----------------------------

SModelS v1.1.1 comes with a simple C++ interface, see the cpp directory.
Obviously, a C++ compiler is needed, alongside with the python developers
(header) files (libpython-dev on ubuntu, python-devel on rpm-based distros).
