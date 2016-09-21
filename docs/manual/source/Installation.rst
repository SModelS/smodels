.. index:: Installation and Deployment

Installation and Deployment
===========================

Standard Installation
---------------------

SModelS is a Python library that requires Python version 2.6 or later
(but not version 3).  Internally, SModelS uses the following tools:

 * `Pythia 6.4.27 <http://arxiv.org/abs/hep-ph/0603175>`_
 * `NLL-fast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_ 1.2 (7 TeV), 2.1 (8 TeV), and 3.1 (13 TeV)

These tools are built into SModelS, they need not be installed separately.
In addition, SModelS depends on the following *external* Python libraries:

 * unum
 * numpy
 * argparse
 * docutils>=0.3
 * scipy>=0.9.0
 * pyslha>=3.1.0

For installation, SModelS makes use of Python's *setuptools*.
Thus ::

  python setup.py install

should install the entire project, compile the internal Pythia and NLL-fast versions
using gfortran. It should also resolve the external dependencies, i.e. install
the Python libraries listed above using e.g. *pip*.
If the python libraries are installed in a system folder (as is the default behavior),
it will be necessary to run the install command with superuser privilege.
Alternatively, one can run setup.py with the argument "--user".

In case the compilation of SModelS fails, it is advised to try to compile
the tools manually, by issuing "make" in the *lib/* directory.
In case the installation of the external libraries fails, you can also try to install
them manually, then rerun setup.py.
For Ubuntu and SL6 machines, a recipe is given below.

There is also a diagnostic tool available: ::

   python smodels/tools/toolBox.py

should list and check all internal tools (Pythia and NLL-fast) and external
(numpy, scipy, unum, ... ) dependencies.

In case everything fails, please contact smodels-users@lists.oeaw.ac.at

Installation on Ubuntu 16.04
----------------------------
sudo apt install python-scipy python-numpy python-docutils python-argparse
pip install unum pyslha
python setup.py install
Note that the last two commands you either run as superuser, or with the
"--user" flag.

Installation on SL6
-------------------

Installation on an SL6 (Scientific Linux 6 or Scientific Linux CERN 6) machine
is tricky, because SModelS requires a more recent version of *scipy* than is provided by SL6.
We succeeded to install SModelS on SL6 by doing:

 * yum install gcc-c++ libstdc++-devel libevent-devel python-devel lapack \
               lapack-devel blas blas-devel libgfortran python-distutils-extra

followed by:

 * pip install nose unum argparse numpy pyslha scipy

Note, that these steps can safely be done within a Python ``virtualenv``.
Pip can also be called with the "--user" option.


Installation on SL5 and similar distributions
---------------------------------------------

In some distributions like SL5, the python default version may be smaller than
2.6.  In these cases, ``virtualenv`` has to be set up for a python version >=         2.6.  E.g. for python 2.6, do ``virtualenv --python=python2.6 <envname>``,            and modify by hand the first line in the executable from ``#!/usr/bin/env python``
to ``#!/usr/bin/env python2.6``.
Then perform the steps listed under ``Installation on SL6``.
