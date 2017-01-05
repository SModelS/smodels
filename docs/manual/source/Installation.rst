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
If the Python libraries are installed in a system folder (as is the default behavior),
it will be necessary to run the install command with superuser privilege.
Alternatively, one can run setup.py with the argument "--user" ::

  python setup.py install --user

In case the compilation of SModelS fails, it is advised to try to compile
the tools manually, by issuing "make" in the *lib/* directory.
In case the installation of the external libraries fails, you can also try to install
them manually, then rerun setup.py.
For Ubuntu, SL6 machines and other platforms, a recipe is given below.

There is also a diagnostic tool available: ::

   python smodels/tools/toolBox.py

should list and check all internal tools (Pythia and NLL-fast) and external
(numpy, scipy, unum, ... ) dependencies.

In case everything fails, please contact smodels-users@lists.oeaw.ac.at



Installation on Ubuntu 16.04
----------------------------

 * sudo apt install python-scipy python-numpy python-docutils python-argparse
 * pip install unum pyslha
 * python setup.py install

Note that the last two commands you either run as superuser, or with the "--user" flag.


Installation on SL6
-------------------

Installation on an SL6 (Scientific Linux 6 or Scientific Linux CERN 6) machine
is tricky, because SModelS requires a more recent version of *scipy* than is provided by SL6.
We succeeded to install SModelS on SL6 by doing:

 * yum install gcc-c++ libstdc++-devel libevent-devel python-devel lapack lapack-devel blas blas-devel libgfortran python-distutils-extra

followed by:

 * pip install nose unum argparse numpy pyslha scipy

Note, that these steps can safely be done within a Python ``virtualenv``.
Pip can also be called with the "--user" option.


Installation on SL5 and similar distributions
---------------------------------------------

In some distributions like SL5, the Python default version may be smaller than
2.6.  In these cases, ``virtualenv`` has to be set up for a Python version >=         2.6.  E.g. for Python 2.6, do ``virtualenv --python=python2.6 <envname>``,            and modify by hand the first line in the executable from ``#!/usr/bin/env python``
to ``#!/usr/bin/env python2.6``.
Then perform the steps listed under ``Installation on SL6``.



Installation on other platforms or without superuser privileges using Anaconda
------------------------------------------------------------------------------

Another easy and platform independent way of installing SModelS
without superuser priviledges is via `Anaconda <https://www.continuum.io/downloads>`_ .
Anaconda provides a local installation of pip as well as several additional python packages.

First download and install Anaconda for Python 2.7. Make sure Anaconda's bin and lib folders
are added to your system and python paths ::

    PATH="<anaconda-folder>/bin:$PATH"
    PYTHONPATH=$PYTHONPATH:"<anaconda-folder>/lib/python2.7/site-packages"

Now you should be able to run the SModelS installation (in the smodels top folder) as a user ::

   python setup.py install --user

In order to make sure all libraries have been correctly installed, you can run ::

   python smodels/tools/toolBox.py


.. _addingFastlim:

Adding fastlim data
-------------------

The official SModelS database can be augmented with data from the
`fastlim <http://cern.ch/fastlim>`_ database.
A tarball with the *properly converted* fastlim efficiency maps can be found in our
`download section <http://smodels.hephy.at/downloads/v1.1>`_.
The tarball then needs to be exploded in the top level directory of the database.

That is, the following steps need to be performed:

 * cd smodels-database
 * wget http://smodels.hephy.at/downloads/v1.1.0/smodels-fastlim-v1.1.0.tgz
 * tar -xzvf smodels-fastlim-v1.1.0.tgz
 * rm smodels-fastlim-v1.1.0.tgz

SModelS auto-detects fastlim results and issues an acknowledgement.

Please make sure, that when using their efficiency maps, fastlim gets proper
acknowledgement, see the bibtex file in the smodels-fastlim tarball.

Adding one's own results
------------------------

As will be explained in :doc:`Database of Experimental Results <DatabaseStructure>`, the database of
experimental results is organized as files in an ordinary directory hierarchy.
Therefore, adding additional experimental results is a matter of copying and
editing text files.  The next time the
:ref:`Database <Database>` class is instantiated, the binary (Pickle) database
file is updated automatically.
