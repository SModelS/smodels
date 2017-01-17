.. index:: Installation and Deployment

Installation and Deployment
===========================

Standard Installation
---------------------

SModelS is a Python library that requires Python version 2.6 or later (but not
version 3) with Python's setuptools and gfortran installed. Internally, SModelS uses the
following tools:


 * `Pythia 6.4.27 <http://arxiv.org/abs/hep-ph/0603175>`_
 * `NLL-fast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_ 1.2 (7 TeV), 2.1 (8 TeV), and 3.1 (13 TeV)

These tools are built into SModelS and require gfortran, but they need not be installed separately.
In addition to setuptools, SModelS depends on the following *external* Python
libraries [*]_:

 * unum
 * numpy
 * argparse
 * docutils>=0.3
 * scipy>=0.9.0
 * pyslha>=3.1.0

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


Installation on Ubuntu 16.04
----------------------------

Installation on Ubuntu machines should be straightforward with superuser privileges
(if you do not have superuser privileges see instructions below):

 * sudo apt install gfortran python-setuptools python-scipy python-numpy python-docutils python-argparse
 * python setup.py install

Note that the last command can be run as superuser, or with the "--user" flag.


Installation on SL6
-------------------

Installation on an SL6 (Scientific Linux 6 or Scientific Linux CERN 6) machine
is tricky, because SModelS requires a more recent version of *scipy* than is provided by SL6.
We succeeded to install SModelS on SL6 by doing:

 * yum install gcc-c++ libstdc++-devel libevent-devel python-devel lapack lapack-devel blas blas-devel libgfortran python-distutils-extra

followed by:

 * pip install nose unum argparse numpy pyslha scipy

Note, that these steps can safely be done within a Python ``virtualenv``.
Pip can also be called with the "--user" flag.


Installation on SL5 and similar distributions
---------------------------------------------

In some distributions like SL5, the Python default version may be smaller than
2.6.  In these cases, ``virtualenv`` has to be set up for a Python version >=         2.6.  E.g. for Python 2.6, do ``virtualenv --python=python2.6 <envname>``,            and modify by hand the first line in the executable from ``#!/usr/bin/env python``
to ``#!/usr/bin/env python2.6``.
Then perform the steps listed under ``Installation on SL6``.



Installation on other platforms or without superuser privileges using Anaconda
------------------------------------------------------------------------------

Another easy and platform independent way of installing SModelS
without superuser priviledges is via Anaconda (https://www.continuum.io/downloads).
Anaconda provides a local installation of pip as well as several additional python packages.
Here we assume a version of gfortran is already installed in your system.

 * download and install Anaconda for Python 2.7 (https://www.continuum.io/downloads)
 * make sure Anaconda's bin and lib folders are added to your system and python paths ::

    PATH="<anaconda-folder>/bin:$PATH"
    PYTHONPATH=$PYTHONPATH:"<anaconda-folder>/lib/python2.7/site-packages"
 
and then install SModelS as a user: ::
 
 python setup.py install --user

In order to make sure all libraries have been correctly installed, you can run ::
   
 python smodels/tools/toolBox.py


Adding results to the database
------------------------------

.. _addingFastlim:

Adding FastLim data
^^^^^^^^^^^^^^^^^^^

The official SModelS database can be augmented with data from the
`fastlim <http://cern.ch/fastlim>`_ database.
A tarball with the *properly converted* fastlim-1.0 efficiency maps can be found in our
`download section <http://smodels.hephy.at/downloads/v1.1>`_.
The tarball then needs to be exploded in the top level directory of the database.

That is, the following steps need to be performed ::

 mv smodels-v1.1-fastlim-1.0.tgz <smodels-database folder>
 cd <smodels-database folder>
 tar -xzvf smodels-v1.1-fastlim-1.0.tgz
 rm smodels-v1.1-fastlim-1.0.tgz

Once the fastlim folders have been added to the database,
SModelS auto-detects fastlim results and issues an acknowledgement.
When using these results, please properly cite the fastlim paper; for
convenience, a bibtex file is provided in the smodels-fastlim tarball.


Adding one's own results
^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`Database of Experimental Results <databaseStruct>`  is
organized as files in an ordinary directory hierarchy. Therefore,
adding additional experimental results is a matter of copying and editing text
files.  
Once the new folders and files have been added following the
:ref:`database structure format <folderStruct>`, SModelS
automatically rebuilds the binary (Pickle) database file.
The added results will then be available for using with the
the SModelS tools.


.. [*] The :ref:`database browser <databaseBrowser>` interface provided by smodelsTools.py also
   requires IPython. However, all the other SModelS functionalities are independent of IPython.



 
