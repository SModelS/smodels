.. index:: Installation and Deployment

.. |EM| replace:: :ref:`EM-type <EMtype>`
.. |UL| replace:: :ref:`UL-type <ULtype>`
.. |EMr| replace:: :ref:`EM-type result <EMtype>`
.. |ULr| replace:: :ref:`UL-type result <ULtype>`
.. |EMrs| replace:: :ref:`EM-type results <EMtype>`
.. |ULrs| replace:: :ref:`UL-type results <ULtype>`
.. |database| replace:: :ref:`database <Database>`

Installation and Deployment
===========================

Standard Installation
---------------------

SModelS is a Python library that requires Python version 2.6 or later, including version 3, which is the default. It depends on the following *external* Python libraries:

 * unum>=4.0.0
 * numpy>=1.13.0
 * argparse
 * requests>=2.0.0
 * docutils>=0.3
 * scipy>=1.0.0
 * pyslha>=3.1.0

In addition, the :ref:`cross section computer <xsecCalc>` provided by :ref:`smodelsTools.py <smodelsTools>`
requires:

 * `Pythia 8.2 <https://arxiv.org/abs/1410.3012>`_ (requires a C++ compiler) or `Pythia 6.4.27 <http://arxiv.org/abs/hep-ph/0603175>`_ (requires gfortran)
 * `NLL-fast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_ 1.2 (7 TeV), 2.1 (8 TeV), and 3.1 (13 TeV) (requires a fortran compiler)

These tools need not be installed separately, as the SModelS build system takes care of that. The current default is that both Pythia6 and Pythia8 are installed together with NLLfast. However, the user can easily adapt the Makefile in the lib/ directory to fit his or her needs.
Finally, the :ref:`database browser <databaseBrowser>` provided by :ref:`smodelsTools.py <smodelsTools>`
requires `IPython <https://ipython.org/>`_, while the :ref:`interactive plotter <interactivePlots>` requires `plotly <https://plot.ly/python/>`_ and `pandas <https://pandas.pydata.org/>`_. 


Installation Methods
^^^^^^^^^^^^^^^^^^^^

.. _phenoInstallation:

1. The first installation method installs SModelS in the source directory, which can then be used for :ref:`running SModelS <runningSModelS>`.
   After downloading the source from the `SModelS releases page <https://github.com/SModelS/smodels/releases>`_
   and extracting it, run::

     make smodels

   in the top-level directory. The installation will remove redundant folders, install the required 
   dependencies (using pip install) and compile Pythia and NLL-fast. If the cross section computer is not needed, one can replace *smodels* with *smodels_noexternaltools* in the above command.
   In case the Python libraries can not be successfully
   installed, the user can install them separately using his/her preferred method. Pythia and NLL-fast can also be compiled separately
   running **make externaltools**.

2. If Python's *setuptools* is installed in your machine, SModelS and its dependencies
   can also be installed running::

     setup.py install

   within the main smodels directory. If the python libraries are installed in a system folder (as is the default behavior),
   it will be necessary to run the install command with superuser privilege.
   Alternatively, one can run setup.py with the "--user" flag: ::

     setup.py install --user

   If *setuptools* is not installed, you can try to install the external libraries
   manually and then rerun setup.py.
   For Ubuntu, SL6 machines and other platforms, a recipe is given below.


3. If *pip3* (or *pip*) is installed in your machine, installing SModelS should be as easy as: ::

     pip3 install smodels

   In this case, *gfortran* and *g++* need to be installed separately, if one
   wishes to compute cross sections with pythia6 and pythia8, respectively.
   Also, it might be necessary to perform::

     sudo smodelsTools.py fixpermissions

   in case of system-wide installs. User-specific installations on the other hand: ::

     pip3 install --user smodels

   will install smodels into the user's ~/.local directory. Depending on your
   platform, the environment variables $PATH, $PYTHONPATH, $LD_LIBRARY_PATH
   (or $DYLD_LIBRARY_PATH) might have to be set appropriately.


There is also a diagnostic tool available: ::

  smodelsTools.py toolbox

should list and check all internal tools (Pythia and NLL-fast) and external
(numpy, scipy, unum, ... ) dependencies.

In case everything fails, please contact smodels-users@lists.oeaw.ac.at


Installation on Ubuntu >= 16.04
-------------------------------

Installation on Ubuntu machines should be straightforward with superuser privileges
(if you do not have superuser privileges see instructions below):

 * sudo apt install gfortran python-setuptools python-scipy python-numpy python-docutils python-argparse
 * setup.py install

Note that the last command can be run as superuser, or with the "--user" flag.

Installation on SL7
-------------------

Installation on an SL7 or CentOS7 is straightforward:

 * yum install gcc-c++ scipy numpy

 * pip install unum pyslha argparse


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
2.6.  In these cases, ``virtualenv`` has to be set up for a Python version >=         2.6.  E.g. for Python 2.6, do ``virtualenv --python=python2.6 <envname>``,            and modify by hand the first line in the executable from ``#!/usr/bin/env python3``
to ``#!/usr/bin/env python2.6``.
Then perform the steps listed under ``Installation on SL6``.



Installation on other platforms or without superuser privileges using Anaconda
------------------------------------------------------------------------------

Another easy and platform independent way of installing SModelS
without superuser priviledges is via Anaconda (https://www.continuum.io/downloads).
Anaconda provides a local installation of pip as well as several additional python packages.
Here we assume a version of gfortran is already installed in your system.

 * download and install Anaconda for Python 3.6 (https://www.continuum.io/downloads)
 * make sure Anaconda's bin and lib folders are added to your system and Python paths ::

    PATH="<anaconda-folder>/bin:$PATH"
    PYTHONPATH=$PYTHONPATH:"<anaconda-folder>/lib/python3.6/site-packages"

and then install SModelS as a user: ::

 setup.py install --user

In order to make sure all libraries have been correctly installed, you can run: ::

    smodelsTools.py toolBox


Installation of the C++ interface
---------------------------------

SModelS v1.1.1 comes with a simple C++ interface, see the cpp directory.
Obviously, a C++ compiler is need, alongside with the python developers
(header) files (libpython-dev on ubuntu, python-devel on rpm-based distros).


Adding results to the database
------------------------------

The installation procedure explained above also installs SModelS'
:ref:`database of experimental results <databaseStruct>`
in the smodels-database subdirectory.
The complete list of analyses and results included in the database can be
consulted at `http://smodels.hephy.at/wiki/ListOfAnalyses <http://smodels.hephy.at/wiki/ListOfAnalyses>`_.
We note that all the results in the official database release have been
carefully validated  and the validation material can be
found at `http://smodels.hephy.at/wiki/Validation <http://smodels.hephy.at/wiki/Validation>`_.

The database can conveniently be updated independently from SModelS code
updates. It suffices to unpack any new database tarball and replace the database
directory. In the same fashion, one can easily add additional results as
explained below.


.. _addingFastlim:

Adding FastLim data
^^^^^^^^^^^^^^^^^^^

The official SModelS database can be augmented with data from the
`fastlim <http://cern.ch/fastlim>`_ database.
A tarball with the *properly converted* fastlim-1.0 efficiency maps can be found in our
download section at `http://smodels.hephy.at <http://smodels.hephy.at/>`_.
The tarball then needs to be exploded in the top level directory of the database.

That is, the following steps need to be performed ::

 mv smodels-v1.1-fastlim-1.0.tgz <smodels-database folder>
 cd <smodels-database folder>
 tar -xzvf smodels-v1.1-fastlim-1.0.tgz
 rm smodels-v1.1-fastlim-1.0.tgz

Efficiencies with a relative statistical uncertainty greater than 25%
we consider to be zero. Also, per default we discard zeroes-only results.
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

