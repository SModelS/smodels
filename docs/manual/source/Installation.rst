.. index:: Installation and Deployment

.. |EM| replace:: :ref:`EM-type <EMtype>`
.. |UL| replace:: :ref:`UL-type <ULtype>`
.. |EMr| replace:: :ref:`EM-type result <EMtype>`
.. |ULr| replace:: :ref:`UL-type result <ULtype>`
.. |EMrs| replace:: :ref:`EM-type results <EMtype>`
.. |ULrs| replace:: :ref:`UL-type results <ULtype>`
.. |database| replace:: :ref:`database <Database>`
.. |parameters| replace:: :ref:`parameters file <parameterFile>`

Installation and Deployment
===========================

Standard Installation
---------------------

SModelS is a Python library that requires Python version 3.6 or later. It depends on the following *external* Python libraries:

.. include:: dependencies.rst

For speed reasons, we moreover recommend pytorch>=1.8.0 as backend for pyhf. This is, however, optional: if pytorch is not available, SModelS will use the default backend. 

In addition, the :ref:`cross section computer <xsecCalc>` provided by :ref:`smodelsTools.py <smodelsTools>`
requires:

 * `Pythia 8.3 <https://arxiv.org/abs/1410.3012>`_ (requires a C++ compiler) or `Pythia 6.4.27 <http://arxiv.org/abs/hep-ph/0603175>`_ (requires gfortran)
 * `NLL-fast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_ 1.2 (7 TeV), 2.1 (8 TeV), and 3.1 (13 TeV) (requires a fortran compiler)

These tools need not be installed separately, as the SModelS build system takes care of that. The current default is that both Pythia6 and Pythia8 are installed together with NLLfast.
Finally, the :ref:`database browser <databaseBrowser>` provided by :ref:`smodelsTools.py <smodelsTools>`
requires `IPython <https://ipython.org/>`_, while the :ref:`interactive plotter <interactivePlots>` requires `plotly <https://plot.ly/python/>`_ and `pandas <https://pandas.pydata.org/>`_. 


Installation Methods
^^^^^^^^^^^^^^^^^^^^

.. _phenoInstallation:

 * The first installation method installs SModelS in the source directory.
   After downloading the source from the `SModelS releases page <https://github.com/SModelS/smodels/releases>`_
   and extracting it, run::

     make smodels

   in the top-level directory. The installation will remove redundant folders, install the required 
   dependencies (using pip install) and compile Pythia and NLL-fast. 
   If the cross section computer is not needed, one can replace *smodels* with *smodels_noexternaltools* in the above command.
   In case the Python libraries can not be successfully
   installed, the user can install them separately using his/her preferred method. Pythia and NLL-fast can also be compiled separately
   running **make externaltools**. In case the Fortran comiler isn't found, 
   try *make FC=<path-to-gfortran> smodels* or *make FC=<path-to-gfortran> externaltools*. 

 * If Python's *setuptools* is installed in your machine, SModelS and its dependencies
   can also be installed without the use of pip.
   After downloading the source from the `SModelS releases page <https://github.com/SModelS/smodels/releases>`_
   and extracting it, run::


     setup.py install

   within the main smodels directory. If the python libraries are installed in a system folder (as is the default behavior),
   it will be necessary to run the install command with superuser privilege.
   Alternatively, one can run setup.py with the "--user" flag: ::

     setup.py install --user

   If *setuptools* is not installed, you can try to install the external libraries
   manually and then rerun setup.py.
   For Ubuntu, SL6 machines and other platforms, a recipe is given below.


   Note that this installation method will install smodels into the default system or user directory (e.g. ~/.local/lib/python3/site-packages/).
   Depending on your platform, the environment variables $PATH, $PYTHONPATH, $LD_LIBRARY_PATH
   (or $DYLD_LIBRARY_PATH) might have to be set appropriately.




 * Finally, if *pip3* (or *pip*) is installed in your machine, it is also possible to install SModelS directly without the need for
   downloading the source code: ::

     pip3 install smodels

   in case of system-wide installs or : ::

     pip3 install --user smodels
   
   for user-specific installations.


   Note that this installation method will install smodels into the default system or user directory (e.g. ~/.local/lib/python3/site-packages/).
   Depending on your platform, the environment variables $PATH, $PYTHONPATH, $LD_LIBRARY_PATH
   (or $DYLD_LIBRARY_PATH) might have to be set appropriately.
   Be aware that the example files and the |parameters| discussed in the manual 
   will also be located in your default system or user directory. Furthermore the database
   folder is not included (see :ref:`database installation <installingDB>` below).


There is also a diagnostic tool available: ::

  smodelsTools.py toolbox

should list and check all internal tools (Pythia and NLL-fast) and external
(numpy, scipy, unum, ... ) dependencies.

In case everything fails, please contact smodels-users@lists.oeaw.ac.at


.. _installingDB:

Installing the SModelS Database
-------------------------------

The simplest way is to **provide the URL** of the official :ref:`database <databaseStruct>` as the
database path when running SModelS (see :ref:`path <parameterFilePath>` in |parameters|).
In this case the corresponding database version binary file will be automatically downloaded
and used.  The available database URLs can be found on 
the `SModelS Database releases page <https://github.com/SModelS/smodels-database-release/releases>`_ .
For using the latest official database, which is compatible with the code version used, 
one can also simply set:: 

     path = official

in the |parameters|. Per default, the database pickle file will be located in the users' .cache/smodels/ directory. 
If you want the pickled database file to be cached in a different location, set the environment variable SMODELS_CACHEDIR 
accordingly, e.g. to '/tmp'. 

For performance reasons, from v2.2.0 onwards, the signal regions (SRs) of some of the |EMrs| are aggregated in the official database. 
(For example for CMS-SUS-19-006, the original 174 SRs have been aggregated to 40; this speeds up the calculation 
without too much loss in precision when combining SRs). In order to use the original, non-aggregated |EMrs|, set::

     path = official+nonaggregated


**Alternatively, one can download the text version** of the :ref:`database <databaseStruct>` and pickle it locally. 
This can be convenient if one wants to add or edit experimental results. 
The source code of the available databases can again be found on 
the `SModelS Database releases page <https://github.com/SModelS/smodels-database-release/releases>`_. 
After download, unpack it to a convenient location (e.g., to a 'smodels-database' folder in the SModelS source directory),  
and then specify the local path to this folder in the |parameters|, e.g.:: 

     path = ./smodels-database/

The first time SModelS is run, a :ref:`binary file <databasePickle>` will be built
using this text database folder, which can then be used in all subsequent runs. 
As above, by default this contains some |EMrs| with aggregated SRs. The non-aggregated versions 
are stored as a tarball on the top level of the database folder; for v2.2.0 this is the file *nonaggregated220.tar.gz*. 
To use this, simply expand this tarball in the directory::

 cd <smodels-database folder>
 tar -xzvf nonaggregated220.tar.gz

The database  :ref:`binary file <databasePickle>` will then be re-built accordingly upon first usage. 

The complete list of analyses and results included in the database can be
consulted at `https://smodels.github.io/wiki/ListOfAnalyses <https://smodels.github.io/wiki/ListOfAnalyses>`_.
We note that all the results in the official database release have been
carefully validated  and the validation material can be
found at `https://smodels.github.io/wiki/Validation <https://smodels.github.io/wiki/Validation>`_.

The database can conveniently be updated independently from SModelS code
updates. It suffices to unpack any new database tarball and replace the database
directory or provide the :ref:`path <parameterFilePath>` 
to the new folder, binary or URL address.
In the same fashion, one can easily add additional results as
explained below.


.. _addingFastlim:

Adding FastLim data
^^^^^^^^^^^^^^^^^^^

The official SModelS database can be augmented with data from the
`fastlim <http://cern.ch/fastlim>`_ results. The simplest way is to set:  :: 

     path = official+fastlim

in the |parameters|. It is also possible to
directly download a database binary file including the fastlim maps; dedicated URLs are provied on 
the `SModelS Database releases page <https://github.com/SModelS/smodels-database-release/releases>`_ for this purpose.

For using this with the text database, 
a tarball with the properly converted fastlim-1.0 efficiency maps (*smodels-v1.1-fastlim-1.0.tgz*) is located in the top level directory of the database   
( it can also be downloaded separately from `Github <https://github.com/SModelS/smodels-database-release/blob/master/smodels-v1.1-fastlim-1.0.tgz>`_.)
As for adding non-aggregated results (see above), this tarball simply needs to be exploded to be added to the database: ::

 cd <smodels-database folder>
 tar -xzvf smodels-v1.1-fastlim-1.0.tgz
 rm smodels-v1.1-fastlim-1.0.tgz

Once the fastlim folders have been added to the database,
SModelS auto-detects fastlim results and issues an acknowledgement.

When using the Fastlim results, please properly cite the fastlim paper; for
convenience, a bibtex file is provided in the smodels-fastlim tarball.

Finally we point out that when converting the Fastlim efficiency maps
efficiencies with a relative statistical uncertainty greater than 25%
were set to zero. Also, per default we discard zeroes-only results.


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


System-specific Installation Instructions
-----------------------------------------


Installation on Ubuntu >= 16.04
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Installation on Ubuntu machines should be straightforward with superuser privileges
(if you do not have superuser privileges see instructions below):

 * sudo apt install gfortran python-setuptools python-scipy python-numpy python-docutils python-argparse
 * setup.py install

Note that the last command can be run as superuser, or with the "--user" flag.

Installation on SL7
^^^^^^^^^^^^^^^^^^^

Installation on an SL7 or CentOS7 is straightforward:

 * yum install gcc-c++ scipy numpy

 * pip install unum pyslha argparse


Installation on SL6
^^^^^^^^^^^^^^^^^^^

Installation on an SL6 (Scientific Linux 6 or Scientific Linux CERN 6) machine
is tricky, because SModelS requires a more recent version of *scipy* than is provided by SL6.
We succeeded to install SModelS on SL6 by doing:

 * yum install gcc-c++ libstdc++-devel libevent-devel python-devel lapack lapack-devel blas blas-devel libgfortran python-distutils-extra

followed by:

 * pip install nose unum argparse numpy pyslha scipy

Note, that these steps can safely be done within a Python ``virtualenv``.
Pip can also be called with the "--user" flag.

Installation on other platforms or without superuser privileges using Anaconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

From version 1.1.1 on, SModelS comes with a simple C++ interface, see the cpp directory.
Obviously, a C++ compiler is need, alongside with the python developers
(header) files (libpython-dev on ubuntu, python-devel on rpm-based distros).



