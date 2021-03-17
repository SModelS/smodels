.. index:: What's New

.. |element| replace:: :ref:`element <element>`
.. |elements| replace:: :ref:`elements <element>`
.. |topology| replace:: :ref:`topology <topology>`
.. |topologies| replace:: :ref:`topologies <topology>`
.. |decomposition| replace:: :doc:`decomposition <Decomposition>`
.. |constraint| replace:: :ref:`constraint <ULconstraint>`
.. |constraints| replace:: :ref:`constraints <ULconstraint>`
.. |runSModelS| replace:: :ref:`runSModelS.py <runSModelS>`
.. |database| replace:: :ref:`database <Database>`
.. |Fastlim| replace:: :ref:`Fastlim <addingFastlim>`
.. |output| replace:: :ref:`output <smodelsOutput>`
.. |results| replace:: :ref:`experimental results <ExpResult>`
.. |txnames| replace:: :ref:`txnames <TxName>`
.. |EM| replace:: :ref:`EM-type <EMtype>`
.. |UL| replace:: :ref:`UL-type <ULtype>`
.. |EMr| replace:: :ref:`EM-type result <EMtype>`
.. |ULr| replace:: :ref:`UL-type result <ULtype>`
.. |EMrs| replace:: :ref:`EM-type results <EMtype>`
.. |ULrs| replace:: :ref:`UL-type results <ULtype>`
.. |ExpRes| replace:: :ref:`Experimental Result<ExpResult>`
.. |ExpRess| replace:: :ref:`Experimental Results<ExpResult>`
.. |expres| replace:: :ref:`experimental result<ExpResult>`
.. |express| replace:: :ref:`experimental results<ExpResult>`
.. |Dataset| replace:: :ref:`DataSet<DataSet>`
.. |Datasets| replace:: :ref:`DataSets<DataSet>`
.. |dataset| replace:: :ref:`data set<DataSet>`
.. |datasets| replace:: :ref:`data sets<DataSet>`
.. |parameters| replace:: :ref:`parameters file <parameterFile>`
.. |ssigBRe| replace:: :math:`\sum \sigma \times BR \times \epsilon`
.. |Cpp| replace:: :ref:`C++ Interface<Cpp>`



What's New
==========
The major novelties of all releases since v1.0 are as follows:

New in Version 2.1.0:
^^^^^^^^^^^^^^^^^^^^^
  
  * ability to merge databases using '+' as a delimiter

New in Version 2.0.0:
^^^^^^^^^^^^^^^^^^^^^

  * Introduction of :ref:`particle class <particleClass>`
  * Introduction of model class (see :ref:`Basic Input <basicInput>`)
  * Input model can now be defined by an SLHA file with :ref:`QNUMBERS blocks <qnumberSLHA>`
  * Unified treatment of SLHA and LHE input files (see :ref:`decomposer <decomp>` and :ref:`LHE-reader <lhereader>`)
  * :ref:`Decomposition <decomposition>`  and |ExpRess| can now handle :ref:`lifetime dependent results <widthGrid>`
  * Added :ref:`field "type" <txnameFile>` to the experimental results in the database
  * Added (optional) :ref:`field "intermediateState" <txnameFile>` to the experimental results in the database
  * Inclusive branches can now describe inclusive vertices
  * Added possibility for analysis specific detector size
  * New :ref:`missing topologies <topCoverage>` algorithm and output
  * Added "latest" and "latest_fastlim" :ref:`Database <parameterFileDatabase>` abbreviations
  * Added support for central database server
  * Small bug fix in :ref:`likelihood computation <likelihoodCalc>`
  * Small fix due to an API change in pyhf 0.6
  * Changes in output: :ref:`width values added <pyOut>`, :ref:`coverage groups <coverageGroups>` and others (see :ref:`output description <outputDescription>` for details)
  * Added option for signal strength multipliers in :ref:`cross section calculator <xsecCalc>`
  * Small bug fixes in :ref:`models <basicInput>`  



New in Version 1.2.4:
^^^^^^^^^^^^^^^^^^^^^
  * added pyhf support
  * pickle path bug fix
  * bug fix for parallel xseccomputers
  * Introduced the SMODELS_CACHEDIR environment variable to allow for a different
    location of the cached database file
  * fixed dataId bug in datasets

New in Version 1.2.3:
^^^^^^^^^^^^^^^^^^^^^
  * |database| updated with results from more than 20 new analyses
  * server for databases is now smodels.github.io, not smodels.hephy.at
  * small bug fix for displaced topologies
  * small fix in slha printer, r_expected was r_observed
  * :ref:`Downloaded database files <parameterFilePath>` now stored in $HOME/.cache/smodels

New in Version 1.2.2:
^^^^^^^^^^^^^^^^^^^^^

  * Updated official |database|, added T3GQ eff maps and a few ATLAS 13 TeV results, see `github database release page <https://github.com/SModelS/smodels-database-release/releases>`_
  * Database "official" now refers to a database without fastlim results, "official_fastlim", to the official database *with* fastlim
  * List displaced signatures in :ref:`missing topologies <topCoverage>`
  * Improved description about lifetime reweighting in doc
  * Fix in :ref:`cluster<ULcluster>` for asymmetric masses
  * Small improvements in the :ref:`interactive plots tool <interactivePlots>`

New in Version 1.2.1:
^^^^^^^^^^^^^^^^^^^^^

  * Fix in particleNames.py for non-MSSM models
  * Fixed the `marginalize <marginalize.html>`_ recipe
  * Fixed the T2bbWWoff 44 signal regions plots in `ConfrontPredictions <ConfrontPredictions.html>`_  in manual

New in Version 1.2.0:
^^^^^^^^^^^^^^^^^^^^^

  * Decomposition and experimental results can include
    non-MET BSM final states (e.g. heavy stable charged particles)
  * Added lifetime reweighting at |decomposition| for meta-stable particles
  * Added finalState property for Elements
  * Introduction of :ref:`inclusive simplified models <inclusiveSMS>`
  * Inclusion of HSCP and R-hadron results in the database

New in Version 1.1.3:
^^^^^^^^^^^^^^^^^^^^^

  * Support for :ref:`covariance matrices <combineSRs>` and combination of signal regions (see :ref:`combineSR <parameterFileCombineSRs>` in |parameters|)
  * New plotting tool added to smodelsTools (see :ref:`Interactive Plots Maker <interactivePlots>`)
  * Path to particles.py can now be specified in parameters.ini file (see :ref:`model <parameterFileModel>` in |parameters|)
  * Wildcards allowed when selecting analyses, datasets, txnames (see :ref:`analyses <parameterFileAnalyses>`, :ref:`txnames <parameterFileTxnames>` and :ref:`dataselector <parameterFileDataselector>`  in |parameters|)
  * Option to show individual contribution from topologies to total theory prediction (see :ref:`addTxWeights <parameterFileAddTxWeights>` in |parameters|)
  * URLs are allowed as database paths (see :ref:`path <parameterFilePath>` in |parameters|)
  * Python default changed from python2 to python3
  * Fixed lastUpdate bug, now giving correct date
  * Changes in pickling (e.g. subpickling, removing redundant zeroes)
  * Added fixpermissions to smodelsTools.py, for system-wide installs (see :ref:`Files Permissions Fixer <permissionsFixer>`)
  * Fixed small issue with pair production of even particles
  * Moved the :ref:`code documentation <CodeDocs>` to the manual
  * Added :ref:`option for installing <phenoInstallation>` within the source folder


New in Version 1.1.2:
^^^^^^^^^^^^^^^^^^^^^

* Database update only, the code is the same as v1.1.1

New in Version 1.1.1:
^^^^^^^^^^^^^^^^^^^^^

* |Cpp|
* Support for pythia8 (see :ref:`Cross Section Calculator <xsecCalc>`)
* improved binary database
* automated SLHA and LHE file detection
* Fix and improvements for missing topologies
* Added SLHA-type output
* Small improvements in interpolation and clustering


New in Version 1.1.0:
^^^^^^^^^^^^^^^^^^^^^

* the inclusion of efficiency maps (see |EMrs|)
* a new and more flexible database format (see :ref:`Database structure <databaseStruct>`)
* inclusion of likelihood and :math:`\chi^2` calculation for |EMrs|
  (see :ref:`likelihood calculation <likelihoodCalc>`)
* extended information on the :ref:`topology coverage <topCoverage>`
* inclusion of a database broswer tool for easy access to the information
  stored in the database (see :ref:`database browser <databaseBrowser>`)
* the database now supports also a more efficient :ref:`binary format <databasePickle>`
* performance improvement for the |decomposition| of the input model
* inclusion of new simplified results to the |database| (including a few 13 TeV results)
* |Fastlim| efficiency maps can now also be used in SModelS
