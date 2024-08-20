.. index:: What's New

.. |decomposition| replace:: :doc:`decomposition <Decomposition>`
.. |constraint| replace:: :ref:`constraint <ULconstraint>`
.. |constraints| replace:: :ref:`constraints <ULconstraint>`
.. |runSModelS| replace:: :ref:`runSModelS.py <runSModelS>`
.. |database| replace:: :ref:`database<databaseDefs>`
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


New in Version 3.0.0:
^^^^^^^^^^^^^^^^^^^^^

  * **Extension to general SMS topologies** (no longer restricted to Z\ :sub:`2` symmetric topologies)
  * Large code refactoring
  * Added option for selecting which quantum numbers to be ignored in prompt results (see  the :ref:`ignorePromptQNumbers<erasePrompt>` option in parameters.ini)
  * Added :ref:`outputFormat <parameterFileOutputFormat>` option to parameters.ini to allow for writing the output using the old format (old bracket notation instead of new string representation of SMS)
  * Output for |EMrs| now reports negative log likelihoods, instead of likelihoods
  * Changes in missing topologies (coverage): the grouping of topologies is now done only by final states and ignores the topology structure (see :ref:`missing topologies <finalStateSMS>`)
  * `model.updateParticles <base.html#base.model.Model.updateParticles>`_ and `crossSection.getXsecFromSLHAFile <base.html#base.crossSection.getXsecFromSLHAFile>`_ can now also supply SLHA strings instead of SLHA filenames as argument
  * Z2parity attribute of particles is no longer needed (in :ref:`QNUMBERS blocks <qnumberSLHA>`)
  * `jsonFiles <DatabaseStructure.html#experimental-result-folder>`_ entries in database now allow to specify also pyhf region names and region types (signal or control region)
  * Enabled EMs for control regions to emulate signal leakage to control regions (see :ref:`pyhf Approach <pyhfllhd>`)
  * Introduced signalUncertainty field in the globalInfo.txt files to quantify signal uncertainties for pyhf statistical models (see :ref:`pyhf Approach <pyhfllhd>`)
  * Introduced centralized database dictionary to decrease redundancies in SMS matching (see :ref:`SMS Dictionary <smsDictionary>`)
  * Clustering of SMS for UL results replaced by a (simple) K-means clustering algorithm (see :ref:`Clustering <cluster>`)
  * Changed the lock file mechanism for downloading the database to work with all file systems, fixes `#37 <https://github.com/SModelS/smodels/issues/37>`_
  * Added CITATION.cff file, closes `#38 <https://github.com/SModelS/smodels/issues/38>`_
  * Added smodels-analyses.json in database

  * |database| extension: added results from 11 ATLAS and 6 CMS analyses (hfm=HistFactory model, cov=covariance matrix for SR combination):

     * results from ATLAS: ATLAS-SUSY-2018-33 (EM), ATLAS-SUSY-2018-16 (EM+hfm), ATLAS-SUSY-2018-13 (EM), ATLAS-SUSY-2018-09 (UL), ATLAS-EXOT-2019-03 (UL), ATLAS-EXOT-2018-48 (UL), ATLAS-EXOT-2018-06 (UL), ATLAS-EXOT-2013-11 (UL)
     * results from CMS: CMS-SUS-21-007 (UL), CMS-EXO-20-008 (UL), CMS-EXO-19-012 (UL), CMS-EXO-16-057 (UL), CMS-EXO-12-059 (UL)
     * EM results from recasts: ATLAS-SUSY-2019-08 (MA5), ATLAS-SUSY-2018-42 (`LLP repo <https://github.com/llprecasting/recastingCodes/tree/main/HSCPs/ATLAS-SUSY-2018-42>`_), ATLAS-SUSY-2018-22 (CM2), CMS-EXO-20-004 (`MonoXSMS <https://doi.org/10.5281/zenodo.13324003>`_)


    Note that the 4 ATLAS-EXOT and 5 CMS-EXO analyses above are resonance searches, while ATLAS-SUSY-2018-13 is an RPV SUSY search. These 10 analyses can only be treated with the new graph-based topology description of SModelS v3.


New in Version 2.3.3:
^^^^^^^^^^^^^^^^^^^^^

  * added :ref:`resummino cross section <xsecResummino>` computer
  * fixed bug in computation of error on muhat, for `pyhf likelihoods <tools.html#tools.pyhfInterface.PyhfUpperLimitComputer.lmax>`_
    (affects mostly the numpy backend)
  * small change in initialisation of gradient descent method for computation of
    `combined mu_hat <tools.html#tools.analysesCombinations.AnaCombLikelihoodComputer.lmax>`_, to increate robustness of method

New in Version 2.3.2:
^^^^^^^^^^^^^^^^^^^^^

  * fixed bug in initialisation of :ref:`analyses combination <analysesCombination>`
  * added smodels version to output of "txt" printer

New in Version 2.3.1:
^^^^^^^^^^^^^^^^^^^^^

  * fixed bug for reading :ref:`QNUMBERS blocks <qnumberSLHA>` from SLHA files
  * small fixes in how pythia6 and pythia8 are built
  * small fix in truncated Gaussian llhd experimental feature
  * small fix in computation of combined upper limits
  * combinationmatrices are now forced to be symmetric
  * added isCombinableWith method also for CombinedDataSets
  * added a recipe for how to use a `combinations matrix <combinationsmatrix.html>`_
  * `runtime.nCPUs() <tools.html#tools.runtime.nCPUs>`_ now returns number of available CPUs, not all CPUs
  * :ref:`xsecComputer <xsecCalc>` now has --tempdir option
  * StatsComputer now has `CLs <tools.html#tools.statsTools.StatsComputer.CLs>`_ method
  * changed default prompt width from 1e-8 to 1e-11 GeV in code

New in Version 2.3.0:
^^^^^^^^^^^^^^^^^^^^^

  * fixed bug for an LHE input only with anti-particles
  * fixed error that truncated signal yields when computing expected upper limits
  * added minMass parameter for setting a minimum mass threshold for BSM masses
  * fixed C++ interface to work with python 3.11
  * bumped up pythia8 from 8.307 to 8.308
  * SModelS can now track inter-analyses combinability at the level of whole `analyses <experiment.html#experiment.expResultObj.ExpResult.isCombinableWith>`_ as well as individual `signal regions <experiment.html#experiment.datasetObj.DataSet.isCombinableWith>`_
  * added support for :ref:`SLv2 <SLV2>` (Gaussian with a skew), arXiv:1809.05548
  * refactored the statistics modules
  * introduced "full_llhds" database add-on (see :ref:`parameter.ini file <parameterFileDatabase>`)
  * |database| extension, added new results from 6 ATLAS and 4 CMS analyses (hfm=HistFactory model, cov=covariance matrix for SR combination):

     * results from ATLAS: ATLAS-SUSY-2018-05 (UL,EM+hfm), ATLAS-SUSY-2018-32 (EM+hfm), ATLAS-SUSY-2018-41 (EM+cov, updated), ATLAS-SUSY-2018-42 (UL,EM), ATLAS-SUSY-2019-02 (UL,EM+cov), ATLAS-SUSY-2013-12 (8TeV, EM)
     * results from CMS: CMS-SUS-19-010 (UL), CMS-SUS-20-004 (UL,EM), CMS-SUS-21-002 (UL,EM+cov)
     * added expected ULs to CMS-SUS-19-009

New in Version 2.2.1:
^^^^^^^^^^^^^^^^^^^^^

  * fixes in :ref:`analyses combinations <analysesCombination>`, simplified and pyhf :ref:`likelihoods <likelihoodCalc>`
  * small fixes for python 3.10
  * bumped up pythia8 from 8.306 to 8.307
  * included :ref:`example on how to plot likelihoods from analysis combination  <Examples>`
  * small bug fix for particle addition

New in Version 2.2.0.post1:
^^^^^^^^^^^^^^^^^^^^^^^^^^^

  * removed dependency on importlib.metadata to make it work with python <= 3.7

New in Version 2.2.0:
^^^^^^^^^^^^^^^^^^^^^

  * introduced (user-defined) :ref:`combinations of analyses <analysesCombination>`
  * changed expected limits computed with pyhf from post-fit to pre-fit
  * a few smaller changes around expected likelihoods and limits
  * changed default value of :ref:`promptWidth parameter <parameterFileModel>` from 1e-8 to 1e-11 GeV
  * allow :ref:`ncpus <parameterFileNcpus>` to take on zero and negative values in
    ini file [meaning use all but this (absolute) number of CPU cores]
  * notion of "nonaggregated" databases introduced
  * small fixes in the :ref:`Howto's <Examples>`
  * updates in references.bib, installation notes
  * more small fixes in unit tests
  * |database| extension, added new results from 4 ATLAS and 11 CMS analyses:

     * results from ATLAS: ATLAS-SUSY-2018-08 (UL+EM), ATLAS-SUSY-2018-40 (UL+EM), ATLAS-SUSY-2018-41 (UL+EM), ATLAS-SUSY-2019-09 (UL+EM, full likelihood)
     * results from CMS: CMS-SUS-16-050 (EM), CMS-SUS-18-004 (UL), CMS-SUS-18-007 (UL), CMS-SUS-19-008 (UL), CMS-SUS-19-011 (UL), CMS-SUS-19-013 (UL), CMS-SUS-20-001 (UL), CMS-SUS-20-002 (UL)
     * recast with MadAnalysis5: CMS-SUS-16-039 (EM), CMS-SUS-16-048 (EM), CMS-SUS-19-006 (EM); all incl. covariance matrices

New in Version 2.1.1:
^^^^^^^^^^^^^^^^^^^^^

  * caching weight matrix in simplified :ref:`likelihoods <likelihoodCalc>`
  * notion of "debug" databases introduced
  * introduced :ref:`reportAllSRs <parameterFileReportAllSRs>` option
  * tiny fix in mybinder link (see https://pypi.org/project/smodels/)
  * small fixes in unit tests
  * improved truncated Gaussians in likelihoodsFromLimits (but kept as experimental feature)
  * :ref:`experimental features <parameterExperimentalFeatures>` can now be turned on via ini file

New in Version 2.1.0:
^^^^^^^^^^^^^^^^^^^^^

  * Ability to merge :ref:`Databases <parameterFileDatabase>` using '+' as a delimiter: "latest_fastlim" and "official_fastlim" are now written as "latest+fastlim", and "official+fastlim".
  * useSuperseded flag in `getExpResults  <experiment.html#experiment.databaseObj.Database.getExpResults>`_ is marked as deprecated, as we now just put superseded results in separate database
  * |Datasets| now have an `.isCombinableWith <experiment.html#experiment.datasetObj.DataSet.isCombinableWith>`_ function
  * Slightly extended output of :ref:`summary printer <parameterFileSummaryprinter>`
  * Added scan summary (:ref:`summary.txt <scanSummary>`) when running over multiple files
  * Added :ref:`expandedOutput <parameterFileSLHAprinter>` option to slha-printer
  * :ref:`Output <outputDescription>` for efficiency-map results now reports :ref:`L, L_max and L_SM <likelihoodCalc>`
  * The :ref:`likelihood <likelihoodCalc>` is now maximized only for positive values of the signal strength
    in the computation of L_max
  * Pythia8 version in :ref:`xsecComputer <xsecCalc>` updated from 8226 to 8306
  * Improved :ref:`interactive plots <interactivePlots>`
  * |database| updated with results from 5 new ATLAS and 1 new CMS analyses: CMS-EXO-19-010 (disappearing tracks) UL,  ATLAS-SUSY-2016-08 (displaced leptons) EM,      ATLAS-SUSY-2018-10 (1l+jets) UL+EM, ATLAS-SUSY-2018-12 (0l+jets) UL+EM, ATLAS-SUSY-2018-22 (0l+jets) UL+EM, ATLAS-SUSY-2018-23 (EWino, WH) UL
  * added EM results for  ATLAS-SUSY-2017-03 (EWino, WZ), ATLAS-SUSY-2018-06 (EWino, WZ),  ATLAS-SUSY-2018-14 (sleptons),  CMS-SUSY-14-021 (stops)
  * created and added THSCPM10 and THSCPM11 EMs for ATLAS-SUSY-2016-32;
  * replaced some 8 TeV ATLAS conf notes with the published results:   (ATLAS-CONF-2013-007 -> ATLAS-SUSY-2013-09, ATLAS-CONF-2013-061 -> ATLAS-SUSY-2013-18, ATLAS-CONF-2013-089  -> ATLAS-SUSY-2013-20)
  * corrected off-shell regions of some existing |EMrs| (in three 13 TeV and eigth 8 TeV analyses).


New in Version 2.0.0:
^^^^^^^^^^^^^^^^^^^^^

  * Introduction of :ref:`particle class <particleClass>`
  * Introduction of model class (see :ref:`Basic Input <basicInput>`)
  * Input model can now be defined by an SLHA file with :ref:`QNUMBERS blocks <qnumberSLHA>`
  * Unified treatment of SLHA and LHE input files (see :ref:`decomposer <decomp>` and :ref:`LHE-reader <lhereader>`)
  * :ref:`Decomposition <decomposition>`  and |ExpRess| can now handle lifetime dependent results
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
  * Fix in cluster for asymmetric masses
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
