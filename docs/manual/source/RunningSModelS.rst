.. index:: Using SModelS

.. |invisible compression| replace:: :ref:`invisible compression <invComp>`
.. |mass compression| replace:: :ref:`mass compression <massComp>`
.. |element| replace:: :ref:`element <element>`
.. |elements| replace:: :ref:`elements <element>`
.. |topology| replace:: :ref:`topology <topology>`
.. |topologies| replace:: :ref:`topologies <topology>`
.. |decomposition| replace:: :doc:`decomposition <Decomposition>`
.. |Decompose| replace:: :doc:`Decompose <Decomposition>`
.. |theory predictions| replace:: :doc:`theory predictions <TheoryPredictions>`
.. |theory prediction| replace:: :doc:`theory prediction <TheoryPredictions>`
.. |constraint| replace:: :ref:`constraint <ULconstraint>`
.. |constraints| replace:: :ref:`constraints <ULconstraint>`
.. |runSModelS| replace:: :ref:`runSModelS.py <runSModelS>`
.. |database| replace:: :ref:`database <Database>`
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


.. _runningSModelS:

Using SModelS
=============

SModelS can take SLHA or LHE files as input (see :doc:`Basic Input <BasicInput>`).
It ships with a command-line tool :ref:`runSModelS.py <runSModelS>`, which reports
on the SMS |decomposition| and |theory predictions| in several :ref:`output formats <smodelsOutput>`.

For users more familiar with Python and the SModelS basics, an example
code :ref:`Example.py <exampleCode>` is provided showing how to access
the main SModelS functionalities: :ref:`decomposition <Decomposition>`, the |database|
and :ref:`computation of theory predictions <TheoryPredictions>`.


The command-line tool (:ref:`runSModelS.py <runSModelS>`) and the example Python
code (:ref:`Example.py <exampleCode>`) are described below.


.. note:: For non-MSSM (incl. non-SUSY) input models the user needs to write their own *model.py*
          file and specify which BSM particles are even or odd under the assumed
          Z\ :sub:`2` symmetry (see :ref:`adding new particles <newParticles>`).
          From version 1.2.0 onwards it is also necessary to define the BSM particle quantum numbers in the same file [#]_.





.. _runSModelS:

runSModelS.py
-------------


*runSModelS.py* covers several different applications of the SModelS functionality,
with the option of turning various features on or off, as well as
setting the :ref:`basic parameters <parameterFile>`.
These functionalities include detailed checks of input SLHA files,
running the |decomposition|,
evaluating the :doc:`theory predictions <TheoryPredictions>` and comparing them to the experimental
limits available in the |database|,
determining :ref:`missing topologies <topCoverage>` and printing the |output|
in several available formats.

Starting on v1.1, *runSModelS.py* is equipped with two additional
functionalities. First, it can process a folder containing a set of SLHA or LHE
file, second, it supports parallelization of this input folder.



**The usage of runSModelS is:**

.. include:: RunSModelS.rst


A typical usage example is: ::

   runSModelS.py -f inputFiles/slha/simplyGluino.slha -p parameters.ini -o ./ -v warning

The resulting |output| will be generated in the current folder, according to the printer options set in the
|parameters|.



.. _parameterFile:


The Parameters File
-------------------

The basic options and parameters used by *runSModelS.py* are defined in the parameters file.
An example parameter file, including all available parameters together
with a short description, is stored in :download:`parameters.ini <images/parameters.ini>`.
If no parameter file is specified, the default parameters stored in
:download:`smodels/etc/parameters_default.ini <images/parameters_default.ini>` are used.
Below we give more detailed information about each entry in the parameters file.


.. _parameterFileOptions:

* *options*: main options for turning SModelS features on or off

.. _parameterFileCheckInput:

  * **checkInput** (True/False): if True, *runSModelS.py* will run the :ref:`file check tool <fileChecks>` on the input file and verify if the input contains all the necessary information.

.. _parameterFileDoInvisible:

  * **doInvisible** (True/False): turns |invisible compression| on or off during the |decomposition|.

.. _parameterFileDoCompress:

  * **doCompress** (True/False): turns |mass compression| on or off during the |decomposition|.

.. _parameterFileComputeStatistics:

  * **computeStatistics** (True/False): turns the likelihood and :math:`\chi^2` computation on or off
    (see :ref:`likelihood calculation <likelihoodCalc>`).
    If True, the likelihood and :math:`\chi^2` values are computed for the |EMrs|.

.. _parameterFileTestCoverage:

  * **testCoverage** (True/False): set to True to run the :ref:`coverage <topCoverage>` tool.

.. _parameterFileCombineSRs:

  * **combineSRs** (True/False): set to True to use, whenever available, covariance matrices to combine signal regions. NB this might take a few secs per point. Set to False to use only the most sensitive signal region (faster!). Available v1.1.3 onwards.

.. _parameterFileParticles:

* *particles*: defines the particle content of the BSM model

.. _parameterFileModel:

  * **model**: pathname to the Python file that defines the particle content of the BSM model, given either in Unix file notation ("/path/to/model.py") or as Python module path ("path.to.model"). Defaults to *share.models.mssm* which is a standard MSSM. See smodels/share/models folder for more examples. Directory name can be omitted; in that case, the current working directory as well as smodels/share/models are searched for.

.. _promptWidth:

  * **promptWidth**: total decay width in GeV above which decays are considered prompt; available v2.0 onwards.

.. _stableWidth:

  * **stableWidth**: total decay width in GeV below which particles are considered as (quasi)stable; available v2.0 onwards.


.. _parameterFileParameters:

* *parameters*: basic parameter values for running SModelS

.. _parameterFileSigmacut:

  * **sigmacut** (float): minimum value for an |element| weight (in fb). :ref:`Elements <element>`
    with a weight below sigmacut are neglected during the |decomposition|
    of SLHA files (see :ref:`Minimum Decomposition Weight <minweight>`).
    The default value is 0.03 fb. Note that, depending on the input model, the running time may increase considerably if sigmacut is too low, while too large values might eliminate relevant |elements|.

.. _parameterFileMinmassgap:

  * **minmassgap** (float): maximum value of the mass difference (in GeV) for
    perfoming :ref:`mass compression <massComp>`. *Only used if doCompress = True*

.. _parameterFileMaxcond:

  * **maxcond** (float): maximum allowed value (in the [0,1] interval) for the violation of :ref:`upper limit conditions <ULconditions>`. A zero value means the conditions are strictly enforced, while 1 means the conditions are never enforced.
    *Only relevant for printing the* :ref:`output summary <fileOut>`.

.. _parameterFileNcpus:

  * **ncpus** (int): number of CPUs. When processing multiple SLHA/LHE files,
    SModelS can run in a parallelized fashion, splitting up the input files in equal chunks.
    *ncpus = -1* parallelizes to as many processes as number of CPU cores of the machine. Default value is 1. Warning: python already parallelizes many tasks internally.

.. _parameterFileDatabase:

* *database*: allows for selection of a subset of :ref:`experimental results <ExpResult>` from the |database|

.. _parameterFilePath:

  * **path**: the absolute (or relative) path to the :ref:`database <databaseStruct>`. The user can supply either the directory name of the database, or the path to the :ref:`pickle file <databasePickle>`. Also http addresses may be given, e.g. https://smodels.github.io/database/official113. The path "official" refers to the official database of your SModelS version -- without fastlim; "official_fastlim" includes fastlim results. In addition, the paths "latest_fastlim" and "latest" refer to the latest databases, with and without fastlim results, respectively. See the `github database release page <https://github.com/SModelS/smodels-database-release/releases>`_ for a list of public database versions.

.. _parameterFileAnalyses:

  * **analyses** (list of results): set to ['all'] to use all available results. If a list of :ref:`experimental analyses <ExpResult>`
    is given, only these will be used. For instance, setting analyses = CMS-PAS-SUS-13-008,ATLAS-CONF-2013-024
    will only use the |results| from `CMS-PAS-SUS-13-008 <https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13008>`_
    and `ATLAS-CONF-2013-024 <https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-024/>`_.
    Wildcards (*, ?, [<list-of-or'ed-letters>]) are expanded in the same way the shell does wildcard expansion for file names.
    So analyses = CMS* leads to evaluation of results from the CMS-experiment only, for example. *SUS* selects everything containining SUS, no matter if from CMS or ATLAS.
    Furthermore selection of analyses can be confined on their centre-of-mass energy with a suffix beginning with a colon and an energy string in unum-style, like :13*TeV. Note that the asterisk behind the colon is not a wildcard. :13, :13TeV and :13 TeV are also understood but discouraged.

.. _parameterFileTxnames:

  * **txnames** (list of topologies): set to ['all'] to use all available simplified model |topologies|. The |topologies| are labeled according to the :ref:`txname convention <TxName>`.
    If a list of |txnames| are given, only the corresponding |topologies| will be considered. For instance, setting txnames = T2 will
    only consider |results| for :math:`pp \to \tilde{q} + \tilde{q} \to  (jet+\tilde{\chi}_1^0) + (jet+\tilde{\chi}_1^0)`
    and the |output| will only contain constraints for this topology.
    *A list of all* |topologies| *and their corresponding* |txnames| *can be found* `here <https://smodels.github.io/wiki/SmsDictionary>`_
    Wildcards (\*, ?, [<list-of-or'ed-letters>]) are expanded in the same way the shell does wildcard expansion for file names.
    So, for example, txnames = T[12]*bb* picks all txnames beginning with T1 or T2 and containg bb as of the time of writing were: T1bbbb, T1bbbt, T1bbqq, T1bbtt, T2bb, T2bbWW, T2bbWWoff

.. _parameterFileDataselector:

  * **dataselector** (list of datasets): set to ['all'] to use all available |datasets|. If dataselector = upperLimit (efficiencyMap), only |ULrs| (|EMrs|) will be used. Furthermore, if
    a list of signal regions (|datasets|) is given, only the |results| containing these datasets will be used. For instance, if dataselector = SRA mCT150,SRA mCT200, only
    these signal regions will be used.
    Wildcards (\*, ?, [<list-of-or'ed-letters>]) are expanded in the same way the shell does wildcard expansion for file names. Wildcard examples are given above.

.. _parameterFileDataTypes:

  * **dataTypes** dataType of the analysis (all, efficiencyMap or upperLimit). Can be wildcarded with usual shell wildcards: * ? [<list-of-or'ed-letters>]. Wildcard examples are given above.

.. _parameterFilePrinter:

* *printer*: main options for the |output| format

.. _parameterFileOutputType:

  * **outputType** (list of outputs): use to list all the output formats to be generated.
    Available output formats are: summary, stdout, log, python, xml, slha.

.. _parameterFileStdoutprinter:

* *stdout-printer*: options for the stdout or log printer

.. _parameterFilePrintDatabase:

  * **printDatabase** (True/False): set to True to print the list of selected |results| to stdout.

.. _parameterFileAddAnaInfo:

  * **addAnaInfo**  (True/False): set to True to include detailed information about the |txnames| tested by each :ref:`experimental result <ExpResult>`. *Only used if printDatabase=True*.

.. _parameterFilePrintDecomp:

  * **printDecomp** (True/False): set to True to print basic information from the |decomposition| (|topologies|, total weights, ...).

.. _parameterFileAddElementInfo:

  * **addElementInfo**  (True/False): set to True to include detailed information about the |elements| generated by the |decomposition|. *Only used if printDecomp=True*.

.. _parameterFilePrintExtendedResults:

  * **printExtendedResults** (True/False): set to True to print extended information about the  |theory predictions|, including the PIDs of the particles
    contributing to the predicted cross section, their masses and the expected upper limit (if available).

.. _parameterFileAddCoverageID:

  * **addCoverageID** (True/False): set to True to print the list of element IDs contributing to each missing topology (see :ref:`coverage <topCoverage>`).
    *Only used if testCoverage = True*. This option should be used along with *addElementInfo = True* so the user can precisely identify
    which elements were classified as missing.

.. _parameterFileSummaryprinter:

* *summary-printer*: options for the summary printer

.. _parameterFileExpandedSummary:

  * **expandedSummary** (True/False): set True to include in the summary output all applicable |results|, False for only the strongest one.

.. _parameterFilePythonprinter:

* *python-printer*: options for the Python printer

.. _parameterFileAddElementList:

  * **addElementList** (True/False): set True to include in the Python output all information about all |elements| generated in the |decomposition|. If set to True the
    output file can be quite large.

.. _parameterFileAddTxWeights:

  * **addTxWeights** (True/False): set True to print the contribution from individual topologies to each theory prediction. Available v1.1.3 onwards.

.. _parameterFileXmlprinter:

* *xml-printer*: options for the xml printer

.. _parameterFileAddElementListXML:

  * **addElementList** (True/False): set True to include in the xml output all information about all |elements| generated in the |decomposition|. If set to True the
    output file can be quite large.

.. _parameterFileAddTxWeightsXML:

  * **addTxWeights** (True/False): set True to print the contribution from individual topologies to each theory prediction. Available v1.1.3 onwards.




.. _smodelsOutput:

The Output
----------

The results of |runSModelS| are printed to the format(s) specified by the **outputType** in the |parameters|.
The following formats are available:

 * a human-readable :ref:`screen output (stdout) <screenOut>` or :ref:`log output <logOut>`. These are intended to
   provide detailed information about the |database|, the |decomposition|, the
   |theory predictions| and the :ref:`missing topologies <topCoverage>`. The
   output complexity can be controlled through several options in the |parameters|. Due to its size, this output
   is not suitable for storing the results from a large scan, being more appropriate for a single file input.

 * a human-readable text file output containing a :ref:`summary of the output <fileOut>`. This format
   contains the main SModelS results: the |theory predictions| and the :ref:`missing topologies <topCoverage>`.
   It can be used for a large scan, since the output can be made quite compact, using the options in the |parameters|.

 * a :ref:`python dictionary <pyOut>` printed to a file containing information about the |decomposition|, the
   |theory predictions| and the :ref:`missing topologies <topCoverage>`. The output can be significantly long, if
   all options in the |parameters| are set to True. However this output can be easily imported to a Python enviroment, making it
   easy to access the desired information. For users familiar with the Python language this is the recommended
   format.

 * a :ref:`xml file <pyOut>` containing information about the |decomposition|, the
   |theory predictions| and the :ref:`missing topologies <topCoverage>`. The output can be significantly long, if
   all options are set to True. Due to its broad usage, the xml output can be easily converted to the
   user's preferred format.

 * a :ref:`SLHA file <slhaOut>` containing information about the
   |theory predictions| and the :ref:`missing topologies <topCoverage>`. The output follows a SLHA-type
   format and contains a summary of the most constraining results and the missed topologies.

A detailed explanation of the information contained in each type of output is given
in :ref:`SModels Output <outputDescription>`.


.. _exampleCode:

Example.py
----------

Although :ref:`runSModelS.py <runSModelS>` provides the main SModelS features with a command line interface,
users more familiar with Python and the SModelS language may prefer to write their own main program.
A simple example code for this purpose is provided in :download:`examples/Example.py`.
Below we go step-by-step through this example code:


* *Import the SModelS modules and methods*. If the example code file is not located in
  the smodels installation folder, simply add "sys.path.append(<smodels installation path>)" before importing smodels. Set SModelS verbosity level.

.. literalinclude:: /examples/Example.py
   :lines: 11-20

* *Set the path to the database URL*. Specify which :ref:`database <databaseStruct>` to use. It can be the path
  to the smodels-database folder, the path to a :ref:`pickle file <databasePickle>` or (starting with v1.1.3) a URL path.

.. literalinclude:: /examples/Example.py
   :lines: 21-22

* *Define the input model*. By default SModelS assumes the MSSM particle content. For using SModelS
  with a different particle content, the user must define the new particle content and set modelFile
  to the path of the model file (see **particles:model** in :ref:`Parameter File <parameterFile>`).

.. literalinclude:: /examples/Example.py
   :lines: 29-30


* *Path to the input file*. Specify the location of the input file. It must be a
  SLHA or LHE file (see :ref:`Basic Input <BasicInput>`).

.. literalinclude:: /examples/Example.py
   :lines: 33-34

* *Set main options for* |decomposition|.
  Specify the values of :ref:`sigmacut <minweight>` and :ref:`minmassgap <massComp>`:

.. literalinclude:: /examples/Example.py
   :lines: 37-38

* |Decompose| *model*. Depending on the type
  of input format, choose either
  the `slhaDecomposer.decompose <theory.html#theory.slhaDecomposer.decompose>`_ or
  `lheDecomposer.decompose <theory.html#theory.slhaDecomposer.decompose>`_ method. The **doCompress** and **doInvisible** options turn the |mass compression| and |invisible compression| on/off.

.. literalinclude:: /examples/Example.py
   :lines: 41-45

* *Access basic information* from decomposition, using the
  `topology list <theory.html#theory.topology.TopologyList>`_
  and `topology  <theory.html#theory.topology.Topology>`_ objects:

.. literalinclude:: /examples/Example.py
   :lines: 47-61

*output:*

.. literalinclude:: /images/ExampleOutput.txt
   :lines: 2-9


* *Load the* |express| to be used to constrain the input model.
  Here, all results are used:

.. literalinclude:: /examples/Example.py
   :lines: 65

Alternatively, the `getExpResults  <experiment.html#experiment.databaseObj.Database.getExpResults>`_ method
can take as arguments specific results to be loaded.

* *Print basic information about the results loaded*.
  Below we show how to count the number of |ULrs| and |EMrs| loaded:

.. literalinclude:: /examples/Example.py
   :lines: 69-76

*output:*

.. literalinclude:: /images/ExampleOutput.txt
   :lines: 11


* *Compute the* |theory predictions| for each |expres|.
  The output is a list of
  `theory prediction objects <theory.html#theory.theoryPrediction.TheoryPrediction>`_
  (for each |expres|):

.. literalinclude:: /examples/Example.py
   :lines: 82-83

* *Print the results*. For each |expres|, loop over the corresponding |theory predictions|
  and print the relevant information:

.. literalinclude:: /examples/Example.py
   :lines: 86-98

*output:*

.. literalinclude:: /images/ExampleOutput.txt
   :lines: 15-22

* *Get the corresponding upper limit*. This value can
  be compared to the |theory prediction| to decide whether a model is excluded or not:

.. literalinclude:: /examples/Example.py
   :lines: 101-102

*output:*

.. literalinclude:: /images/ExampleOutput.txt
   :lines: 23

* *Print the r-value*, i.e. the ratio |theory prediction|/upper limit.
  A value of :math:`r \geq 1` means that an experimental result excludes the input model.
  For |EMrs| also compute the :math:`\chi^2` and :ref:`likelihood <likelihoodCalc>`.
  Determine the most constraining result:

.. literalinclude:: /examples/Example.py
   :lines: 105-113

*output:*

.. literalinclude:: /images/ExampleOutput.txt
   :lines: 24-25

* *Print the most constraining experimental result*. Using the largest *r*-value,
  determine if the model has been excluded or not by the selected |express|:

.. literalinclude:: /examples/Example.py
   :lines: 116-120


*output:*

.. literalinclude:: /images/ExampleOutput.txt
   :lines: 588-589


* *Identify missing topologies*. Using the output from decomposition, identify
  the :ref:`missing topologies <topCoverage>` and print some basic information:

.. literalinclude:: /examples/Example.py
   :lines: 125-145


*output:*

.. literalinclude:: /images/ExampleOutput.txt
   :lines: 591-602,612-617


It is worth noting that SModelS does not include any statistical treatment for
the results, for instance, correction factors like the "look elsewhere effect".
Due to this, the results are claimed to be "likely excluded" in the output.


**Notes:**
 * For an SLHA :ref:`input file <BasicInput>`, the decays of SM :ref:`particles <particleClass>`
   (or BSM Z\ :sub:`2`-even particles) are always ignored during
   the decomposition. Furthermore, if there are two cross sections at different
   calculation order (say LO and NLO) for the same process, only the highest order is used.
 * The list of |elements| can be extremely long. Try setting **addElementInfo** = False
   and/or **printDecomp** = False to obtain a smaller output.
 * A comment of caution is in order regarding naively using the highest :math:`r`-value
   reported by SModelS, as this does not necessarily come from the most sensitive analysis.
   For a rigorous statistical interpretation, one should use the  :math:`r`-value of
   the result with the highest *expected* :math:`r` (:math:`r_{exp}`).
   Unfortunately, for |ULrs|, the expected limits are often not available;
   :math:`r_{exp}` is then reported as N/A in the SModelS output.

.. [#] We note that SLHA files including decay tables and cross sections, together with the corresponding *model.py*, can conveniently be generated via the SModelS-micrOMEGAS interface, see `arXiv:1606.03834 <http://www.arXiv.org/abs/1606.03834>`_
