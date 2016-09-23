.. index:: Running SModelS

.. |invisible compression| replace:: :ref:`invisible compression <invComp>`
.. |mass compression| replace:: :ref:`mass compression <massComp>`
.. |element| replace:: :ref:`element <element>`
.. |elements| replace:: :ref:`elements <element>`
.. |topology| replace:: :ref:`topology <topology>`
.. |topologies| replace:: :ref:`topologies <topology>`
.. |analysis| replace:: :ref:`analysis <ULanalysis>`
.. |analyses| replace:: :ref:`analyses <ULanalysis>`
.. |decomposition| replace:: :doc:`decomposition <Decomposition>`
.. |theory predictions| replace:: :doc:`theory predictions <TheoryPredictions>`
.. |theory prediction| replace:: :doc:`theory prediction <TheoryPredictions>`
.. |constraint| replace:: :ref:`constraint <ULconstraint>`
.. |constraints| replace:: :ref:`constraints <ULconstraint>`
.. |runSModelS| replace:: :ref:`runSModelS.py <runSModelS>`

.. _runningSModelS:

Running SModelS
===============

For the first-time user, SModelS ships with a command-line tool :ref:`runSModelS.py <runSModelS>`, which
takes an SLHA or LHE file as input (see :doc:`Basic Input <BasicInput>`), and reports on the SMS
|decomposition| and |theory predictions| in a simple :ref:`output text file <fileOut>`.

For users more familiar with Python and the SModelS basics, an example
code :ref:`Example.py <exampleCode>` is provided showing how to access
the main SModelS functionalities: :doc:`decomposition <Decomposition>`, :doc:`analysis database <Database>`
and :doc:`computation of theory predictions <TheoryPredictions>`.

The commandline tool (:ref:`runSModelS.py <runSModelS>`) and the example Python
code (:ref:`Example.py <exampleCode>`) are described below.

.. _runSModelS:

runSModelS.py
-------------


*runSModelS.py* covers several different applications of the SModelS functionality,
with the option of turning various features on or off, as well as
setting the :ref:`basic parameters <parameterFile>`.

These functionalities include detailed checks of input SLHA files,
running the |decomposition| and printing the :ref:`output <output>`,
evaluating the :doc:`theory predictions <TheoryPredictions>` and comparing them to the experimental
limits available in the :doc:`database <Database>`,
determining :ref:`missing topologies <topCoverage>` and printing a :ref:`summary text file <output>`.

These settings may be changed from their default values using the :ref:`parameter file <parameterFile>`.


**usage:** 
		runSModelS.py [-h] -f FILENAME [-p PARAMETERFILE] [-o OUTPUTDIR] [-d]
                     [-t] [-V] [-c] [-v VERBOSE] [-T TIMEOUT]

optional arguments:
  -h, --help            show this help message and exit
  -f FILENAME, --filename FILENAME
                        name of SLHA or LHE input file, necessary input, if
                        directory is given, loop over all files in the directory
  -p PARAMETERFILE, --parameterFile PARAMETERFILE
                        name of parameter file, optional argument, if not set,
                        use all parameters from etc/parameters_default.ini
  -o OUTPUTDIR, --outputDir OUTPUTDIR
                        name of output directory, optional argument, default
                        is: results
  -d, --development     enable development output
  -t, --force_txt       force loading the text database
  -V, --version         show program's version number and exit
  -c, --run-crashreport
                        parse crash report file and use its contents for a
                        SModelS run.Supply the crash file simply via '--
                        filename myfile.crash'
  -v VERBOSE, --verbose VERBOSE
                        verbosity level. accepted values are: debug, info,
                        warning, error.
  -T TIMEOUT, --timeout TIMEOUT
                        define a limit on the running time (in secs).If not

.. _parameterFile:


The Parameters File
^^^^^^^^^^^^^^^^^^^

The basic options and parameters used by *runSModelS.py* are defined in the parameter file.
An example parameter file, including all available parameters together
with a short description, is stored in :download:`parameters.ini <images/parameters.ini>`.
If no parameter file is specified the default parameters stored in 
:download:`/etc/parameters_default.ini <images/parameters_default.ini>` are used.
Below we give more detailed information about each entry in the parameters file.

* *path*: relevant folder paths
   * **databasePath** (path to database): the absolute (or relative) path to the :doc:`SModelS database <Database>`, can supply either the directory name of the database, or the file name of a pickle file (see :doc:`Database of Experimental Results <DatabaseStructure>`)

*default values*:

.. literalinclude:: /images/parameters_default.ini
   :lines: 1-2
   
* *options*: main options for turning SModelS features on and off

  * **inputType** (SLHA/LHE): determines the type of input file (see :doc:`Basic Input <BasicInput>`). 
    Must be SLHA for a SLHA input file or LHE for a LHE input file.
  * **checkInput** (True/False): if True, *runSModelS.py* will run several :ref:`consistency checks <fileChecks>` on the input file.
    Makes sure that the file contains all the necessary information. For a SLHA file input :ref:`further checks <slhaChecks>`
    are performed, such as if the file contain charged stable particles, if there are inconsistent decays,...
  * **doInvisible** (True/False): turns |invisible compression| on and off during the |decomposition|.
    Set to False to turn |invisible compression| off.
  * **doCompress** (True/False): turns |mass compression| on and off during the |decomposition|.
    Set to False to turn |mass compression| off.
  * **findMissingTopos** (True/False): set to True to run the :ref:`missing topologies <topCoverage>` tool.

*default values*:

.. literalinclude:: /images/parameters_default.ini
   :lines: 3-8
   
* *parameters*: basic parameter values for running SModelS

  * **sigmacut** (float): minimum value for the |element| weight (in fb). During |decomposition| 
    all elements with weights below sigmacut are neglected. Too Small values of 
    sigmacut (see :ref:`Minimum Decomposition Weight <minweight>`) will result in longer running time, while too large values might eliminate relevant |elements|.  
  * **minmassgap** (float): maximum mass difference value (in GeV) for perfoming :ref:`mass compression <massComp>`.
    *Only used if doCompress = True*
  * **maxcond** (float): maximum allowed value (in the [0,1] interval) for the violation of :ref:`analysis conditions <ULconditions>`.
    A zero value means the conditions are strictly enforced, while 1 means the conditions
    are never enforced. 
    *Only relevant for printing the* :ref:`output summary <fileOut>`.
  * **ncpus** (int): number of CPUs. When processing multiple SLHA/LHE files, SModelS can run in a parallelized fashion, splitting up the input files in equal chunks. "-1" is equal to the number of CPU cores of the machine.

*default values*:

.. literalinclude:: /images/parameters_default.ini
   :lines: 9-13  
   
* *database*: select a subset of the available :doc:`database analyses <Database>`

  * **analyses** (list of analyses name): set to all to use all available analyses. If a list of analyses names
    are given, only these analyses will be applied. For instance, setting analyses = CMS-PAS-SUS-13-008, ATLAS-CONF-2013-024
    will only use the |analyses| from `CMS-PAS-SUS-13-008 <https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13008>`_
    and `ATLAS-CONF-2013-024 <https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-024/>`_.
  * **txnames** (list of constraints): set to all to use all available :ref:`analyses constraints <ULconstraint>`.
    If a list of |constraints| are given (using the :ref:`Tx name shorthand <Txname>`,
    only these |constraints| will be considered. For instance, setting txnames = T2 will
    only consider the analyses containing upper limits for :math:`[[[jet]],[[jet]]]`.
    A list of all |constraints| and their :ref:`Tx names <Txname>` can be found `here <http://smodels.hephy.at/wiki/SmsDictionary>`_)

*default values*:

.. literalinclude:: /images/parameters_default.ini
   :lines: 14-16  
   
   
* *stdout*: basic options for the :ref:`screen output <screenOut>`

  * **printDecomp** (True/False): set to True to print basic information from the |decomposition| (|topologies|, total weights, ...)
  * **addElmentInfo** (True/False): set to True to include detailed information about the |elements| generated by the |decomposition|
    *Only used if printDecomp=True*.
  * **printAnalyses** (True/False): set to True to print basic information about the analyses being used 
    (list of :doc:`analysis names <AnalysesNames>`, luminosities,...)
  * **addAnaInfo** (True/False): set to True to include detailed information about the |elements| tested by each |analysis| (|elements|
    appearing in the :ref:`analysis constraint <ULconstraint>`). *Only used if printAnalyses=True*.
  * **printResults** (True/False): set to True to print information about the |theory predictions| computed and the
    respective analysis upper limits

*default values*:

.. literalinclude:: /images/parameters_default.ini
   :lines: 17-22
   
* *file*: basic options for the :ref:`file output <fileOut>`

  * **expandedSummary** (True/False): set to True to print to file all applicable analyses. Set to False to print only
    the most constraining analysis.
  * **addConstraintInfo** (True/False): set to True to include the analysis |constraint| in :ref:`bracket notation <BracketNotation>`
    for each analysis printed to file.

*default values*:

.. literalinclude:: /images/parameters_default.ini
   :lines: 23-25             





.. _output:

The Output
^^^^^^^^^^

The results of |runSModelS| are printed both to the screen and to the output (summary) file. If no
output file is specified when calling |runSModelS|, the file output will be printed to ./summary.txt.
Below we explain in detail the information contained in the :ref:`screen output <screenOut>` and
the :ref:`output file <fileOut>`.  

.. _screenOut:

Screen Output
*************

If all the options in *stdout* (**printDecomp**,  **addElmentInfo**, **printAnalyses**, **addAnaInfo** and
**printResults**)  are set to True (see :ref:`parameter file <parameterFile>`), the screen output contains the following information:
  
* a full list of the |topologies| generated by the |decomposition| [*]_ (if **printDecomp** = True). Each |topology| entry
  contains basic information about the |topology| as well as the number of |elements| with this |topology|
  and the sum over all the |elements| weights. If **addElmentInfo** = True the |elements| belonging to each
  |topology| are also explicitly shown, as well as the |element|'s mass, :ref:`final states <final states>`
  and weight:

.. literalinclude:: /images/screenoutput.txt
   :lines: 13-28
   

* a list of all the |analyses| considered (if **printAnalyses** = True). Note that this list correspond to all the analyses
  selected in the *database* options (see :ref:`parameters file <parameterFile>`). If **addAnaInfo** = True,
  for each |analysis| entry a list of all the |elements| appearing in the :ref:`analysis constraint <ULconstraint>`
  is also shown using the :ref:`bracket notation <bracketNotation>`:

.. literalinclude:: /images/screenoutput.txt
   :lines: 6945-6951   

* a list of all the |theory predictions| obtained and the corresponding |analysis| upper limit (if **printResults** = True).
  For each |theory prediction| entry, the :doc:`analysis name <AnalysesNames>`, *sqrts*, the :ref:`cluster <ULcluster>` average mass
  in each branch, the cross-section value and the list of the :ref:`condition values <ULconditions>` for the :ref:`cluster <ULcluster>`
  and the corresponding |analysis| upper limit are shown:

.. literalinclude:: /images/screenoutput.txt
   :lines: 7188-7194
   
* possible (mostly harmless) warnings [*]_

.. literalinclude:: /images/screenoutput.txt
   :lines: 1,3,7   


.. _fileOut:

Summary File Output
*******************

If both **expandedSummary** and **addConstraintInfo** are set to True
(see :ref:`parameter file <parameterFile>`), the file output contains the following information:

* a status flag for the input file and the decomposition indicating possible problems. 
  These flags should be consulted in case of unexpected/missing results:

.. literalinclude:: /images/summary.txt
   :lines: 1-2

* the name of the input file and a few important :ref:`input parameters <parameterFile>`:


.. literalinclude:: /images/summary.txt
   :lines: 3-6

* the version of the :doc:`database <Database>` used to obtain the results:


.. literalinclude:: /images/summary.txt
   :lines: 7


* the list of |analyses| which constrain the input model.
  For each |analysis|, the :doc:`analysis name <AnalysesNames>`, the value of the |analysis|
  center-of-mass energy (*sqrts*), the amount of :ref:`condition <ULconditions>` violation,
  the |theory prediction| value (value for the relevant signal cross-section),
  the experimental upper limit and the ratio (*r*) of the signal cross-section and the
  upper limit (:math:`r = theory\, prediction/upper\, limit`) are printed. 
  A value of :math:`r\ge 1` indicates that the model is likely excluded by the corresponding |analysis|.
  
.. literalinclude:: /images/summary.txt
   :lines: 8-11
   
* if  **addConstraintInfo** = True, the :ref:`analysis constraint <ULconstraint>` in :ref:`bracket notation <bracketNotation>`
  is also included just below the |analysis| entry:
  
.. literalinclude:: /images/summary.txt
   :lines: 12
     
* The last line of this block gives the maximum value of *r*, :math:`R = max(r)`.  

.. literalinclude:: /images/summary.txt
   :lines: 45-46

* if **findMissingTopos** = True,  a list of the :ref:`missing topologies <topCoverage>` and their cross sections at the given
  value of *sqrts* is also included. This list represents the |elements| or sum of |elements| 
  (shown using the :ref:`bracket notation <bracketNotation>`) with the highest
  weights (:math:`\sigma \times BR`) which are not tested by any |analysis|:


.. literalinclude:: /images/summary.txt
   :lines: 48-53

.. _exampleCode:

Example.py
----------

Although :ref:`runSModelS.py <runSModelS>` provides the main SModelS features with a command line interface, 
users more familiar with Python and the SModelS language may prefer to write their own main program. 
A simple example code for this purpose is provided in :download:`examples/Example.py`.
Below we go step-by-step through this example code:

* *Import the SModelS methods*. Import the methods to be used later. If the file is not located in the smodels
  installation folder simply add "sys.path.append(<smodels installation path>)" before importing smodels

.. literalinclude:: /examples/Example.py
   :lines: 11-19

* *Set the address to the dabase*. Specify where the SModelS :doc:`database <Database>` has been installed

.. literalinclude:: /examples/Example.py
   :lines: 21-22
   
* *Path to the input file*. Specify the location of the input file. It must be a SLHA or LHE file (see :ref:`Basic Input <BasicInput>`)

.. literalinclude:: /examples/Example.py
   :lines: 32
   
* *Define the basic parameters for* |decomposition|. Specify the values of :ref:`sigmacut <minweight>` and :ref:`minmassgap <massComp>`:   

.. literalinclude:: /examples/Example.py
   :lines: 36-37
   
* *Perform the* |decomposition|. Depending on the type
  of input format, choose either the `slhaDecomposer.decompose <../../../documentation/build/html/theory.html#theory.slhaDecomposer.decompose>`_ or
  `lheDecomposer.decompose <../../../documentation/build/html/theory.html#theory.slhaDecomposer.decompose>`_ method. The **doCompress** and **doInvisible** options turn on/off the |mass compression| and |invisible compression|, respectively
  
.. literalinclude:: /examples/Example.py
   :lines: 40
   
* *Print the decomposition output*. Set outputLevel = 0 (no output), 1 (basic output) or 2 (extended output)   

.. literalinclude:: /examples/Example.py
   :lines: 44
   
* *Load the the analyses* :ref:`database <Database>`. Load the experimental |analyses| and store the list of analyses   

.. literalinclude:: /examples/Example.py
   :lines: 47
   
* *Compute the* |theory predictions|. For each analysis in list of analyses compute the |theory predictions|. The output
  is a list of `theory prediction objects <../../../documentation/build/html/theory.html#theory.theoryPrediction.TheoryPrediction>`_
  (for each analysis) with results for each :ref:`cluster <ULcluster>`
  
.. literalinclude:: /examples/Example.py
   :lines: 50
   
* *Print the output*. Loop over all analyses and results and print the |theory predictions| information

.. literalinclude:: /examples/Example.py
   :lines: 53-61
   
* *Get analysis upper limit*. For each anlysis and |theory prediction|, obtain the experimental upper limit. This value can
  be compared to the |theory prediction| value to decide if a model is excluded or not.

.. literalinclude:: /examples/Example.py
   :lines: 64

.. [*] For an SLHA :ref:`input file <BasicInput>`, the decay of :ref:`final states <final states>` (or Z\ :sub:`2`-even particles
       such as the Higgs, W,...) are always ignored during the decomposition. Furthermore, if there are two cross-sections
       at different calculation order (say LO and NLO) for the same process, only the highest order is used.
.. [*] The list of |elements| can be extremely long. Try setting **addElmentInfo** = False and/or **printDecomp** = False to obtain
       a smaller output.       
