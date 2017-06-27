.. index:: SModelS Tools

.. |element| replace:: :ref:`element <element>`
.. |elements| replace:: :ref:`elements <element>`
.. |topology| replace:: :ref:`topology <topology>`
.. |topologies| replace:: :ref:`topologies <topology>`
.. |decomposition| replace:: :doc:`decomposition <Decomposition>`
.. |theory predictions| replace:: :doc:`theory predictions <TheoryPredictions>`
.. |theory prediction| replace:: :doc:`theory prediction <TheoryPredictions>`
.. |constraint| replace:: :ref:`constraint <ULconstraint>`
.. |constraints| replace:: :ref:`constraints <ULconstraint>`
.. |intermediate states| replace:: :ref:`intermediate states <odd states>`
.. |final states| replace:: :ref:`final states <final states>`
.. |database| replace:: :ref:`database <Database>`
.. |bracket notation| replace:: :ref:`bracket notation <bracketNotation>`
.. |ExpRes| replace:: :ref:`Experimental Result<ExpResult>`
.. |ExpRess| replace:: :ref:`Experimental Results<ExpResult>`
.. |Database| replace:: :ref:`Database <Database>`
.. |Dataset| replace:: :ref:`DataSet<DataSet>`
.. |Datasets| replace:: :ref:`DataSets<DataSet>`
.. |results| replace:: :ref:`experimental results <ExpResult>`
.. |branches| replace:: :ref:`branches <branch>`
.. |branch| replace:: :ref:`branch <branch>`
.. |EMrs| replace:: :ref:`EM-type results <EMtype>`
.. |ULrs| replace:: :ref:`UL-type results <ULtype>`

.. _smodelsTools:

SModelS Tools
=============

Inside SModelS there is a number of tools that may be convenient for the user:

* a :ref:`cross section calculator <xsecCalc>` based on `Pythia8 <http://home.thep.lu.se/~torbjorn/Pythia.html>`_ (or `Pythia6 <http://pythia6.hepforge.org>`_) and 
  `NLLfast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_,
* :ref:`SLHA and LHE file checkers <fileChecks>` to check your input files for completeness and sanity,
* a :ref:`database Browser <databaseBrowser>` to provide easy access to the |database| of experimental results.

.. _xsecCalc:

Cross Section Calculator
------------------------

This tool computes LHC production cross sections for *MSSM particles*
and writes them out in :ref:`SLHA convention <xsecblock>`. This can in particular be 
convenient for adding cross sections to SLHA input files, see :doc:`Basic Input <BasicInput>`. 
The calculation is done at LO with `Pythia8 <http://home.thep.lu.se/~torbjorn/Pythia.html>`_ or `Pythia6.4 <http://pythia6.hepforge.org>`_ ; K-factors 
for colored particles are computed with `NLLfast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_ .


**The usage of the cross section calculator is:**

.. include:: XSecComputer.rst


Further Pythia parameters are defined in :download:`smodels/etc/pythia8.cfg </images/pythia8.cfg>` (for Pythia 8)
or :download:`smodels/etc/pythia.card </images/pythia.card>` (for Pythia 6).
.

A typical
usage example is: ::

   smodelsTools.py xseccomputer -s 8 13 -e 10000 -p -f inputFiles/slha/higgsinoStop.slha

which will compute 8 TeV and 13 TeV LO cross sections (at the LHC) for all MSSM processes using 10k MC events.
If, *after* the LO cross sections have been computed, one wants to add the NLO+NLL cross sections for gluinos and squarks: ::

   smodelsTools.py xseccomputer -s 8 13 -p -N -O -f inputFiles/slha/higgsinoStop.slha

The resulting file will then contain LO cross sections for all MSSM processes and NLO+NLL cross sections for 
the available processes in `NLLfast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_  
(gluino and squark production).
When reading the input file, SModelS will then use only the highest order cross sections available for each process.

* **The cross section calculation is implemented by the** `computeXSec function <../../../documentation/build/html/tools.html#tools.xsecComputer.computeXSec>`_


.. _fileChecks:

Input File Checks
-----------------

As discussed in :doc:`Basic Input <BasicInput>`,
SModelS accepts both SLHA and LHE input files. It can be convenient to perform certain sanity checks on these files as described below.

* **The input file checks are implemented by the** `FileStatus class <../../../documentation/build/html/tools.html#tools.ioObjects.FileStatus>`_

.. _lheChecks:

LHE File Checker
^^^^^^^^^^^^^^^^

For a LHE input file only very basic checks are performed, namely that

- the file exists,

- it contains at least one event,

- the information on the total cross section and the center of mass energy can be found.


**The usage of the LHE checker is simply:**

.. include:: LheChecker.rst

A typical
usage example is: ::

   smodelsTools.py lhechecker -f inputFiles/slha/gluino_squarks.lhe

.. _slhaChecks:

SLHA File Checker
^^^^^^^^^^^^^^^^^

The SLHA file checker allows to perform quite rigorous checks of SLHA input files. Concretely, it verifies that

* the file exists and is given in SLHA format,

* the file contains masses and decay branching ratios in standard SLHA format,

* the file contains cross sections according to the :ref:`SLHA format for cross sections <xsecSLHA>`,


* the lightest :ref:`Z2-odd state <odd states>` (the LSP in supersymmetric models) is neutral,

* there are no stable charged particles nor displaced vertices (no non-prompt visible decays), as currently all the analyses considered by SModelS require a prompt MET signature.

In addition, one can ask that

* all decays listed in the DECAY block are kinematically allowed, *i.e.* the sum of masses of the decay products may not exceed the mother mass. *This check for "illegal decays" is turned off by default.*

If any of the above tests fail (return a negative result), an error message is shown.

Some more comments are in order.
In order to check that the lightest Z\ :sub:`2`-odd state has zero electric and color charges, the quantum numbers of the BSM particles must be given in the
``qNumbers`` dictionary in :download:`particles.py <images/particles.py>`. The format is

``[2*spin, 3*electric charge, dimension of SU(3) representation]``

The list of quantum numbers is also required to check for displaced vertices or heavy charged particles.
The check for long-lived (or stable) particles first verifies if these
appear in one of the cross section blocks and their cross section
exceeds the minimum cross section value defined by :ref:`sigmacut <parameterFile>` (see  :ref:`Minimum Decomposition Weight <minweight>`).
If the cross section is larger than sigmacut and the particle is stable,
the checker verifies if it is neutral (both electric and color charges
are zero). On the other hand, if the particle is unstable, but its lifetime (times *c*)
is larger than a minimum value (*default = 10 mm*), the particle is considered
as a non-prompt decay.
For non-prompt decays, all channels are then checked for visible decay products.
If the branching ratio to visible decays times the maximum production cross section
for the particle exceeds :ref:`sigmacut <parameterFile>`, the particle's decay
is considered as a displaced vertex.


**The usage of the SLHA checker is:**

.. include:: SlhaChecker.rst


A typical
usage example is: ::

   smodelsTools.py slhachecker -m 0.001 -s 0.01 -f inputFiles/slha/lightSquarks.slha

Running this will print the status flag and a message with potential warnings
and error messages.

.. _databaseBrowser:

Database Browser
----------------

In several cases the user might be interested in an easy way to directly access the |database| of |ExpRess|.
This can be conveniently done using the database browser. The browser owns several methods to select  |ExpRess|
or |Datasets| satisfying some user-defined conditions as well as to access the meta data and data inside each
|ExpRes|.

**The usage of the browser interface is:**


.. include:: DatabaseBrowser.rst


A typical usage example is: ::

    smodelsTools.py database-browser -p ./smodels-database

Loading the database may take a few seconds if the :ref:`binary database file <databasePickle>` exists.
Otherwise the :ref:`pickle file <databasePickle>` will be created.
Starting the browser opens an IPython session, which can be used 
to select specific experimental results (or groups of experimental results),
check upper limits and/or efficiencies for specific masses/topologies and access all the available
information in the database.
A simple example is given below:

.. code-block:: IPython

   In [1]: print browser  #Print all experimental results in the browser
   ['ATLAS-SUSY-2015-01', 'ATLAS-SUSY-2015-01', 'ATLAS-SUSY-2015-02', 'ATLAS-SUSY-2015-02', ...
   
   In [2]: browser.selectExpResultsWith(txName = 'T1tttt', dataType = 'upperLimit') #Select only the UL results with the topology T1tttt
   
   In [3]: print browser #Print all experimental results in the browser (after selection)
   ['ATLAS-SUSY-2015-09', 'CMS-PAS-SUS-15-002', 'CMS-PAS-SUS-16-014', 'CMS-PAS-SUS-16-015', ...
   
   In [4]: gluinoMass, LSPmass = 800.*GeV, 100.*GeV  #Define masses for the T1tttt topology
   
   In [5]: browser.getULFor('CMS-PAS-SUS-15-002','T1tttt',[[gluinoMass,LSPmass],[gluinoMass,LSPmass]]) #Get UL for a specific experimental result
   Out[5]: 5.03E-02 [pb]
     
   In [6]: for expResult in browser[:5]:  #Get the upper limits for the first five of the selected results for the given topology and mass
      ...:     print expResult.getValuesFor('id'),'UL = ',expResult.getUpperLimitFor(txname='T1tttt',mass=[[gluinoMass,LSPmass],[gluinoMass,LSPmass]])
      ...:     
      ['ATLAS-SUSY-2015-09'] UL =  None
      ['CMS-PAS-SUS-15-002'] UL =  5.03E-02 [pb]
      ['CMS-PAS-SUS-16-014'] UL =  4.10E-02 [pb]
      ['CMS-PAS-SUS-16-015'] UL =  1.80E-02 [pb]
      ['CMS-PAS-SUS-16-016'] UL =  5.76E-02 [pb]

      
   In [7]: for expResult in browser[:5]:  #Print the luminosities for the first five of the selected experimental results
      ...:     print expResult.getValuesFor('id'),expResult.getValuesFor('lumi')
      ...:     
      ['ATLAS-SUSY-2015-09'] [3.20E+00 [1/fb]]
      ['CMS-PAS-SUS-15-002'] [2.20E+00 [1/fb]]
      ['CMS-PAS-SUS-16-014'] [1.29E+01 [1/fb]]
      ['CMS-PAS-SUS-16-015'] [1.29E+01 [1/fb]]
      ['CMS-PAS-SUS-16-016'] [1.29E+01 [1/fb]]


Further Python example codes using the functionalities of the browser
can be found in :ref:`Howto's <Examples>`.

* **The Database browser tool is implemented by the**  `Browser class <../../../documentation/build/html/tools.html#tools.databaseBrowser.Browser>`_


