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
.. |database| replace:: :doc:`database <Database>`
.. |bracket notation| replace:: :ref:`bracket notation <bracketNotation>`
.. |ExpRes| replace:: :ref:`Experimental Result<ExpResult>`
.. |ExpRess| replace:: :ref:`Experimental Results<ExpResult>`
.. |Database| replace:: :ref:`Database <Database>`
.. |database| replace:: :ref:`database <Database>`
.. |Dataset| replace:: :ref:`DataSet<DataSet>`
.. |Datasets| replace:: :ref:`DataSets<DataSet>`
.. |results| replace:: :ref:`experimental results <ExpResult>`
.. |branches| replace:: :ref:`branches <branch>`
.. |branch| replace:: :ref:`branch <branch>`
.. |EMrs| replace:: :ref:`EM-type results <EMtype>`
.. |ULrs| replace:: :ref:`UL-type results <ULtype>`
.. |br| raw:: html
   <br />

SModelS Tools
=============

Inside SModelS there is a number of tools that may be convenient for the user:

* a :ref:`cross section calculator <xsecCalc>` based on `Pythia8 <http://home.thep.lu.se/~torbjorn/Pythia.html>`_ (or `Pythia6 <http://pythia6.hepforge.org>`_) and 
  `NLLfast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_,
* :ref:`SLHA and LHE file checkers <fileChecks>` to check your input files for completeness and sanity,
* a :ref:`database Browser <databaseBrowser>` to provide easy access to the |database| of experimental results,
* a module for identifying :ref:`missing topologies <topCoverage>`.

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



*In some more detail*:
  -s SQRTS, --sqrts SQRTS 
                        (int) an integer (or integers) with the value (in TeV) of the LHC center-of-mass energy for computing the cross sections
  -e NEVENTS, --nevents NEVENTS 
                        (int) the number of Monte Carlo events when running Pythia
  -c NCPUS, --ncpus NCPUS 
                        (int) number of cpu cores to be used.
                        It is only used when cross sections are computed for multiple SLHA files.
  -p, --tofile          if set, the cross sections will be written back to the file. 
    If in the input file already contains cross sections, only the non-overlapping ones will be written. 
    If not set, the cross sections will be written to the screen.
  -q, --query           if set, will only check if the input file already contains cross sections.  
  -k, --keep            if set, keep the temporary directory containing the Pythia run output. This option is only relevant when checking for errors when running Pythia.
  -n, --NLO             if set, use Pythia and NLLfast to compute NLO cross sections. Note that since NLLfast only contains results for production of squarks and gluinos, only these cross sections will be generated
  -N, --NLL             if set, use Pythia and NLLfast to compute NLO+NLL cross sections. 
                        Note that since NLLfast only contains results for production of squarks and gluinos, only these cross sections will be generated
  -O, --LOfromSLHA      if set, SModelS will read the LO cross sections from the input file and use NLLfast to compute the NLO or NLO+NLL cross sections for squarks and gluinos
  -f FILENAME, --filename FILENAME
                        name of input SLHA file or a folder containing SLHA files

Further Pythia parameters are defined in :download:`etc/pythia.card </images/pythia.card>`.

A typical
usage example is: ::

   smodelsTools.py xseccomputer -s 8 13 -e 10000 -p -f inputFiles/slha/compressedSpec.slha

which will compute 8 TeV and 13 TeV LO cross sections (at the LHC) for all MSSM processes using 10k MC events.
If, *after* the LO cross sections have been computed, one wants to add the NLO+NLL cross sections for gluinos and squarks: ::

   smodelsTools.py xseccomputer -s 8 13 -p -N -O -f inputFiles/slha/compressedSpec.slha

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
   ['ATLAS-SUSY-2015-09', 'CMS-SUS-PAS-15-002', 'ATLAS-CONF-2012-105', 'ATLAS-CONF-2012-166', 'ATLAS-CONF-2013-001', 'ATLAS-CONF-2013-007', 'ATLAS-CONF-2013-024', 'ATLAS-CONF-2013-024', 'ATLAS-CONF-2013-025', 'ATLAS-CONF-2013-028', 'ATLAS-CONF-2013-035', 'ATLAS-CONF-2013-035', 'ATLAS-CONF-2013-036', 'ATLAS-CONF-2013-037', 'ATLAS-CONF-2013-037', 'ATLAS-CONF-2013-047', 'ATLAS-CONF-2013-047', 'ATLAS-CONF-2013-048', 'ATLAS-CONF-2013-048', 'ATLAS-CONF-2013-049', 'ATLAS-CONF-2013-049', 'ATLAS-CONF-2013-053', 'ATLAS-CONF-2013-053', 'ATLAS-CONF-2013-054', 'ATLAS-CONF-2013-061', 'ATLAS-CONF-2013-061', 'ATLAS-CONF-2013-062', 'ATLAS-CONF-2013-062', 'ATLAS-CONF-2013-065', 'ATLAS-CONF-2013-089', 'ATLAS-CONF-2013-093', 'ATLAS-CONF-2013-093', 'ATLAS-SUSY-2013-02', 'ATLAS-SUSY-2013-02', 'ATLAS-SUSY-2013-04', 'ATLAS-SUSY-2013-04', 'ATLAS-SUSY-2013-05', 'ATLAS-SUSY-2013-05', 'ATLAS-SUSY-2013-08', 'ATLAS-SUSY-2013-09', 'ATLAS-SUSY-2013-09', 'ATLAS-SUSY-2013-11', 'ATLAS-SUSY-2013-11', 'ATLAS-SUSY-2013-12', 'ATLAS-SUSY-2013-14', 'ATLAS-SUSY-2013-15', 'ATLAS-SUSY-2013-15', 'ATLAS-SUSY-2013-16', 'ATLAS-SUSY-2013-16', 'ATLAS-SUSY-2013-18', 'ATLAS-SUSY-2013-18', 'ATLAS-SUSY-2013-19', 'ATLAS-SUSY-2013-21', 'ATLAS-SUSY-2013-23', 'ATLAS-SUSY-2014-03', 'CMS-PAS-SUS-12-022', 'CMS-PAS-SUS-12-026', 'CMS-PAS-SUS-13-015', 'CMS-PAS-SUS-13-015', 'CMS-PAS-SUS-13-016', 'CMS-PAS-SUS-13-016', 'CMS-PAS-SUS-13-018', 'CMS-PAS-SUS-13-023', 'CMS-PAS-SUS-14-011', 'CMS-SUS-12-024', 'CMS-SUS-12-024', 'CMS-SUS-12-028', 'CMS-SUS-13-002', 'CMS-SUS-13-004', 'CMS-SUS-13-006', 'CMS-SUS-13-006', 'CMS-SUS-13-007', 'CMS-SUS-13-007', 'CMS-SUS-13-011', 'CMS-SUS-13-011', 'CMS-SUS-13-012', 'CMS-SUS-13-013', 'CMS-SUS-13-013', 'CMS-SUS-13-019', 'CMS-SUS-14-010', 'CMS-SUS-14-021', 'CMS-SUS-14-021']
   
   In [2]: browser.selectExpResultsWith(txName = 'T1tttt', dataType = 'upperLimit') #Select only the UL results with the topology T1tttt
   
   In [3]: print browser #Print all experimental results in the browser (after selection)
   ['ATLAS-SUSY-2015-09', 'CMS-SUS-PAS-15-002', 'ATLAS-CONF-2012-105', 'ATLAS-CONF-2013-007', 'ATLAS-CONF-2013-061', 'ATLAS-SUSY-2013-04', 'ATLAS-SUSY-2013-09', 'ATLAS-SUSY-2013-18', 'CMS-PAS-SUS-12-026', 'CMS-PAS-SUS-13-016', 'CMS-PAS-SUS-14-011', 'CMS-SUS-12-024', 'CMS-SUS-12-028', 'CMS-SUS-13-002', 'CMS-SUS-13-004', 'CMS-SUS-13-007', 'CMS-SUS-13-012', 'CMS-SUS-13-013', 'CMS-SUS-13-019', 'CMS-SUS-14-010']
   
   In [4]: gluinoMass, LSPmass = 800.*GeV, 100.*GeV  #Define masses for the T1tttt topology
   
   In [5]: browser.getULFor('CMS-SUS-PAS-15-002','T1tttt',[[gluinoMass,LSPmass],[gluinoMass,LSPmass]]) #Get UL for a specific experimental result
   Out[5]: 5.03E-02 [pb]
     
   In [6]: for expResult in browser:  #Get the upper limits for all the selected results for the given topology and mass
      ...:     print expResult.getValuesFor('id'),'UL = ',expResult.getUpperLimitFor(txname='T1tttt',mass=[[gluinoMass,LSPmass],[gluinoMass,LSPmass]])
      ...:     
   ['ATLAS-SUSY-2015-09'] UL =  None
   ['CMS-SUS-PAS-15-002'] UL =  5.03E-02 [pb]
   ['ATLAS-CONF-2012-105'] UL =  6.70E-02 [pb]
   ['ATLAS-CONF-2013-007'] UL =  2.40E-02 [pb]
   ['ATLAS-CONF-2013-061'] UL =  1.25E-02 [pb]
   ['ATLAS-SUSY-2013-04'] UL =  1.40E-02 [pb]
   ['ATLAS-SUSY-2013-09'] UL =  1.73E-02 [pb]
   ['ATLAS-SUSY-2013-18'] UL =  4.30E-03 [pb]
   ['CMS-PAS-SUS-12-026'] UL =  4.60E-02 [pb]
   ['CMS-PAS-SUS-13-016'] UL =  3.55E-02 [pb]
   ['CMS-PAS-SUS-14-011'] UL =  2.47E-02 [pb]
   ['CMS-SUS-12-024'] UL =  3.62E-02 [pb]
   ['CMS-SUS-12-028'] UL =  5.31E-02 [pb]
   ['CMS-SUS-13-002'] UL =  3.48E-02 [pb]
   ['CMS-SUS-13-004'] UL =  2.47E-02 [pb]
   ['CMS-SUS-13-007'] UL =  6.00E-03 [pb]
   ['CMS-SUS-13-012'] UL =  2.14E-02 [pb]
   ['CMS-SUS-13-013'] UL =  1.90E-02 [pb]
   ['CMS-SUS-13-019'] UL =  1.35E-02 [pb]
   ['CMS-SUS-14-010'] UL =  4.82E-03 [pb]
      
   In [7]: for expResult in browser:  #Print the luminosities for the selected experimental results
      ...:     print expResult.getValuesFor('id'),expResult.getValuesFor('lumi')
      ...:     
   ['ATLAS-SUSY-2015-09'] [3.20E+00 [1/fb]]
   ['CMS-SUS-PAS-15-002'] [2.20E+00 [1/fb]]
   ['ATLAS-CONF-2012-105'] [5.80E+00 [1/fb]]
   ['ATLAS-CONF-2013-007'] [2.07E+01 [1/fb]]
   ['ATLAS-CONF-2013-061'] [2.01E+01 [1/fb]]
   ['ATLAS-SUSY-2013-04'] [2.03E+01 [1/fb]]
   ['ATLAS-SUSY-2013-09'] [2.03E+01 [1/fb]]
   ['ATLAS-SUSY-2013-18'] [2.01E+01 [1/fb]]
   ['CMS-PAS-SUS-12-026'] [9.20E+00 [1/fb]]
   ['CMS-PAS-SUS-13-016'] [1.97E+01 [1/fb]]
   ['CMS-PAS-SUS-14-011'] [1.93E+01 [1/fb]]
   ['CMS-SUS-12-024'] [1.94E+01 [1/fb]]
   ['CMS-SUS-12-028'] [1.17E+01 [1/fb]]
   ['CMS-SUS-13-002'] [1.95E+01 [1/fb]]
   ['CMS-SUS-13-004'] [1.93E+01 [1/fb]]
   ['CMS-SUS-13-007'] [1.93E+01 [1/fb]]
   ['CMS-SUS-13-012'] [1.95E+01 [1/fb]]
   ['CMS-SUS-13-013'] [1.95E+01 [1/fb]]
   ['CMS-SUS-13-019'] [1.95E+01 [1/fb]]
   ['CMS-SUS-14-010'] [1.95E+01 [1/fb]]



Further Python example codes using the functionalities of the browser
can be found in :ref:`Howto's <Examples>`.

* **The Database browser tool is implemented by the**  `Browser class <../../../documentation/build/html/tools.html#tools.databaseBrowser.Browser>`_


.. _topCoverage:

Topology Coverage
-----------------

Unlike the :ref:`database browser <databaseBrowser>`, the :ref:`file checks <fileChecks>` and the :ref:`cross section calculator <xsecCalc>`, 
the topology coverage tool can not be independently accessed.
It requires the output from the SMS |decomposition| and |theory predictions|.
Given the |decomposition| output (list of |elements|), as well as the |database|
information, it finds and classifies the |elements| which are
not tested by any of the |results| in the |database|.
These elements are grouped into the following classes:

* *missingTopos*: |elements| which are not tested by any of the |results| in the |database| (independent of the element mass).
  The missing topologies are further classified as:
   * *longCascade*: |elements| with long cascade decays (more than one intermediate particle in one of the |branches|);
   * *asymmetricBranches*: |elements| where the first |branch| differs from the second |branch| (but that are not considered as long cascade decays).

* *outsideGrid*: |elements| which could be tested by one or more experimental result, but are not constrained because the mass array is outside the mass grid;

In order to classify the |elements|, the tool loops over all the |elements| found in the
|decomposition| and checks if they are tested by one or more |results| in the |database| [*]_.
All the |elements| which are not tested by any of the |results| in the |database| (independent of their masses)
are added to the *missingTopos* class.
The remaining |elements| which do appear in one or more of the |results|, but have
not been tested because their masses fall outside the efficiency or upper limit grids (see |EMrs| and |ULrs|),
are added to the *outsideGrid* class.


Usually the list of  *missing* or *outsideGrid* elements is considerably long.
Hence, to compress this list, all |elements| differing only by their
masses (with the same |final states|) or electric charges are combined. Moreover, by default, electrons and muons
are combined to light leptons (denoted "l"): gluons and light quarks are combined into jets.
The *missing* topologies are then further classified (if applicable) into *longCascade* or *asymmetricBranches* topologies.


The topologies for each of the four categories are then grouped according to the final state (for the *missingTopos* and
*outsideGrid* classes) or according to the PDG ids of the initially produced motherparticles (for the *longCascade* and
*asymmetricBranches* classes). 
We note that for the latter the |elements| deriving from different mother particles, but with the same |final states| and mass configuration cannot be distinguished, and are therefore combined in this grouping.
The full list of mother PDG id pairs can be accessed in the python printout or the comment of the text printout.


The topology coverage tool is normally called from within SModelS (e.g. when running :ref:`runSModelS.py <runSModelS>`) by setting **testCoverage=True**
in the :ref:`parameters file <parameterFile>`.
In the output, contributions in each category are ordered by cross section. 
By default only the ones with the ten largest cross sections are shown.

* **The topology coverage tool is implemented by the** `Uncovered class <../../../documentation/build/html/tools.html#tools.coverage.Uncovered>`_ 


.. [*] If :ref:`mass <massComp>` or :ref:`invisible compression <invComp>` are turned on, elements which can be :ref:`compressed <elementComp>` are not considered, to avoid double counting.
