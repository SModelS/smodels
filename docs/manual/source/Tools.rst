.. index:: SModelS Tools

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
.. |intermediate states| replace:: :ref:`intermediate states <odd states>`
.. |final states| replace:: :ref:`final states <final states>`
.. |database| replace:: :doc:`database <Database>`
.. |bracket notation| replace:: :ref:`bracket notation <bracketNotation>`
.. |ExpRes| replace:: :ref:`Experimental Result<ExpResult>`
.. |ExpRess| replace:: :ref:`Experimental Results<ExpResult>`
.. |Database| replace:: :ref:`Database <Database>`
.. |Dataset| replace:: :ref:`Data Set<DataSet>`
.. |Datasets| replace:: :ref:`Data Sets<DataSet>`

SModelS Tools
=============

Inside SModelS there are a number of tools that may be convenient for the user:

* a :ref:`cross section calculator <xsecCalc>` based on `Pythia6 <http://home.thep.lu.se/~torbjorn/Pythia.html>`_  and 
  `NLLfast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_,
* :ref:`SLHA and LHE file checkers <fileChecks>` to check your input files for completeness and sanity,
* a :ref:`Database Browser <databaseBrowser>` to provide easy access to the |Database| of experimental results,
* a module for identifying :ref:`missing topologies <topCoverage>`.

.. _xsecCalc:

Cross-Section Calculator
------------------------

This little tool computes LHC production cross-sections for *MSSM particles*
and writes them out in :ref:`SLHA convention <xsecblock>`. This can in particular be convenient for adding cross-sections to SLHA input files, see :doc:`Basic Input <BasicInput>`. The calculation is done at LO with `Pythia6.4 <http://home.thep.lu.se/~torbjorn/Pythia.html>`_ ; K-factors for colored particles are computed with `NLLfast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_ .


**The usage of the cross-section calculator is:**

smodelsTools.py xseccomputer [-h] -f FILENAME [-s SQRTS [SQRTS ...]] [-e NEVENTS] [-p] [-k] [-n] [-N] [-O]

*arguments*:
  -h, --help            show this help message and exit
  -s SQRTS, --sqrts SQRTS
                        sqrt(s) TeV. Can supply more than one value.
  -e NEVENTS, --nevents NEVENTS
                        number of events to be simulated.
  -p, --tofile          write cross sections to file
  -k, --keep            do not unlink temporary directory
  -n, --NLO             compute at the NLO level (default is LO)
  -N, --NLL             compute at the NLO+NLL level (takes precedence over NLO,
                        default is LO)
  -O, --LOfromSLHA      use LO cross-sections from file to compute the NLO or
                        NLO+NLL cross-sections
  -f FILENAME, --filename FILENAME
                        SLHA file to compute cross sections for

Some more explanations:

* *-s* (int): an integer with the value (in TeV) of the LHC center-of-mass energy for computing the cross-sections
* *-e* (int): the number of Monte Carlo events when running Pythia
* *-p*: if set, the cross-sections will be written back to the file. If in the input file already
  contains cross-sections, only the non-overlapping ones will be written. If not set, the cross-sections
  will be written to the screen.
* *-k*: if set, keep the temporary directory containing the Pythia run output. This option is only
  relevant when checking for errors when running Pythia.
* *-n*: if set, use Pythia and NLLfast to compute NLO cross-sections. Note that since NLLfast only contains
  results for production of squarks and gluinos, only these cross-sections will be generated
* *-N*: if set, use Pythia and NLLfast to compute NLO+NLL cross-sections. Note that since NLLfast only contains
  results for production of squarks and gluinos, only these cross-sections will be generated
* *-O*: if set, SModelS will read the LO cross-sections from the input file
  and use NLLfast to compute the NLO or NLO+NLL cross-sections for squarks and gluinos
* *-f*: name of input SLHA file

Further Pythia parameters are defined in :download:`etc/pythia.card </images/pythia.card>`.

A typical
usage example is: ::

   smodelsTools.py xseccomputer -s 8 -e 10000 -p -f compressedSpec.slha

which will compute 8 TeV LO cross-sections (at the LHC) for all MSSM processes using 10k MC events.
If, *after* the LO cross-sections have been computed, one wants to add the NLO+NLL cross-sections for gluinos and squarks: ::

   smodelsTools.py xseccomputer -s 8 -p -N -O -f compressedSpec.slha

The resulting file will then contain LO cross-sections for all MSSM processes and NLO+NLL cross-sections for gluinos and squarks.
When reading the input file, SModelS will then use only the highest order cross-sections available for each process.


.. _fileChecks:

Input File Checks
-----------------

As discussed in :doc:`Basic Input <BasicInput>`, SModelS accepts both SLHA and LHE input files. It can be convenient to perform certain sanity checks on these files as described below.

* **The input file checks are implemented by the** `FileStatus class <../../../documentation/build/html/tools.html#tools.ioObjects.FileStatus>`_

.. _lheChecks:

LHE File Checker
^^^^^^^^^^^^^^^^

For a LHE input file only very basic checks are performed, namely that

- the file exists,

- it contains at least one event,

- the information on the total cross section and the center of mass energy can be found.


**The usage of the LHE checker is simply:**

smodelsTools.py lhechecker [-h] -f FILENAME

*arguments*:

  -h, --help            show this help message and exit
  
  -f FILENAME, --filename FILENAME
  

A typical
usage example is: ::

   smodelsTools.py lhechecker -f gluino_squarks.lhe

.. _slhaChecks:

SLHA File Checker
^^^^^^^^^^^^^^^^^

The SLHA file checker allows to perform quite rigorous checks of SLHA input files. Concretely, it verifies that

* the file exists and is given in SLHA format,

* the file contains masses and decay branching ratios in standard SLHA format

* the file contains cross-sections according to the :ref:`SLHA format for cross-sections <xsecSLHA>`,


* the lightest :ref:`Z2-odd state <odd states>` (the LSP in supersymmetric models) is neutral,

* there are no stable charged particles nor displaced vertices (no non-prompt visible decays), as currently all the analyses considered by SModelS require a prompt MET signature.

In addition, one can ask that

* all decays listed in the DECAY block are kinematically allowed, *i.e.* the sum of masses of the decay products may not exceed the mother mass. *NB This check for "illegal decays" is turned off by default.*

If any of the above tests fail (return a negative result), an error message is shown.

Some more comments are in order.
In order to check that the lightest Z\ :sub:`2`-odd state has zero electric and color charges, the quantum numbers of the BSM particles must be given in the
``qNumbers`` dictionary in :download:`particles.py <images/particles.py>`. The format is

``[2*spin, 3*electric charge, dimension of SU(3) representation]``

The list of quantum numbers is also required to check for displaced vertices or heavy charged particles.
The check for long-lived (or stable) particles first verifies if these
appear in one of the cross-section blocks and their cross-section
exceeds the minimum cross-section value defined by :ref:`sigmacut <parameterFile>` (see  :ref:`Minimum Decomposition Weight <minweight>`).
If the cross-section is larger than sigmacut and the particle is stable,
the checker verifies if it is neutral (both electric and color charges
are zero). On the other hand, if the particle is unstable, but its lifetime (times *c*)
is larger than a minimum value (*default = 10 mm*), the particle is considered
as a non-prompt decay.
For non-prompt decays, all channels are then checked for visible decay products.
If the branching ratio to visible decays times the maximum production cross-section
for the particle exceeds :ref:`sigmacut <parameterFile>`, the particle's decay
is considered as a displaced vertex.


**The usage of the SLHA checker is:**

smodelsTools.py slhachecker [-h] [-xS] [-lsp] [-longlived] [-m DISPLACEMENT] [-s SIGMACUT] [-illegal] -f FILENAME

*arguments*:
  -h, --help            show this help message and exit
  -xS, --xsec           turn off the check for xsection blocks
  -lsp, --lsp           turn off the check for charged lsp
  -longlived, --longlived
                        turn off the check for stable charged particles and
                        visible displaced vertices
  -m DISPLACEMENT, --displacement DISPLACEMENT
                        give maximum displacement of secondary vertex in m
  -s SIGMACUT, --sigmacut SIGMACUT
                        give sigmacut in fb
  -illegal, --illegal   turn on check for kinematically forbidden decays
  -dB, --decayBlocks    turn off the check for missing decay blocks
  -f FILENAME, --filename FILENAME
                        name of input SLHA file


In some more detail:

* *-f*: path to the input file
* *-xS*: if this flag is set, the check for a cross section block will not be performed
* *-lsp*: if this flag is set, the check for a neutral LSP will not be performed
* *-longlived*: if this flag is set, check for non-prompt visible decays or stable charged particles will not be performed
* *-m* (float): use this to set the value of c*tau (in meters) where a decay is no longer considered prompt
* *-s* (float): use this to set the value of sigmacut, that is used as a cutoff for relevant non-promt decays or long lived charged particle production
* *-illegal*: if this flag is set, the check for illegal (kinematically forbidden) decays will be performed
* *-dB*: if this flag is set, the check for missing decay blocks will not be performed

A typical
usage example is: ::

   smodelsTools.py slhachecker -m 0.001 -s 0.01 -f lightSquarks.slha

Running this will print the status flag and a message with potential warnings
and error messages.

.. _databaseBrowser:

Database Browser
----------------

In several cases the user might be interested in an easy way to directly access the |Database| of |ExpRess|.
This can be conveniently done using the database browser. The browser owns several methods to select  |ExpRess|
or |Datasets| satisfying some user-defined conditions as well as to access the meta data and data inside each
|ExpRes|. 

Unlike most of the other SModelS tools, the browser can not be directly accessed from the command line.
However, several python example codes using the functionalities of the browser
can be found in :ref:`More Examples <Examples>`.
Below we will quickly describe how to instantiate the browser and use its main functionalities.

* **The Database browser tool is implemented by the**  `Browser class <../../../documentation/build/html/tools.html#tools.databaseBrowser.Browser>`_


.. _topCoverage:

Topology Coverage
-----------------

Unlike the :ref:`file checks <fileChecks>` and the :ref:`cross-section calculator <xsecCalc>`, the missing topologies tool can be called only *after* the SMS |decomposition| and |theory predictions| have been computed.
Given the |decomposition| output (list of |elements|), as well as the |database|
information, it finds the |elements| which are
not tested by any of the |analyses| in the |database|.

To this end, the tool loops over all the |elements| found in the
|decomposition| and checks if they are tested by one or more |analyses| in the |database|.
If :ref:`mass <massComp>` or :ref:`invisible compression <invComp>`
are turned on, elements which can be :ref:`compressed <elementComp>` are not considered, to avoid double counting.
All the |elements| not appearing in any of the |constraints| in the |database| are then marked
as "missing". A missing topology is then characterized
by a sum over the missing |elements| differing only by their
masses (with the same |final states|) or electric charges.

The missing topologies tool is normally called from within SModelS (e.g. when running :ref:`runSModelS.py <runSModelS>`) by setting **findMissingTopos=True**
in the :ref:`parameters file <parameterFile>` .
In the output, the missing topologies are ordered by cross section. By default only the ones with the ten largest cross-sections are shown.

* **The missing topologies tool is implemented by the** `MissingTopoList class <../../../documentation/build/html/tools.html#tools.missingTopologies.MissingTopoList>`_ 

