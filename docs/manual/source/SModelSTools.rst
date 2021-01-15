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
.. |final states| replace:: :ref:`final states <final statesEven>`
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
* a :ref:`database browser <databaseBrowser>` to provide easy access to the |database| of experimental results, 
* a plotting tool to make :ref:`interactive plots <interactivePlots>` based on `plotly <https://plot.ly/python/>`_ (v1.1.3 onwards),
* a :ref:`file permissions fixer <permissionsFixer>` to fix a problem with file permissions for the cross section computers in system-wide installs, and 
* a :ref:`toolbox <toolBox>` to quickly show the state of the external tools.

.. _xsecCalc:

Cross Section Calculator
------------------------

This tool computes LHC production cross sections for *MSSM particles*
and writes them out in :ref:`SLHA convention <xsecblock>`. This can in particular be 
convenient for adding cross sections to SLHA input files, see :doc:`Basic Input <BasicInput>`. 
The calculation is done at LO with `Pythia8 <http://home.thep.lu.se/~torbjorn/Pythia.html>`_ or `Pythia6.4 <http://pythia6.hepforge.org>`_ ; K-factors 
for colored particles are computed with `NLLfast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_. Signal strength multipliers can optionally be supplied for each "mother" particle. 

**The usage of the cross section calculator is:**

.. include:: XSecComputer.rst


Further Pythia parameters are defined in :download:`smodels/etc/pythia8.cfg </images/pythia8.cfg>` (for Pythia 8)
or :download:`smodels/etc/pythia.card </images/pythia.card>` (for Pythia 6).

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


An example using signal strength multipliers (*available from SModelS v2.0 onwards*) is: ::

   smodelsTools.py xseccomputer -s 8 13 -e 10000 --ssmultipliers "{ (1000021,1000021): 4.0, (1000001,-10000001): 2.0 }" -p -f inputFiles/slha/higgsinoStop.slha

This will compute 8 TeV and 13 TeV LO cross sections as above, but the cross section for gluino-pair production (pid 1000021) gets enhanced by a factor of 4, and squark-antisquark production gets enhanced by a factor of 2. For the pids, strings can be supplied instead of integers. Unix filename wildcard syntax is also supported. E.g. '100000?' matches all left-handed squarks but no anti-squarks, '\*1000001' matches both (left-handed) down and anti-down. Multiple signal strength multipliers may be applicable to a single theory prediction. 

Note that signal strength multipliers get applied only to LO cross sections. This means they are propagated to NLO and NLL level only iff the LO cross sections are computed first and the NLO/NLL corrections added afterwards. In other words, if the xseccomputer is called with -n or -N argument but without -O (--LOfromSLHA), the --ssmultipliers argument will be ignored. 


* **The cross section calculation is implemented by the** `xsecComputer function <tools.html#tools.xsecComputer.XSecComputer>`_


.. _fileChecks:

Input File Checks
-----------------

As discussed in :doc:`Basic Input <BasicInput>`,
SModelS accepts both SLHA and LHE input files. It can be convenient to perform certain sanity checks on these files as described below.

* **The input file checks are implemented by the** `FileStatus class <tools.html#tools.ioObjects.FileStatus>`_

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

   smodelsTools.py lhechecker -f inputFiles/lhe/gluino_squarks.lhe

.. _slhaChecks:

SLHA File Checker
^^^^^^^^^^^^^^^^^

The SLHA file checker allows to perform quite rigorous checks of SLHA input files. Concretely, it verifies that

* the file exists and is given in SLHA format,

* the file contains masses and decay branching ratios in standard SLHA format,

* the file contains cross sections according to the :ref:`SLHA format for cross sections <xsecSLHA>`,

In addition, one can ask that

* all decays listed in the DECAY block are kinematically allowed, *i.e.* the sum of masses of the decay products may not exceed the mother mass. *This check for "illegal decays" is turned off by default.*

If any of the above tests fail (return a negative result), an error message is shown.


**The usage of the SLHA checker is:**

.. include:: SlhaChecker.rst


A typical
usage example is: ::

   smodelsTools.py slhachecker -f inputFiles/slha/gluino_squarks.slha

Running this will print the status flag and a message with potential warnings
and error messages.

.. note:: In SModelS versions prior to 1.2, the SLHA file checker also
          checked for the existence of displaced vertices or heavy stable charged
          particles in the input file. Since the inclusion of long lived signatures in
          SModelS, these checks are no longer done by the SLHA file checker.
          
          

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

   In [1]: print ( browser )  #Print all experimental results in the browser
   ['ATLAS-SUSY-2015-01', 'ATLAS-SUSY-2015-01', 'ATLAS-SUSY-2015-02', 'ATLAS-SUSY-2015-02', ...
   
   In [2]: browser.selectExpResultsWith(txName = 'T1tttt', dataType = 'upperLimit') #Select only the UL results with the topology T1tttt
   
   In [3]: print ( browser ) #Print all experimental results in the browser (after selection)
   ['ATLAS-SUSY-2015-09', 'CMS-SUS-15-002', 'CMS-PAS-SUS-16-014', 'CMS-PAS-SUS-16-015', ...
   
   In [4]: gluinoMass, LSPmass = 800.*GeV, 100.*GeV  #Define masses for the T1tttt topology
   
   In [5]: browser.getULFor('CMS-SUS-15-002','T1tttt',[[gluinoMass,LSPmass],[gluinoMass,LSPmass]]) #Get UL for a specific experimental result
   Out[5]: 2.23E+01 [fb]
     
   In [6]: for expResult in browser[:5]:  #Get the upper limits for the first five of the selected results for the given topology and mass
      ...:     print ( expResult.getValuesFor('id'),'UL = ',expResult.getUpperLimitFor(txname='T1tttt',mass=[[gluinoMass,LSPmass],[gluinoMass,LSPmass]]) )
      ...:     
      ['ATLAS-SUSY-2015-09'] UL =  None
      ['CMS-PAS-SUS-16-014'] UL =  4.10E+01 [fb]
      ['CMS-PAS-SUS-16-015'] UL =  1.80E+01 [fb]
      ['CMS-PAS-SUS-16-016'] UL =  5.76E+01 [fb]
      ['CMS-PAS-SUS-16-019'] UL =  1.37E+01 [fb]

      
   In [7]: for expResult in browser[:5]:  #Print the luminosities for the first five selected experimental results
      ...:     print ( expResult.getValuesFor('id'),expResult.getValuesFor('lumi') )
      ...:     
      ['ATLAS-SUSY-2015-09'] [3.20E+00 [1/fb]]
      ['CMS-PAS-SUS-16-014'] [1.29E+01 [1/fb]]
      ['CMS-PAS-SUS-16-015'] [1.29E+01 [1/fb]]
      ['CMS-PAS-SUS-16-016'] [1.29E+01 [1/fb]]
      ['CMS-PAS-SUS-16-019'] [1.29E+01 [1/fb]]


Further Python example codes using the functionalities of the browser
can be found in :ref:`Howto's <Examples>`.

* **The Database browser tool is implemented by the**  `Browser class <tools.html#tools.databaseBrowser.Browser>`_


.. _interactivePlots:

Interactive Plots Maker
-----------------------

This tool allows to easily produce interactive plots which relate the SModelS output with information on the user's model stored in the SLHA files. It gives 2d plots in the parameter space defined by the user, with additional user-defined information appearing in hover boxes. The output is in html format for viewing in a web browser. The aim is not to make publication-ready plots but to facilitate getting an overview of e.g. the properties of points in a scan. NB: this needs SLHA model input and SModelS python output!

**Required python packages are:** plotly, pandas, pyslha, os, decimal

**The usage of the interactive plots tool is:**

.. include:: InteractivePlots.rst


A typical
usage example is: ::

   smodelsTools.py interactive-plots -f inputFiles/scanExample/smodels-output/ -s inputFiles/scanExample/slha -p iplots_parameters.py -o results/iplots/

which will produce 3x9 plots in the gluino vs squark mass plane from a small scan example, viewable in a web browser.


iplots parameters file
^^^^^^^^^^^^^^^^^^^^^^

The options for the interactive plots tool are defined in a parameters file, *iplots_parameters.py* in the above example.  
An example file, including all available parameters together with a short description, is stored in :download:`iplots_parameters.py <images/iplots_parameters.py>`.
Since the plotting information is model dependent, there is no default setting -- the iplots parameters file is mandatory input. 
Below we give more detailed information about each entry in this file.

* *plot_title*: main overall title for your plots, typically the model name.

* *x and y axes*: SLHA block and PDG code number of the variables you want to plot, e.g. 'm_gluino': ['MASS', 1000021].

  * **variable_x**: In a dictionary form, give the name of the x-axis variable, and the block and PDG code number to find it in the SLHA file. Example: variable_x = {'m_gluino[GeV]': ['MASS', 1000021]}. 
  * **variable_y**: same for the y-axis. Example: variable_y = {'m_suR[GeV]': ['MASS', 2000002]}

* *spectrum hover information*: defines which information from the input SLHA file will appear in the hover box. The syntax is again a python dictonary.

  * **slha_hover_information**: information from the input SLHA file, e.g. model parameters or masses. Example: slha_hover_information = {'m_gluino': ['MASS', 1000021], 'm_suR': ['MASS', 2000002], 'm_LSP': ['MASS', 1000022]} 

  * **ctau_hover_information**: displays the mean decay length in meter for the listed particle(s). Example: ctau_hover_information = {'ctau_chi1+': 1000024}

  * **BR_hover_information**: defines for which particle(s) to display decay channels and branching ratios. Example: BR_hover_information = {'BR_gluino': 1000021}. **WARNING:** Lists of branching ratios can be very long, so the may not fit in the hover box. One can define the number of entries with **BR_get_top**, e.g. BR_get_top = 5 (default: BR_get_top = 'all').

* *SModelS hover information*: defines, as a list of keywords, which information to display from the SModelS output. Example: smodels_hover_information = ['SmodelS_status', 'r_max', 'Tx', 'Analysis', 'file']. The options are:

  * **SmodelS_status**: prints whether the point is excluded or not by SModelS
 
  * **r_max**: shows the highest r-value for each parameter point

  * **chi2**: shows the chi^2 value, if available (if not, the output is 'none')

  * **Tx**: shows the topology/ies which give r_max

  * **Analysis**: shows the experimental analysis from which the strongest constraint (r_max) comes from

  * **MT_max**: shows the missing topology with the largest cross section (in SModelS bracket notation)

  * **MT_max_xsec**: shows the cross section of MT_max

  * **MT_total_xsec**: shows the total missing cross section (i.e. the sum of all missing topologies cross sections)  

  * **MT_long_xsec**: shows the total missing cross section in long cascade decays  

  * **MT_asym_xsec**: shows the total missing cross section in decays with asymmetric branches 

  * **MT_outgrid_xsec**: shows the total missing cross section outside the mass grids of the experimental results

  * **file**: shows the name of the input spectrum file 

* *Choice of plots to make*

  * **plot_data**: which points you want to plot; the options are: all, non-excluded, excluded points. Example: plot_data = ['all', 'non-excluded', 'excluded'] 

  * **plot_list**: which quantities to plot in the x,y plane; the same options as for SModels hover information apply. Example: plot_list = ['SmodelS_status','r_max', 'chi2', 'Tx', 'Analysis', 'MT_max', 'MT_max_xsec', 'MT_total_xsec', 'MT_long_xsec', 'MT_asym_xsec']  


.. _permissionsFixer:

File Permissions Fixer
----------------------

In case the software was installed under a different user than it is used
(as is the case for system-wide installs), we ship a simple tool that fixes 
the file permissions for the cross section calculation code.

**The usage of the permissions fixer is:**

.. include:: FixPermissions.rst

Execute the command as root, i.e.: ::

   sudo smodelsTools.py fixpermissions

.. _toolBox:

ToolBox
-------

As a quick way to show the status of all external tools, use 
**the toolbox:**

.. include:: ToolBox.rst

