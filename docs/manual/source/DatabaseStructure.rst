.. index:: Database Structure

.. _databaseStruct:

.. |constraint| replace:: :ref:`constraint <ULconstraint>`
.. |conditions| replace:: :ref:`conditions <ULconditions>` 
.. |fb-1| replace:: :math:`\mathrm{fb}^{-1}`
.. |sqrts| replace:: :math:`\sqrt{s}`
.. |EM| replace:: :ref:`EM-type <EMtype>`
.. |UL| replace:: :ref:`UL-type <ULtype>`
.. |EMr| replace:: :ref:`EM-type result <EMtype>`
.. |ULr| replace:: :ref:`UL-type result <ULtype>`
.. |EMrs| replace:: :ref:`EM-type results <EMtype>`
.. |ULrs| replace:: :ref:`UL-type results <ULtype>`
.. |ExpRes| replace:: :ref:`Experimental Result<ExpResult>`
.. |ExpRess| replace:: :ref:`Experimental Results<ExpResult>`
.. |Dataset| replace:: :ref:`DataSet<DataSet>`
.. |Datasets| replace:: :ref:`DataSets<DataSet>`
.. |Database| replace:: :ref:`Database <Database>`
.. |element| replace:: :ref:`element <element>`
.. |elements| replace:: :ref:`elements <element>`
.. |bracket notation| replace:: :ref:`bracket notation <bracketNotation>`

Database of Experimental Results
================================

SModelS stores all the information about the experimental results in the 
|Database|. 
Below we describe both the :ref:`directory <folderStruct>` and :ref:`object <objStruct>` structure of the  |Database|.

.. _folderStruct:

Database: Directory Structure
-----------------------------

The :ref:`Database <Database>` is organized as files in an ordinary (UNIX)
directory hierarchy, with a thin Python layer serving as the access to the
database.  The overall structure of the directory hierarchy and its contents is
depicted in the scheme below (click to enlarge):

.. image:: images/DatabaseFolders.png
   :width: 80%

As seen above, the top level of the SModelS database categorizes the analyses
by LHC center-of-mass energies, |sqrts|:

* 8 TeV
* 13 TeV

Also, the top level directory contains a file called ``version`` with the
version string of the database.
The second level splits the results up between the different experiments:

* 8TeV/CMS/
* 8TeV/ATLAS/

The third level of the directory hierarchy encodes the |ExpRess|:

* 8TeV/CMS/CMS-SUS-12-024
* 8TeV/ATLAS/ATLAS-CONF-2013-047
* ...


* **The Database folder is described by the** `Database Class <experiment.html#experiment.databaseObj.Database>`_

Experimental Result Folder
^^^^^^^^^^^^^^^^^^^^^^^^^^

Each |ExpRes| folder contains: 

* a folder for each |Dataset| (e.g. ``data``)
* a ``globalInfo.txt`` file

The ``globalInfo.txt`` file contains the meta information about the |ExpRes|.
It defines the center-of-mass energy |sqrts|, the integrated luminosity, the id
used to identify the result and additional information about the source of the
data.  Here is the content of CMS-SUS-12-024/globalInfo.txt as an example:
      
.. literalinclude:: /literals/globalInfo.txt
   :lines: 1-11

* **Experimental Result folder is described by the** `ExpResult Class <experiment.html#experiment.expResultObj.ExpResult>`_
* **globalInfo files  are descrived by the** `Info Class <experiment.html#experiment.infoObj.Info>`_

Data Set Folder
^^^^^^^^^^^^^^^

Each |Dataset| folder (e.g. ``data``) contains:

* the Upper Limit maps for |ULrs| or Efficiency maps for |EMrs| (``TxName.txt`` files)
* a ``dataInfo.txt`` file containing meta information about the |Dataset|

* **Data Set folders are  described by the** `DataSet Class <experiment.html#experiment.datasetObj.DataSet>`_
* **TxName files are described by the** `TxName Class <experiment.html#experiment.txnameObj.TxName>`_
* **dataInfo files are described by the** `Info Class <experiment.html#experiment.infoObj.Info>`_

Data Set Folder: Upper Limit Type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since |ULrs| have a single dataset (see |Datasets|), the info file only holds
some trivial information, such as the type of |ExpRes| (UL) and the dataset id
(None for UL-type results). Here is the content of CMS-SUS-12-024/data/dataInfo.txt as an
example:

.. literalinclude:: /literals/dataInfo.txt
   :lines: 1-2

For |ULrs|, each ``TxName.txt`` file contains the UL map for a given simplified model
(|element| or sum of |elements|) as well as some meta information,
including the corresponding |constraint| and |conditions|.  The
first few lines of CMS-SUS-12-024/data/T1tttt.txt read:

.. literalinclude:: /literals/T1tttt.txt
   :lines: 1-8

If the finalState property is not provided, the simplified model is assumed to 
contain neutral BSM final states in each branch, leading to a MET signature.
However, if this is not the case, the non-MET final states must be explicitly listed
in the  ``TxName.txt`` file (see :ref:`final state classes <final stateOdd>` for more details).
An example from the CMS-EXO-12-026/data/THSCPM1b.txt file is shown below:

.. literalinclude:: /literals/THSCPM1b.txt
   :lines: 1,2,7,9,10

   
The second block of data in the  ``TxName.txt`` file contains the upper limits as a function of the BSM masses:

.. literalinclude:: /literals/T1tttt.txt
   :lines: 9-19

As we can see, the UL map is given as a Python array with the structure: 
:math:`[[\mbox{masses},\mbox{upper limit}], [\mbox{masses},\mbox{upper limit}],...]`.




Data Set Folder: Efficiency Map Type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For |EMrs| the ``dataInfo.txt`` contains relevant information, such as an id to
identify the |dataset| (signal region), the number of observed and expected
background events for the corresponding signal region and the respective signal
upper limits.  Here is the content of
CMS-SUS-13-012-eff/3NJet6_1000HT1250_200MHT300/dataInfo.txt as an example:

.. literalinclude:: /literals/dataInfo-eff.txt
   :lines: 1-7

For |EMrs|, each ``TxName.txt`` file contains the efficiency map for a given
simplified model (|element| or sum of |elements|) as well as some meta
information.
Here is the first few lines of CMS-SUS-13-012-eff/3NJet6_1000HT1250_200MHT300/T2.txt:

.. literalinclude:: /literals/T2.txt
   :lines: 1-8
   
As seen above, the first block of data in the ``T2.txt`` file contains
information about the |element| (:math:`[[[\mbox{jet}]],[[\mbox{jet}]]]`) 
in |bracket notation| for which the
efficiencies refers to as well as reference to the original data source and
some additional information.
As in the Upper Limit case, the simplified
model is assumed to contain neutral BSM final states (MET signature).
For non-MET final states the  finalState field must list
the :ref:`final state signatures <final stateOdd>`.
The second block of data contains the efficiencies as a function of the BSM masses:

.. literalinclude:: /literals/T2.txt
   :lines: 9-15

As we can see the efficiency map is given as a Python array with the structure: 
:math:`[[\mbox{masses},\mbox{efficiency}], [\mbox{masses},\mbox{efficiency}],...]`.


.. _dbReweighting:

Data Set Folder: Lifetime reweighting factor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to the efficiency maps, there are also maps containing the probability for each simpilified model to decay in this way.
Since v2.0.0 SModelS is able to handle different lifetime dependent signatures. Therefore whenever an element is matched 
to an experimental result it has to be multiplied by the probability to decay in the same combination of prompt or displaced decays and long-lived particles.

The probabilities are computed by: 

.. math::
	\mathcal{F}(\tau) = \frac{(A \times \epsilon)_{\tau}}{(A \times \epsilon)_{\infty}}
	
	
The lifetime reweighting map is given as a Python array with the structure: 
:math:`[[\mbox{masses},\mbox{lifetime},\mbox{efficiency}], [\mbox{masses},\mbox{lifetime},\mbox{efficiency}],...]`.


.. _inclusiveSMS:

Inclusive Simplified Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~


If the analysis signal efficiencies are insensitive to
some of the simplified model final states, it might be convenient to define
*inclusive* simplified models. A typical case are some of the heavy stable charged
particle searches, which only rely on the presence of a non-relativistic charged
particle, which leads to  an anomalous charged track signature.
In this case the signal efficiencies are highly insensitive to the remaining event
activity and the corresponding simplified models can be very inclusive.
In order to handle this inclusive cases in the database we allow for wildcards
when specifying the constraints.
For instance, the constraint for the CMS-EXO-13-006 eff/c000/THSCPM3.txt
reads:

.. literalinclude:: /literals/THSCPM3.txt
   :lines: 1-2
   
and represents the (inclusive) simplified model:

.. image:: images/elementInclusive.png
   :width: 35%

Note that although the final state represented by "\*" is any Z\ :sub:`2`-even :ref:`final states <final statesEven>`,
it must still correspond to a single particle, since the topology specifies a 2-body 
decay for the initially produced BSM particle.
Finally, it might be useful to define even more inclusive simplified models, such
as the one in  CMS-EXO-13-006 eff/c000/THSCPM4.txt:

.. literalinclude:: /literals/THSCPM4.txt
   :lines: 1-2,11

In the above case the simplified model corresponds to an HSCP being initially produced
in association with any BSM particle which leads to a MET signature.
Notice that the notation "[\*]" corresponds to *any `branch*, while ["\*"] means *any particle*:


.. image:: images/elementInclusive2.png
   :width: 35%

In such cases the mass array for the arbitrary branch must also be specified as
using wildcards:

.. literalinclude:: /literals/THSCPM4.txt
   :lines: 12-14



.. _objStruct:

Database: Object Structure
--------------------------

The :ref:`Database  folder structure <folderStruct>` is mapped to Python
objects in SModelS.
The mapping is almost one-to-one, except for a few exceptions.
Below we show the overall object structure  as well as the folders/files the objects
represent (click to enlarge):

.. image:: images/DatabaseObjects.png
   :width: 80%
   
The type of Python object (Python class, Python list,...) is shown in brackets.
For convenience, below we explicitly list the main database folders/files and
the Python objects they are mapped to:

* |Database| folder :math:`\rightarrow` `Database Class <experiment.html#experiment.databaseObj.Database>`_
* |ExpRes| folder :math:`\rightarrow` `ExpResult Class <experiment.html#experiment.databaseObj.ExpResult>`_
* |Dataset| folder :math:`\rightarrow` `DataSet Class <experiment.html#experiment.datasetObj.DataSet>`_
* ``globalInfo.txt`` file  :math:`\rightarrow` `Info Class <experiment.html#experiment.infoObj.Info>`_
* ``dataInfo.txt`` file  :math:`\rightarrow` `Info Class <experiment.html#experiment.infoObj.Info>`_
* ``Txname.txt`` file  :math:`\rightarrow` `TxName Class <experiment.html#experiment.txnameObj.TxName>`_


.. _databasePickle:

Database: Binary (Pickle) Format
--------------------------------

At the first time of instantiating the 
`Database <experiment.html#experiment.databaseObj.Database>`_
class, the text files in *<database-path>*.
are loaded and parsed, and the corresponding                                                    
data objects are built. The efficiency and upper limit maps themselves are                                                                    
subjected to standard preprocessing steps such as a principal component                                                                       
analysis and Delaunay triangulation (see Figure below).
The simplices defined during triangulation are then used for linearly interpolating the data grid,                                            
thus allowing SModelS to compute efficiencies or upper limits for arbitrary                                                                   
mass values (as long as they fall inside the data grid).                                                                                      
This procedure provides an efficient and numerically robust way of                                                                            
dealing with generic data grids, including arbitrary parametrizations of the mass parameter space,                                            
irregular data grids and asymmetric branches.                                                                                                 
                                                                                                                                              
.. image:: images/delaunay.png

..
 %\caption{Delaunay triangulation of an upper limit map with three mass                                                                        %parameters. The colors show the upper limit values.}

For the sake of efficiency, the entire database -- including the Delaunay                                                                     
triangulation -- is then serialized into a pickle                                                                                             
file (*<database-path>/database.pcl*), which will be read directly the next time the database is loaded.                                                                       
If any changes in the database folder structure are detected, the python or the SModelS                                                       
version has changed, SModelS will automatically re-build the pickle file. This                                                                
action may take a few minutes, but it is again performed only once.                                                                           
If desired, the pickling process can be skipped using the option *force_load = `txt'*
in the constructor of
`Database <experiment.html#experiment.databaseObj.Database>`_ .  


..
 Due to the large number of experimental results contained in the SModelS
 |Database|, parsing the :ref:`database folders <folderStruct>` and building the
 corresponding :ref:`database objects <objStruct>` may require a non-negligible
 CPU time. In some cases this may be the most time consuming task when
 testing a single input file.  Furthermore this procedure does not have to be
 repeated every time SModelS is run.
 In order to avoid these issues, SModelS serializes the 
 `database object <experiment.html#experiment.databaseObj.Database>`_
 into a pickle file (*<database-path>/database.pcl*), which can then be read
 directly when loading the database.
 Since reading the pickle file is much faster than parsing the :ref:`database folders <folderStruct>`,
 there is a considerable speed improvement when using the pickle file.
 If any changes in the :ref:`database folder structure <folderStruct>` 
 are detected or the SModelS version has changed,
 SModelS will automatically re-build the pickle file.
 This action may take a few minutes, but it is only performed once.
 SModelS automatically builds (if necessary) and loads the binary database when a 
 `Database object <experiment.html#experiment.databaseObj.Database>`_
 is created. Nonetheless, the user can enforce loading (parsing) the *text
 database* using the option *force_load = 'txt'* in the constructor of
 `Database <experiment.html#experiment.databaseObj.Database>`_ .  

* The pickle file is created by the `createBinaryFile method <experiment.html#experiment.databaseObj.Database.createBinaryFile>`_ 
