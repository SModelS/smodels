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


* **The Database folder is described by the** `Database Class <../../../documentation/build/html/experiment.html#experiment.databaseObj.Database>`_

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

* **Experimental Result folder is described by the** `ExpResult Class <../../../documentation/build/html/experiment.html#experiment.expResultObj.ExpResult>`_
* **globalInfo files  are descrived by the** `Info Class <../../../documentation/build/html/experiment.html#experiment.infoObj.Info>`_

Data Set Folder
^^^^^^^^^^^^^^^

Each |Dataset| folder (e.g. ``data``) contains:

* the Upper Limit maps for |ULrs| or Efficiency maps for |EMrs| (``TxName.txt`` files)
* a ``dataInfo.txt`` file containing meta information about the |Dataset|

* **Data Set folders are  described by the** `DataSet Class <../../../documentation/build/html/experiment.html#experiment.datasetObj.DataSet>`_
* **TxName files are described by the** `TxName Class <../../../documentation/build/html/experiment.html#experiment.txnameObj.TxName>`_
* **dataInfo files are described by the** `Info Class <../../../documentation/build/html/experiment.html#experiment.infoObj.Info>`_

Data Set Folder: Upper Limit Type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since |ULrs| have a single dataset (see |Datasets|), the info file only holds
some trivial information, such as the type of |ExpRes| (UL) and the dataset id
(None). Here is the content of CMS-SUS-12-024/data/dataInfo.txt as an
example:

.. literalinclude:: /literals/dataInfo.txt
   :lines: 1-2

For |ULrs|, each ``TxName.txt`` file contains the UL map for a given simplified model
(|element| or sum of |elements|) as well as some meta information,
including the corresponding |constraint| and |conditions|.  The
first few lines of CMS-SUS-12-024/data/T1tttt.txt read:

.. literalinclude:: /literals/T1tttt.txt
   :lines: 1-8
   
The second block of data contains the upper limits as a function of the BSM masses:

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
The second block of data contains the efficiencies as a function of the BSM masses:

.. literalinclude:: /literals/T2.txt
   :lines: 9-15

As we can see the efficiency map is given as a Python array with the structure: 
:math:`[[\mbox{masses},\mbox{efficiency}], [\mbox{masses},\mbox{efficiency}],...]`.

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

* |Database| folder :math:`\rightarrow` `Database Class <../../../documentation/build/html/experiment.html#experiment.databaseObj.Database>`_
* |ExpRes| folder :math:`\rightarrow` `ExpResult Class <../../../documentation/build/html/experiment.html#experiment.databaseObj.ExpResult>`_
* |Dataset| folder :math:`\rightarrow` `DataSet Class <../../../documentation/build/html/experiment.html#experiment.datasetObj.DataSet>`_
* ``globalInfo.txt`` file  :math:`\rightarrow` `Info Class <../../../documentation/build/html/experiment.html#experiment.infoObj.Info>`_
* ``dataInfo.txt`` file  :math:`\rightarrow` `Info Class <../../../documentation/build/html/experiment.html#experiment.infoObj.Info>`_
* ``Txname.txt`` file  :math:`\rightarrow` `TxName Class <../../../documentation/build/html/experiment.html#experiment.txnameObj.TxName>`_


.. _databasePickle:

Database: Binary (Pickle) Format
--------------------------------

Due to the large number of experimental results contained in the SModelS
|Database|, parsing the :ref:`database folders <folderStruct>` and building the
corresponding :ref:`database objects <objStruct>` may require a non-negligible
CPU time. In some cases this may be the most time consuming task when
testing a single input file.  Furthermore this procedure does not have to be
repeated every time SModelS is run.

In order to avoid these issues, SModelS serializes the 
`database object <../../../documentation/build/html/experiment.html#experiment.databaseObj.Database>`_
into a pickle file (*<database-path>/database.pcl*), which can then be read
directly when loading the database.
Since reading the pickle file is much faster than parsing the :ref:`database folders <folderStruct>`,
there is a considerable speed improvement when using the pickle file.
If any changes in the :ref:`database folder structure <folderStruct>` 
are detected or the SModelS version has changed,
SModelS will automatically re-build the pickle file.
This action may take a few minutes, but it is only performed once.

SModelS automatically builds (if necessary) and loads the binary database when a 
`Database object <../../../documentation/build/html/experiment.html#experiment.databaseObj.Database>`_
is created. Nonetheless, the user can enforce loading (parsing) the *text
database* using the option *force_load = 'txt'* in the constructor of
`Database <../../../documentation/build/html/experiment.html#experiment.databaseObj.Database>`_ .  


* The pickle file is created by the `createBinaryFile method <../../../documentation/build/html/experiment.html#experiment.databaseObj.Database.createBinaryFile>`_ 
