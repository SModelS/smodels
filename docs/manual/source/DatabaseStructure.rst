.. index:: Database Structure

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
.. |ExpRes| replace:: :ref:`Experimental Result <ExpResult>`
.. |ExpRess| replace:: :ref:`Experimental Results <ExpResult>`
.. |Dataset| replace:: :ref:`DataSet<DataSet>`
.. |Datasets| replace:: :ref:`DataSets<DataSet>`
.. |Database| replace:: :ref:`Database <Database>`
.. |element| replace:: :ref:`element <element>`
.. |elements| replace:: :ref:`elements <element>`
.. |particles| replace:: :ref:`particles <particleClass>`
.. |particle| replace:: :ref:`particle <particleClass>`
.. |bracket notation| replace:: :ref:`bracket notation <bracketNotation>`


.. _databaseStruct:

Database of Experimental Results
================================

SModelS stores all the information about the experimental results in the
|Database|.
Below we describe both the :ref:`directory <folderStruct>`, :ref:`object <objStruct>` structure of the  |Database|
and :ref:`how the information in stored in the database is used within SModelS <interpolationDB>`.

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

.. _datasetUL:

Data Set Folder: Upper Limit Type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since |ULrs| have a single dataset (see |Datasets|), the info file only holds
some trivial information, such as the type of |ExpRes| (UL) and the dataset id
(None for UL-type results). Here is the content of CMS-SUS-12-024/data/dataInfo.txt as an
example:

.. literalinclude:: /literals/dataInfo.txt
   :lines: 1-2

Data Set Folder: Efficiency Map Type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For |EMrs| the ``dataInfo.txt`` contains relevant information, such as an id to
identify the |dataset| (signal region), the number of observed and expected
background events for the corresponding signal region and the respective signal
upper limits.  Here is the content of
CMS-SUS-13-012-eff/3NJet6_1000HT1250_200MHT300/dataInfo.txt as an example:

.. literalinclude:: /literals/dataInfo-eff.txt
   :lines: 1-7

.. _txnameFile:

TxName Files
^^^^^^^^^^^^

Each |DataSet| contains one or more ``TxName.txt`` file storing
the bulk of the experimental result data.
For |ULrs|, the TxName file contains the UL maps for a given simplified model
(|element| or sum of |elements|), while for |EMrs| the file contains
the simplified model efficiencies.
In addition, the TxName files also store some meta information, such
as the source of the data and the *type* of result (*prompt* or *displaced*).
If not specified, the type will be assumed to be prompt.\ [#f1]_
For instance, the first few lines of CMS-SUS-12-024/data/T1tttt.txt read:

.. literalinclude:: /literals/T1tttt.txt
   :lines: 1-8

As seen above, the first block of data in the 
file contains information about the |element|
or simplified model ([[['t','t']],[['t','t']]])
in |bracket notation| for which the data refers to as well as 
reference to the original data source and some additional information.
The simplified model is assumed to contain neutral BSM final states (MET signature)
and arbitrary intermediate BSM states.
For non-MET final states the finalState field must list the type of
BSM :ref:`particles <particleClass>` (see |UL| for more details).
An example from the CMS-EXO-12-026/data/THSCPM1b.txt file is shown below:

.. literalinclude:: /literals/THSCPM1b.txt
   :lines: 1,2,7,9,10


The second block of data in the  ``TxName.txt`` file contains the upper limits or efficiencies
as a function of the relevant simplified model parameters:

.. literalinclude:: /literals/T1tttt.txt
   :lines: 9-19
   
.. _widthGrid:   

As we can see, the data grid is given as a Python array with the structure:
:math:`[[\mbox{masses},\mbox{upper limit}], [\mbox{masses},\mbox{upper limit}],...]`.
For prompt analyses, the relevant parameters are usually the BSM masses, since
all decays are assumed to be prompt. On the other hand, results for long-lived
or meta-stable particles may depend on the BSM widths as well.
The width dependence can be easily included through the
following generalization:

.. math::
   [[M_1,M_2...],[M_A,M_B,...]] \to [[(M_1,\Gamma_1),(M_2,\Gamma_2)...],[(M_A,\Gamma_A),(M_B,\Gamma_B),...]]

In order to make the notation more compact, whenever the width dependence is not included,
the corresponding decay will be assumed to be prompt and an effective :ref:`lifetime reweigthing factor <dbReweighting>`
will be applied to the upper limits. For instance, a *mixed type* data grid is also allowed:

.. math::
   [\; [[M_1,(M_2,\Gamma_2)],[M_1,(M_2,\Gamma_2)]],\mbox{UL}\; ],\;\; [\; [[M_1',(M_2',\Gamma_2')],[M_1',(M_2',\Gamma_2')]],\mbox{UL'}\; ], \;\; ...

The example above represents a simplified model where the decay of the mother is prompt,
while the daughter does not have to be stable, hence the dependence on :math:`\Gamma_2`.
In this case, the :ref:`lifetime reweigthing factor <dbReweighting>`
is applied only for the mother decay.


.. _inclusiveSMS:

Inclusive Simplified Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~


If the analysis signal efficiencies or upper limits are insensitive to
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

Note that although the final state represented by "\*" is any Z\ :sub:`2`-even |particle|,
it must still correspond to a single particle, since the topology specifies a 2-body
decay for the initially produced BSM particle.
Finally, it might be useful to define even more inclusive simplified models, such
as the one in  CMS-EXO-13-006 eff/c000/THSCPM4.txt:

.. literalinclude:: /literals/THSCPM4.txt
   :lines: 1-2,11

In the above case the simplified model corresponds to an HSCP being initially produced
in association with any BSM particle which leads to a MET signature.
Note that "[\*]" corresponds to *any branch*, while ["\*"] means *any particle*:


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
The mapping is almost one-to-one, with a few exceptions.
Below we show the overall object structure  as well as the folders/files the objects
represent (click to enlarge):

.. image:: images/DatabaseObjects.png
   :width: 80%

The type of Python object (Python class, Python list,...) is shown in brackets.
For convenience, below we explicitly list the main database folders/files and
the Python objects they are mapped to:

* |Database| folder :math:`\rightarrow` `Database Class <experiment.html#experiment.databaseObj.Database>`_
* |ExpRes| folder :math:`\rightarrow` `ExpResult Class <experiment.html#experiment.expResultObj.ExpResult>`_
* |Dataset| folder :math:`\rightarrow` `DataSet Class <experiment.html#experiment.datasetObj.DataSet>`_
* ``globalInfo.txt`` file  :math:`\rightarrow` `Info Class <experiment.html#experiment.infoObj.Info>`_
* ``dataInfo.txt`` file  :math:`\rightarrow` `Info Class <experiment.html#experiment.infoObj.Info>`_
* ``Txname.txt`` file  :math:`\rightarrow` `TxName Class <experiment.html#experiment.txnameObj.TxName>`_


.. _databasePickle:


Database: Binary (Pickle) Format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the first time of instantiating the
`Database <experiment.html#experiment.databaseObj.Database>`_
class, the text files in *<database-path>* are loaded and parsed, and the
corresponding data objects are built. The efficiency and upper limit maps
themselves are subjected to standard preprocessing steps such as a principal
component analysis and Delaunay triangulation (see :ref:`below <interpolationDB>`).
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
 %\caption{Delaunay triangulation of an upper limit map with three mass                                                                        %parameters. The colors show the upper limit values.}

* The pickle file is created by the `createBinaryFile method <experiment.html#experiment.databaseObj.Database.createBinaryFile>`_

.. _interpolationDB:

Database: Data Processing
-------------------------

All the information contained in the :ref:`database files <folderStruct>`
is stored in the :ref:`database objects <objStruct>`.
Within SModelS the information in the |Database| is mostly used for constraining
the simplified models generated by the :ref:`decomposition <decomposition>` of the input model.
Each simplified model (or :ref:`element <element>`) generated is compared to the
simplified models contrained by the database and specified by the *constraint* and *finalStates* entries
in the  :ref:`TxName files <txnameFile>`.
The comparison allows to identify which results can be used to test the input model.
Once a matching result is found the upper limit or efficiency must be computed
for the given input |element|. As :ref:`described above <txnameFile>`, the
upper limits or efficiencies are provided as function of masses and widths in the form
of a discrete grid.
In order to compute values for any given input |element|, the data has to be
processed as decribed below.



The efficiency and upper limit maps are subjected to a few
standard preprocessing steps.
First all the units are removed, the shape of the grid is stored and
the relevant width dependence is identified (see :ref:`discussion above <txnameFile>`).
Then the masses and widths are transformed into a flat array: 

.. _dataTransf:

.. math::
   [[M_1,(M_2,\Gamma_2)],[M_A,(M_B,\Gamma_B)]] \to [M_1,M_2,M_A,M_B,\log(1+\Gamma_2),\log(1+\Gamma_B)]

 
Finally a principal component analysis and Delaunay triangulation (see :ref:`figure below <delaunay>`)
is applied over the new coordinates.
The simplices defined during triangulation are then used for linearly interpolating
the transformed data grid, thus allowing SModelS to compute efficiencies or upper limits
for arbitrary mass and width values (as long as they fall inside the data grid).
As seen above, 
the width parameters are taken logarithmically before interpolation, which
effectively corresponds to an exponential interpolation.
If the data grid does not explicitly provide a dependence on all the widths
(as in the :ref:`example above <dataTransf>`), the computed upper limit or efficiency
is then reweighted imposing the requirement of prompt
decays (see :ref:`lifetime reweighting <dbReweighting>` for more details). 
This procedure provides an efficient and numerically robust way of dealing with
generic data grids, including arbitrary parametrizations of the mass parameter
space, irregular data grids and asymmetric branches.

.. _delaunay:

.. image:: images/delaunay.png

..
 %\caption{Delaunay triangulation of an upper limit map with three mass                                                                        %parameters. The colors show the upper limit values.}



.. _dbReweighting:

Lifetime Reweighting
^^^^^^^^^^^^^^^^^^^^

From v2.0 onwards SModelS allows to include width dependent efficiencies and upper limits.
However most experimental results do not provide upper limits (or efficiencies) as a function
of the BSM particles' widths, since usually all the decays are assumed to be prompt
and the last BSM particle appearing in the cascade decay is assumed to be stable.\ [#f2]_
In order to apply these results to models which may contain meta-stable
particles, it is possible to approximate the dependence on the widths for the case in which
the experimental result requires all BSM decays to be prompt and the last BSM particle to be stable or decay *outside* the dector. 
In SModelS this is done through a reweighting factor which corresponds to the fraction
of prompt decays (for intermediate states) and decays *outside* the detector (for final BSM states)
for a given set of widths.
For instance, asumme an |EMr| only provides efficiencies (:math:`\epsilon_{prompt}`)
for prompt decays:

.. _widthExample:

.. image:: images/elementC.png
   :width: 45%


Then, for other values of the widths, an effective efficiency (:math:`\epsilon_{eff}`) can be
approximated by:

.. math::

    \epsilon_{eff} = \xi \times \epsilon_{prompt} \mbox{ , where }\xi = \mathcal{F}_{prompt} \left( \Gamma_{X_1} \right) \times \mathcal{F}_{prompt} \left( \Gamma_{X_2} \right) \times \mathcal{F}_{long} \left( \Gamma_{Y_1} \right) \times \mathcal{F}_{long} \left( \Gamma_{Y_2} \right)

In the expression above :math:`\mathcal{F}_{prompt}(\Gamma)` is the probability for the decay to be prompt 
given a width :math:`\Gamma` and :math:`\mathcal{F}_{long}(\Gamma)` is the probability for the decay to
take place *outside* the detector.
The precise values of :math:`\mathcal{F}_{prompt}` and :math:`\mathcal{F}_{long}` 
depend on the relevant detector size (:math:`L`), particle mass (:math:`M`), boost
(:math:`\beta`) and width (:math:`\Gamma`), thus
requiring a Monte Carlo simulation for each input model. Since this is not
within the spirit of the simplified model approach, we approximate the prompt and
long-lived probabilities by:

.. math::
   \mathcal{F}_{long} = \exp\left(- \frac{\Gamma L_{outer}}{\langle \gamma \beta \rangle}\right) \mbox{ and } 
   \mathcal{F}_{prompt} = 1 - \exp\left(- \frac{\Gamma L_{inner}}{\langle \gamma \beta \rangle}\right),

where :math:`L_{outer}` is the effective size of the detector (which we take to be 10 m for both ATLAS
and CMS), :math:`L_{inner}` is the effective radius of the inner detector (which we take to be 1 mm for both ATLAS
and CMS). Finally, we take the effective time dilation factor to be  :math:`\langle \gamma \beta \rangle = 1.3` when
computing :math:`\mathcal{F}_{prompt}` and :math:`\langle \gamma \beta \rangle = 1.43` when computing :math:`\mathcal{F}_{long}`.
We point out that the above approximations are irrelevant if :math:`\Gamma` is very large (:math:`\mathcal{F}_{prompt} \simeq 1`
and :math:`\mathcal{F}_{long} \simeq 0`) or close to zero (:math:`\mathcal{F}_{prompt} \simeq 0`
and :math:`\mathcal{F}_{long} \simeq 1`). Only elements containing particles which have a considerable fraction of displaced
decays will be sensitive to the values chosen above.
Also, a precise treatment of lifetimes is possible if the experimental result
(or a theory group) explicitly provides the efficiencies as a function of the widths, as :ref:`discussed above <widthGrid>`.



The above expressions allows the generalization of the efficiencies computed assuming
prompt decays to models with meta-stable particles. 
For |ULrs| the same arguments apply with one important distinction.
While efficiencies are reduced for displaced decays (:math:`\xi < 1`), upper limits are enhanced, since they
are roughly inversely proportional to signal efficiencies. Therefore, for |ULrs|, we have:

.. math::

    \sigma_{eff}^{UL} = \sigma_{prompt}^{UL}/\xi


Finally, we point out that for the experimental results which provide 
efficiencies or upper limits as a function of some (but not all) BSM widths appearing
in the simplified model (see the :ref:`discussion above <widthGrid>`), 
the reweighting factor :math:`\xi` is computed using only the widths not present
in the grid.


.. [#f1] Prompt results are all those which assumes all decays to be prompt and the last BSM particle to be stable (or decay outside the detector).
       Searches for heavy stable charged particles (HSCPs), for instance, are classified as *prompt*, since the HSCP is assumed to decay
       outside the detector. Displaced results on the other hand require at least one decay to take place inside the detector.

.. [#f2] An obvious exception are searches for long-lived particles with displaced decays.

