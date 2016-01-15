.. index:: Database Definitions

.. _databaseDefs:

Database Definitions
====================

The SModelS database collects experimental results from both ATLAS and CMS, which are
used to compute the experimental constraints on specific models.
Starting with version 1.2, the SModelS database includes two types of experimental constraints:

*  Upper Limit (UL) constraints: upper limit constraints on :math:`\sigma \times BR` of simplified models, 
provided by the experimental collaborations (see :ref:`UL Constraints<ULconstraint>`);
* Efficiency Map (EM) constraints: constraints on :math:`\sigma \times BR \times \epsilon` of simplified models,
either provided by the experimental collaborations or computed by theory 
collaborations (see :ref:`EM Constraints<EMconstraint>`); 

Although the two types of constraints above are very distinct,
both the folder structure and the object structure of SModelS are sufficiently flexible to
simutaneously handle both UL and EM types.
Each publication or conference note can be included in the database 
as a :ref:`Experimental Result<ExpResult>` and are organized according to the folder structure described in
:ref:`Database of Experimental Results<database>`. Hence, the database is simply a collection of experimental results,
as described below. 


.. _Database

Database
--------

The database 


* **The Database is described by the** `ExpResult Class <../../../documentation/build/html/experiment.html#experiment.databaseObjects.Database>`_

.. _ExpResult:

Experimental Result
-------------------
An experimental result contains all the relevant information corresponding to an experimental
publication or conference note. In particular it holds  general information about the experimental analysis,
such as the corresponding luminosity, center of mass energy, publication reference, etc.
The current version allows for two possible types of experimental results: one containing upper limit maps
and one containing efficiency maps.


.. _ULtype:

Upper Limit Type
^^^^^^^^^^^^^^^^

Upper Limit (UL) experimental results contains experimental constraints on the cross-section times branching ratio
( :math:`\sigma \times BR` ) from a specific experimental publication or conference note.
Each UL  result corresponds to the 95% upper limit constraints on :math:`\sigma \times BR` for a given 
:ref:`element <element>` or sum of :ref:`elements <element>`.
 In some cases it is useful to define a shorthand notation for the  :ref:`element <element>` 
 or sum of :ref:`elements <element>` constrained by the  result. SModelS uses TxNames

For illustration, consider this CMS example:

.. _ULplot:

.. image:: images/ULexample.png
   :height: 480px

In this case the UL result constrains the element :math:`[[[jet]],[[jet]]]`, where we are using the notation
defined in :ref:`Bracket Notation <bracketnotation>`.


Each individual UL result holds the upper limit values (*upper limit map*)
on :math:`\sigma \times BR` as a function of the respective parameter space (usually BSM masses or slices over mass planes). 
Furthermore, the corresponding :ref:`constraints <ULconstraint>`
and :ref:`conditions <ULconditions>` must also be specified.
Usually one :ref:`Experimental Result<ExpResult>` contains upper limit maps for more than one :ref:`element <element>`.


As mentioned abovem, UL results also contain information about the analysis luminosity, center-of-mass, 
publication reference and others.
*We also point out that the exclusion curve is never used by SModelS*.

Note that a given experimental publication (or conference note) may contain several UL analyses, since a single
publication may contain upper limits for several different :ref:`elements <element>` (or :ref:`constraints <ULconstraint>`).

.. _EMtype:

Efficiency Map Type
^^^^^^^^^^^^^^^^^^^


For illustration, consider this CMS example:

.. _EMplot:

.. image:: images/ULexample.png
   :height: 480px


* **Experimental Results are described by the** `ExpResult Class <../../../documentation/build/html/experiment.html#experiment.databaseObjects.ExpResult>`_

.. _DataSet:

Data Set
--------



* **Data Sets are described by the** `DataSet Class <../../../documentation/build/html/experiment.html#experiment.datasetObject.DataSet>`_






.. _ULconstraint:

Analysis Constraints
^^^^^^^^^^^^^^^^^^^^

Constraints are defined as the :ref:`element <element>` or sum over :ref:`elements <element>`
which is constrained by the experimental upper limits. The analysis constraints can also be simply expressed in the 
:ref:`bracket notation <bracketnotation>` as a sum of individual elements.

As an example, consider the `ATLAS analysis <https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-049/>`_ shown below:

.. _constraintplot:

.. image:: images/constraintExample.png
   :height: 580px

As we can see, the upper limits apply to the sum of the cross-sections:

.. math::
    \sigma = \sigma([[[e^+]],[[e^-]]]) + \sigma([[[\mu^+]],[[\mu^-]]])
    
In this case the analysis constraint is simply:

.. math::
    [[[e^+]],[[e^-]]] + [[[\mu^+]],[[\mu^-]]]
    
where it is understood that the sum is over the weights of the respective elements
and not over the elements themselves.    
    
    
Note that the sum can be over particle charges, flavors or more complex combinations of elements.
However, *almost all analyses sum only over elements sharing a common* :ref:`topology <topology>`.

.. _ULconditions:

Analysis Conditions
^^^^^^^^^^^^^^^^^^^

When the analysis :ref:`constraints <ULconstraint>` are non-trivial (refer to a sum of elements), it is often the case
that there are implicit (or explicit) assumptions about the contribution of each element. For instance,
in the :ref:`figure above <constraintplot>`, it is implicitly assumed that each lepton flavor contributes equally
to the summed cross-section:

.. math::    
    \sigma([[[e^+]],[[e^-]]]) = \sigma([[[\mu^+]],[[\mu^-]]])           \;\;\; \mbox{(condition)}
    

Therefore, when applying these constraints to general models, one must also verify if
these conditions are satisfied. Once again we can express these conditions in 
:ref:`bracket notation <bracketnotation>`:

.. math::    
    [[[e^+]],[[e^-]]] = [[[\mu^+]],[[\mu^-]]]           \;\;\; \mbox{(condition)}

where it is understood that the condition refers to the weights of the respective elements
and not to the elements themselves.

In several cases it is desirable to relax the analysis conditions, so the analysis
upper limits can be applied to a broader spectrum of models. Once again, for the example mentioned
above, it might be reasonable to impose instead:

.. math::
    [[[e^+]],[[e^-]]] \simeq [[[\mu^+]],[[\mu^-]]]           \;\;\; \mbox{(fuzzy condition)}

The *departure* from the exact condition can then be properly quantified and one can decide whether the analysis upper limits are applicable or not to the model being considered.
Concretely, for each condition a number between 0 and 1 is returned, 
where 0 means the condition is exactly satisfied and 1 means it is maximally violated.
Allowing for a :math:`20\%` violation of a condition corresponds approximately to 
a ''condition violation value'' (or simply condition value) of 0.2.
The condition values  are given as an output of SModelS, so the user can decide what are the
maximum acceptable values.



.. _EManalysis:

Efficiency Map Analyses
-----------------------

Efficiency Map (EM) analyses are more fundamental than :ref:`UL analyses <ULanalysis>`. Instead of holding cross-section upper limits, they correspond to one or more :ref:`efficiency maps <effmap>` together with the 
information about the expected and observed data for the relevant signal region(s).::


   Note: Efficiency Map analyses are not yet functional in the public release!!!
   


.. _effmap:

Efficiency Maps
^^^^^^^^^^^^^^^

Efficiency maps correspond to a grid of simulated acceptance times efficiency 
( :math:`A \times \epsilon` ) values for specific signal region(s). In the following we will refer to :math:`A \times \epsilon` simply as *efficiency*.  

The signal is assumed to correspond to a single element, which characterizes the basic signal kinematics
and hence its efficiency.
The efficiency grid is usually a function of the BSM masses appearing in the element, as shown by the example below:

.. _EMplot:

.. image:: images/EMexample.png
   :height: 480px

Although efficiency maps are most useful for :ref:`EM analyses <EManalysis>`, they can also be constructed for
:ref:`UL analyses <ULanalysis>`. For the latter, the efficiencies for a given element are either 1, if the element
belongs to the :ref:`UL analysis constraint <ULconstraint>`, or 0, if the element
does not belong to the :ref:`UL analysis constraint <ULconstraint>`.
