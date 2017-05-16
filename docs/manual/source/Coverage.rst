.. index:: Missing Topologies

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


.. _topCoverage:

Topology Coverage
=================


The constraints provided by SModelS are obviously limited
by its |database| and the available set of simplified model interpretations
provided by the experimental collaborations or computed by theory groups.
Therefore it is interesting to identify classes of missing simplified models
(or missing topologies) which are relevant for a given input model, but are
not constrained by the SModelS |database|. This task is performed
as a last step in SModelS, once the |decomposition| and the |theory predictions|
have been computed.

Given the |decomposition| output (list of |elements|), as well as the |database|
information, it finds and classifies the |elements| which are
not tested by any of the |results| in the |database|.
These elements are grouped into the following classes:

* *missingTopos*: |elements| which are not tested by any of the |results| in the |database| (independent of the element mass). The missing topologies are further classified as:

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
