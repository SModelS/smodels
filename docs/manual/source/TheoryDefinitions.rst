.. index:: Theory Definitions

.. |element| replace:: :ref:`element <element>`
.. |elements| replace:: :ref:`elements <element>`
.. |topology| replace:: :ref:`topology <topology>`
.. |topologies| replace:: :ref:`topologies <topology>`
.. |bracket notation| replace:: :ref:`bracket notation <bracketNotation>`
.. |final states| replace:: :ref:`final states <final statesEven>`
.. |parameters| replace:: :ref:`parameters file <parameterFile>`


.. _theoryDefs:

Simplified Model Definitions
============================

The so-called `theory module <theory.html#theory>`_ contains the basic tools necessary for decomposing the input model
(either in LHE or SLHA format) into simplified model |topologies| and using the output of the decomposition
to compute the :ref:`theoretical prediction <theoryPredictions>` for a given :ref:`experimental result <ExpResult>`.


The applicability of SModelS is currently restricted to models which contain a Z\ :sub:`2` 
symmetry (R-parity in SUSY, K-parity in UED, ...). This is required in
order to provide a clear structure for the simplified model topologies appearing
during the :ref:`decomposition <decomposition>` of the input model.
Below we describe the basic concepts and language used in SModelS
to describe the simplified model topologies.

.. _element:

Elements
--------

A simplified model topology representing a specific cascade decay of a pair of BSM states produced in
the hard scattering is called an element in the SModelS language.
Elements contain the Z\ :sub:`2`-even particles appearing in
the cascade decay and the BSM (Z\ :sub:`2`-odd) states
which have decayed or appear in the last step of the decay.
Furthermore, the last BSM (Z\ :sub:`2`-odd) particle is classified
according to its quantum numbers as a specific *final state*
class: *MET*, *HSCP*, *R-hadron*,etc. 
A representation of an element is shown below:


.. _elementscheme:

.. image:: images/elementB.png
   :width: 40%
   
An element may also hold information about its corresponding 
weight (cross section times branching ratio times efficiency). [#f1]_
The overall properties of an element are illustrated in the scheme below:

.. _topscheme:

.. image:: images/topSchemeB.png
   :width: 40%

SModelS works under the inherent assumption that, for collider purposes,
all the essential properties of a BSM model can be encapsulated by its
elements.
Such an assumption is extremely helpful to cast the theoretical predictions of a
specific BSM model in a model-independent framework, which can then be compared
against the corresponding experimental limits.
For instance, as shown in the :ref:`scheme above <elementscheme>`, only the
masses and the widths of the BSM states and the quantum numbers (color and electric charge) of the last BSM state are used, while
other properties, such as their spins are ignored (all quantum numbers are, however, stored for book-keeping).

Below we describe in more detail the element properties and their implementation
in SModelS.


* **Elements are described by the** `Element Class <theory.html#theory.element.Element>`_    


.. _vertex:

Vertices
^^^^^^^^
Each Z\ :sub:`2`-odd decay is represented by a vertex containing the outgoing states (one Z\ :sub:`2`-odd
state and the Z\ :sub:`2`-even particles), as shown in the :ref:`scheme above <topscheme>`.


.. _particleClass:

Particles
^^^^^^^^^

A particle represents any Z\ :sub:`2`-even and Z\ :sub:`2`-odd state appearing in an element.
It can hold any of the following information: Z\ :sub:`2`-parity, label (a string describing the particle, e.g. 'e-'), 
pdg number, mass, electric charge, color charge, spin, decay width and decays to other particles.
However, different sub sets of these properties are used for different particles as outlined below.

.. _final statesEven:

Z\ :sub:`2`-even Final States
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Z\ :sub:`2`-even final states coming out of a vertex (see :ref:`scheme above <topscheme>`) usually
correspond to Standard Model particles (electrons, gauge bosons, Higgs,...).
However, if the input model contains  Z\ :sub:`2`-even BSM states (such as additional Higgs bosons),
these also appear as final states.
The only information used from the final states are their labels.
In contrast, stable or long-lived Z\ :sub:`2`-odd particles which might appear in the detector (either as MET or charged tracks)
are *not* classified as final states [#f2]_ .


* Z\ :sub:`2`-even **states are defined in** smodels/share/default_particles.py 

.. _odd states:

Z\ :sub:`2`-odd Intermediate States
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The intermediate Z\ :sub:`2`-odd states are always assumed to consist of BSM particles with Z\ :sub:`2`
conserving decays of the form: (Z\ :sub:`2`-odd state) :math:`\rightarrow`  (Z\ :sub:`2`-odd state') + |final states|.
These decays can either be prompt or displaced.
The only information used from the intermediate states are their masses and widths(see :ref:`scheme above <topscheme>`).

* Z\ :sub:`2`-odd **states are defined by the input model file** (see :ref:`model <parameterFileModel>` in |parameters|)

.. _final stateOdd:

Z\ :sub:`2`-odd Final State Particles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Besides the intermediate Z\ :sub:`2`-odd BSM states, due to the assumed Z\ :sub:`2` symmetry,
the element must also contain one stable Z\ :sub:`2`-odd final state (at least
in collider scales). The quantum numbers of this BSM final state are essential for defining which
type of signature this element represents. 
The only information used from the final state particles are their electric and color charge.
In an element the  Z\ :sub:`2`-odd final state  quantum numbers are mapped to a final state particle,
as defined in the `finalStateParticles module <experiment.html#module-experiment.finalStateParticles>`_ .
Some examples of final state classes are: 'MET', 'HSCP' and 'RHadronQ'.
New final state classes can also be easily defined in this module.  


.. _branch:

Branches
^^^^^^^^

A branch is the basic substructure of an |element|.
It represents a series of cascade decays of a single initial Z\ :sub:`2`-odd
state.
The diagram below illustrates an example of a branch.

.. image:: images/branchTopB.png
   :width: 25%

The structure of each branch is fully defined by its number of vertices and the number of 
|final states| coming out of each vertex. 
Furthermore,  the branch also holds the information about the Z\ :sub:`2`-even |final states|
coming out of each vertex, the Z\ :sub:`2`-odd states
and the Z\ :sub:`2`-odd :ref:`final state particle <final stateOdd>` (e.g. 'MET'), as shown below.


.. image:: images/branchElB.png
   :width: 35%
   
* **Branches are described by the** `Branch Class <theory.html#theory.branch.Branch>`_   


.. _notation:

Element Representation: Bracket Notation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The structure and final states of |elements| are represented in textual form using a nested brackets
notation. The scheme below shows how to convert between the graphical and bracket representations of an element:


.. _bracketnotation:

.. image:: images/bracketNotationB.png
   :width: 50%

The brackets are ordered and nested in the following way. 
The outermost brackets correspond to the :ref:`branches <branch>` of the |element|.
The branches are sorted according to their size (see :ref:`element sorting <elementsorting>`) 
and each branch contains an *ordered* list of :ref:`vertices <vertex>`.
Each vertex contains a list of the  Z\ :sub:`2`-even |final states| (sorted alphabetically) coming out of the vertex.
Schematically, for the example in the :ref:`figure above <bracketnotation>`, we have::

   element = [branch1, branch2]
      branch1 = [vertex1]
         vertex1 = [l+,l-]
      branch2 = [vertex1,vertex2]
         vertex1 = [l+]
         vertex2 = [nu]

Using the above scheme it is possible to unambiguously describe each |element| with a simple list of nested brackets.
However, in order to fully specify all the information relative to a single |element|, we must
also include the list of masses and widths for the Z\ :sub:`2`-odd states, the list of Z\ :sub:`2`-odd
:ref:`final state particles <final stateOdd>` and the element weight.
The masses for the BSM (Z\ :sub:`2`-odd) states can also be represented by a mass array
for each branch, as shown below:

.. _massnotation:

.. image:: images/massNotationB.png
   :width: 65%
   
Finally the Z\ :sub:`2`-odd :ref:`final state particles <final stateOdd>` can also
be represented as a list in addition to the bracket notation:

.. _bracketnotationFull:

.. image:: images/bracketNotationB.png
   :width: 70%
   
.. _topology:

Topologies
----------

It is often useful to classify |elements| according to their
overall structure or topology.
Each topology corresponds to an *undressed*
|element|, removed of its  Z\ :sub:`2`-even
|final states|,  Z\ :sub:`2`-odd final state particles and Z\ :sub:`2`-odd masses.
Therefore the topology is fully determined by its number of
branches, number of vertices in each :ref:`branch <branch>` and number of
 Z\ :sub:`2`-even |final states| coming out of each :ref:`vertex <vertex>`.
An example of a topology is shown below:

.. image:: images/globTopB.png
   :width: 25%

Within SModelS, elements are grouped according to their
topology. Hence  topologies represent a list of elements sharing a
common basic structure (same number of branches, vertices and
final states in each vertex).

* **Topologies are described by the** `Topology Class <theory.html#theory.topology.Topology>`_   

.. [#f1] In order to treat the UL and EM map results on the same footing,
   SModelS applies a trivial binary efficiency to elements for UL-type
   results as will be explained in detail later.
   
.. [#f2] In order to shorten the notation we sometimes refer to  Z\ :sub:`2`-even final states
   simply as ''final states''. This should not be confused with the Z\ :sub:`2`-odd :ref:`final state
   particles <final stateOdd>`.
