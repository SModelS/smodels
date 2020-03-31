.. index:: Theory Definitions

.. |particle| replace:: :ref:`particle <particleClass>`
.. |particles| replace:: :ref:`particles <particleClass>`
.. |element| replace:: :ref:`element <element>`
.. |elements| replace:: :ref:`elements <element>`
.. |topology| replace:: :ref:`topology <topology>`
.. |topologies| replace:: :ref:`topologies <topology>`
.. |bracket notation| replace:: :ref:`bracket notation <bracketNotation>`
.. |parameters| replace:: :ref:`parameters file <parameterFile>`


.. _theoryDefs:

Simplified Model Definitions
============================

The so-called `theory module <theory.html#theory>`_ contains the basic tools necessary for decomposing the input model
into simplified model |topologies| and using the output of the decomposition
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
A representation of an element is shown below:


.. _elementscheme:

.. image:: images/elementB.png
   :width: 40%

An element may also hold information about its corresponding
weight (cross section times branching ratio times efficiency).\ [#f1]_
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
From v2.0 onwards elements hold |particles| and two elements
will be considered equal if both their topological structure
and |particles| are equal.
Below we describe in more detail the element properties and their implementation
in SModelS.


* **Elements are described by the** `Element Class <theory.html#theory.element.Element>`_

.. _particleClass:

Particles
^^^^^^^^^

The basic building block of simplified model |elements| are particles,
which can be both SM (e.g. :math:`l^+,l^-,\nu`  in the :ref:`figure above <elementscheme>`)
or BSM states (e.g. :math:`X1,X2,Y1,Y2,Z1` in the :ref:`figure above <elementscheme>`).
The BSM particles are defined by the input model (see :ref:`model <parameterFileModel>` in |parameters|),
while the SM particles are defined in `SMparticles.py <share.html#share.models.SMparticles>`_ .
All particles must be assigned a Z\ :sub:`2` parity and can have a flexible
number of attributes, such as mass, spin, electric charge, etc.
Two particles are considered equal if all their shared properties
are equal. *Inclusive* or *generic* particles can then be defined if some
of its properties are left undefined. For instance, a particle with electric
charge -1, spin 1/2 but without a defined mass will be matched
to electrons, muons and taus. This is useful when defining generic simplified models
(|elements|) in the :ref:`Database <databaseDefs>`. All *generic* particles
used by the :ref:`Database <databaseDefs>` are separately defined in
`databaseParticles.py <experiment.html#experiment.databaseParticles>`_
Examples for such inclusive definitions are:

 - 'l' for electrons, and muons,
 - 'L' for electrons, muons, and taus,
 - 'q' for u-, d-, and s-quarks,
 - 'jet' for u-, d-, s-, c-quarks and gluons
 - 'anyOdd' for any Z\ :sub:`2`-odd particle
 - '*' for any \ :sub:`2`-even particle


* **Particles are described by the** `Particle Class <theory.html#theory.particle.Particle>`_

.. _vertex:

Vertices
^^^^^^^^

Each Z\ :sub:`2`-odd decay is represented by a vertex containing the outgoing states (one Z\ :sub:`2`-odd
state and the Z\ :sub:`2`-even particles), as shown in the :ref:`scheme above <topscheme>`.

* **Vertices are described by the** `ParticleList Class <theory.html#theory.particle.ParticleList>`_


.. _branch:

Branches
^^^^^^^^

A branch is the basic substructure of an |element|.
It represents a series of cascade decays of a single initial Z\ :sub:`2`-odd
state. The diagram below illustrates an example of a branch.

.. image:: images/branchTopB.png
   :width: 25%

The structure of each branch is fully defined by its number of vertices and the number of
|particles| coming out of each vertex.

* **Branches are described by the** `Branch Class <theory.html#theory.branch.Branch>`_


.. _notation:

Element Representation: Bracket Notation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The structure and Z\ :sub:`2`-even states of |elements| can be represented in a compact and textual form using a nested bracket
notation. The scheme below shows how to convert between the graphical and bracket representations of an element:


.. _bracketnotation:

.. image:: images/bracketNotationB.png
   :width: 50%

The brackets are ordered and nested in the following way.
The outermost brackets correspond to the :ref:`branches <branch>` of the |element|.
The branches are sorted according to their size (see :ref:`element sorting <elementsorting>`)
and each branch contains an *ordered* list of :ref:`vertices <vertex>`.
Each vertex contains a list of the  Z\ :sub:`2`-even particles (represented by their label and sorted alphabetically) coming out of the vertex.
Schematically, for the example in the :ref:`figure above <bracketnotation>`, we have::

   element = [branch1, branch2]
      branch1 = [vertex1]
         vertex1 = [l+,l-]
      branch2 = [vertex1,vertex2]
         vertex1 = [l+]
         vertex2 = [nu]

Although the above scheme can be useful and provides a simplified representation of an |element|,
it provides no information about the Z\ :sub:`2`-odd (BSM) states appearing in the |element|.
However, information about a specific property of Z\ :sub:`2`-odd states can also be represented in a nested bracket notation.
For instance, all the masses of the BSM states in a given |element| can be represented as shown below:

.. _massnotation:

.. image:: images/massNotationB.png
   :width: 65%


Similar arrays can be built with any property (width, charge, spin, etc) of the Z\ :sub:`2`-odd particles in an |element|.


.. _topology:

Topologies
----------

It is often useful to classify |elements| according to their
overall structure or topology.
Each topology corresponds to an *undressed*
|element|, removed of its specific |particle| states.
Therefore the topology is fully determined by its number of
branches, number of vertices in each :ref:`branch <branch>` and number of
|particles| coming out of each :ref:`vertex <vertex>`.
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
