.. index:: What's New in Version 1.1

What's New in Version 1.1
=========================
Since the publication of SModelS v1.0 in December 2014, the code base
has undergone significant structural changes. Version 1.1 comes
with many new features. The major novelties of this release are 
as follows:

Efficiency Maps
---------------
While v1.0 was written to deal only with "upper limit"-type results,
this release now fully supports efficiency map-type results.
Blah blah refer also to TheoryPrediction.

Database
--------
The database has majorly increased in size.
Whereas v1.0 incorporated 36 analyses and 140 individual results, 
the database that ships with v1.1 features 





SModelS is a Python library that requires Python version 2.6 or later
(but not version 3).  Internally, SModelS uses the following tools:

 * `Pythia 6.4.27 <http://arxiv.org/abs/hep-ph/0603175>`_
 * `NLL-fast <http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast>`_ 1.2 (7 TeV), 2.1 (8 TeV), and 3.1 (13 TeV)


Likelihoods
-----------
Thus ::

  python setup.py install

should install the entire project, compile the internal Pythia and NLL-fast versions
using gfortran. It should also resolve the external dependencies, i.e. install

Speedups
--------

Yahooo
