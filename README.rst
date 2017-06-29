.. image:: docs/manual/source/images/banner720.png

==============
SModelS v1.1
==============

**SModelS -- A tool for interpreting simplified-model results from the LHC.**

SModelS is an automatic, public tool for interpreting simplified-model results
from the LHC. It is based on a general procedure to decompose Beyond the
Standard Model (BSM) collider signatures presenting a Z\ :sub:`2` symmetry into
Simplified Model Spectrum (SMS) topologies. Our method provides a way to cast
BSM predictions for the LHC in a model independent framework, which can be
directly confronted with the relevant experimental constraints. Our concrete
implementation currently focusses on supersymmetry searches with missing
energy, for which a large variety of SMS results from ATLAS and CMS are
available. 


For instructions on how to install SModelS, see INSTALLATION.rst.


Running SModelS
===============

SModelS provides a command-line tool (runSModelS.py) for the basic functionalities,
which can be executed as:

*./runSModelS.py -p <parameter file> -f <input file or directory> -o <output directory>*

For help instructions:

*./runSModelS.py -h*

An example file on how to call the SModelS libraries from your own
Python code can be found in *Example.py*.

Detailed explanations on how to use SModelS, including explanations of the
output, can be found in the section ''Running SModelS'' of the `SModelS online manual`_.

A few example input files are provided in the inputFiles folder and can be
used to test *runSModelS.py*.


Citation
========

If you use this software please cite the SModelS v1.1 manual_, 
the original SModelS publication_, as well as the programs
it makes use of (pythia8_/pythia6_, NLL-fast_ and pyslha_). 
For your convenience, the relevant
citations are provided in bibtex format in *references.bib*.

For citing the experimental analyses in the database, you can use
*smodels-database/database.bib*.

.. _manual: https://arxiv.org/abs/1701.06586
.. _publication: https://inspirehep.net/record/1269436
.. _pythia6: https://pythia6.hepforge.org/
.. _pythia8: http://home.thep.lu.se/~torbjorn/Pythia.html
.. _pyslha: http://www.insectnation.org/projects/pyslha.html
.. _NLL-fast: http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast 
.. _SModelS online manual: http://smodels.hephy.at/docs/current/manual/build/html/index.html
