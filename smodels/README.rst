.. image:: https://smodels.github.io/pics/banner.png

.. |PyPI version| image:: https://badge.fury.io/py/smodels.svg
   :target: https://badge.fury.io/py/smodels

.. |GitHub Project| image:: https://img.shields.io/badge/GitHub--blue?style=social&logo=GitHub
   :target: https://github.com/SModelS

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1169739.svg
   :target: https://doi.org/10.5281/zenodo.116973

.. |CodeFactor| image:: https://www.codefactor.io/repository/github/smodels/smodels/badge/master
   :target: https://www.codefactor.io/repository/github/smodels/smodels/overview/master

.. |Docs| image:: https://img.shields.io/badge/docs-master-blue.svg                    
   :target: https://smodels.readthedocs.io

.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/SModelS/pyhep2020/master?filepath=index.ipynb

|GitHub Project| |PyPI version| |CodeFactor| |Binder| |Docs|

==============
SModelS v2
==============

**SModelS -- A tool for interpreting simplified-model results from the LHC.**

SModelS is an automatic, public tool for interpreting simplified-model results
from the LHC. It is based on a general procedure to decompose Beyond the
Standard Model (BSM) collider signatures presenting a Z\ :sub:`2` symmetry into
Simplified Model Spectrum (SMS) topologies. Our method provides a way to cast
BSM predictions for the LHC in a model independent framework, which can be
directly confronted with the relevant experimental constraints.


Installation
============

For instructions on how to install SModelS, see
the section `Installation <http://smodels.readthedocs.io/en/latest/Installation.html>`_ of the `SModelS online manual`_.


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
output, can be found in the section `Using SModelS <http://smodels.readthedocs.io/en/latest/RunningSModelS.html>`_ of the `SModelS online manual`_.

A few example input files are provided in the inputFiles folder and can be
used to test *runSModelS.py*.


Citation
========

If you use this software please cite both the SModelS v1.1_ and v1.2_ manuals,
the original_ SModelS publication, as well as the programs
it makes use of (pythia8_/pythia6_, NLL-fast_ and pyslha_).
If you use specifically the long-lived particles implementation, please cite also this_ paper.
For the Run2 database, cite the v1.2_ manual and the databaseUpdate_ paper.

For your convenience, the relevant
citations are provided in bibtex format in *references.bib*.

For citing the experimental analyses in the database, you can use
*smodels-database/database.bib*.

.. _v1.2: https://inspirehep.net/record/1705426
.. _v1.1: https://inspirehep.net/record/1510436
.. _original: https://inspirehep.net/record/1269436
.. _this: https://inspirehep.net/record/1687820
.. _databaseUpdate: https://inspirehep.net/record/1658765
.. _pythia6: https://pythia6.hepforge.org/
.. _pythia8: http://home.thep.lu.se/~torbjorn/Pythia.html
.. _pyslha: http://www.insectnation.org/projects/pyslha.html
.. _NLL-fast: http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast
.. _SModelS online manual: http://smodels.readthedocs.io/
