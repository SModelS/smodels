.. image:: https://smodels.github.io/pics/banner.png

.. |PyPI version| image:: https://badge.fury.io/py/smodels.svg
   :target: https://badge.fury.io/py/smodels

.. |Anaconda version| image:: https://anaconda.org/conda-forge/smodels/badges/version.svg
   :target: https://anaconda.org/conda-forge/smodels/

.. |GitHub Project| image:: https://img.shields.io/badge/GitHub--blue?style=social&logo=GitHub
   :target: https://github.com/SModelS

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1169739.svg
   :target: https://doi.org/10.5281/zenodo.116973

.. |CodeFactor| image:: https://www.codefactor.io/repository/github/smodels/smodels/badge/main
   :target: https://www.codefactor.io/repository/github/smodels/smodels/overview/main

.. |Docs| image:: https://img.shields.io/badge/docs-main-blue.svg                    
   :target: https://smodels.readthedocs.io

.. |Colab| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/SModelS/tutorials/blob/main/index.ipynb

|GitHub Project| |PyPI version| |Anaconda version| |CodeFactor| |Colab| |Docs|

==============
SModelS v3
==============

**SModelS -- A tool for interpreting simplified-model results from the LHC.**

SModelS is an automatic, public tool for interpreting simplified-model results
from the LHC. It is based on a general procedure to decompose Beyond the
Standard Model (BSM) collider signatures into
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

If you use this software please cite the SModelS v1-v3 manuals, the original
SModelS publication, as well as the programs it makes use of.  For your
convenience, the relevant citations are provided in bibtex format in
`references.bib <https://github.com/SModelS/smodels/blob/main/references.bib>`_.

For citing the experimental analyses in the database, you can use
`database.bib <https://github.com/SModelS/smodels-database-release/blob/main/database.bib>`_.

.. _SModelS online manual: http://smodels.readthedocs.io/
