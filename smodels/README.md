![Banner](https://smodels.github.io/pics/banner.png)

[![GitHub Project](https://img.shields.io/badge/GitHub--blue?style=social&logo=GitHub)](https://github.com/SModelS)
[![PyPI version](https://badge.fury.io/py/smodels.svg)](https://badge.fury.io/py/smodels)
[![Anaconda version](https://anaconda.org/conda-forge/smodels/badges/version.svg)](https://anaconda.org/conda-forge/smodels/)
[![CodeFactor](https://www.codefactor.io/repository/github/smodels/smodels/badge/main)](https://www.codefactor.io/repository/github/smodels/smodels/overview/main)
[![Docs](https://img.shields.io/badge/docs-main-blue.svg)](https://smodels.readthedocs.io)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1169739.svg)](https://doi.org/10.5281/zenodo.116973)

# SModelS v3

**SModelS — A tool for interpreting simplified-model results from the LHC.**

SModelS is an automatic, public tool for interpreting simplified-model results from the LHC. It is based on a general procedure to decompose Beyond the Standard Model (BSM) collider signatures into Simplified Model Spectrum (SMS) topologies. Our method provides a way to cast BSM predictions for the LHC in a model independent framework, which can be directly confronted with the relevant experimental constraints.

## Installation

For instructions on how to install SModelS, see the [Installation](https://smodels.readthedocs.io/en/latest/Installation.html) section of the 
[SModelS online manual](https://smodels.readthedocs.io/).

## Running SModelS

SModelS provides a command-line tool (`runSModelS.py`) for the basic functionalities, which can be executed as:

```bash
./runSModelS.py -p <parameter file> -f <input file or directory> -o <output directory>
```

For help instructions:

```bash
./runSModelS.py -h
```

An example file on how to call the SModelS libraries from your own Python code
can be found in `Example.py`.

Detailed explanations on how to use SModelS, including explanations of the
output, can be found in the [Using
SModelS](https://smodels.readthedocs.io/en/latest/RunningSModelS.html) section
of the [SModelS online manual](https://smodels.readthedocs.io/).

A few example input files are provided in the `inputFiles` folder and can be
used to test `runSModelS.py`.

## Citation

If you use this software please cite the SModelS v1–v3 manuals, the original
SModelS publication, as well as the programs it makes use of. For your
convenience, the relevant citations are provided in bibtex format in 
[references.bib](https://github.com/SModelS/smodels/blob/main/references.bib).

For citing the experimental analyses in the database, you can use:
[database.bib](https://github.com/SModelS/smodels-database-release/blob/main/database.bib).
