.. index:: Database of Experimental Analyses

.. _database:

.. |constraint| replace:: :ref:`constraint <ULconstraint>`
.. |conditions| replace:: :ref:`conditions <ULconditions>` 
.. |fb-1| replace:: :math:`\mathrm{fb}^{-1}`
.. |sqrts| replace:: :math:`\sqrt{s}`
.. |analyses| replace:: :ref:`analyses <ULanalysis>`

Database of Experimental Analyses
=================================

SModelS stores all the information about the experimental results in an |analyses|
database [*]_. The database is organized as files in an ordinary (UNIX) directory hierarchy,
with a thin python layer serving as the access to the database.


The top level of the SModelS database categorizes the analyses by LHC center-of-mass energies, |sqrts|:

* 7TeV
* 8TeV

Also, the top level directory contains a file called ``version`` with the version string
of the database.

The second level splits the results up between the different experiments:

* 8TeV/CMS/
* 8TeV/ATLAS/

The third level of the directory hierarchy encodes the publications:

* 8TeV/CMS/CMS-PAS-SUS-12-026
* 8TeV/ATLAS/ATLAS-SUSY-2013-12
* ...

For each publication, there are two files:

* ``sms.py`` contains the experimental results for the cross-section upper limits as a function of the mass parameters, as python code. 
  The code will fill a nested ``Dict`` dictionary with the upper limits. The dictionary keys are the corresponding *Txnames* for the
  analysis |constraint| (see :doc:`Analyses Names <AnalysesNames>`).
    
* ``info.txt`` contains all the meta information around the publication. Here is the content of CMS-SUS-13-006 as an example:

.. literalinclude:: /literals/info.txt
   :lines: 1-8

In particular, it defines the center-of-mass energy |sqrts| (in TeV) and the integrated luminosity, |fb-1|.
Optionally, units (e.g. "* TeV" or "/ fb") can be added to the number.

The next block:

.. literalinclude:: /literals/info.txt
   :lines: 9-10,15-16

holds the constraint and conditions, as described in |constraint|, and |conditions|, respectively.

Finally,

.. literalinclude:: /literals/info.txt
   :lines: 21

describes the mass planes and how for analyses with intermediate masses, the
two-dimensional histograms map onto the three-dimensional SMS parameter space.



.. [*] Currently only cross section upper limits are included (see :ref:`UL Analyses <ULanalysis>`). 
