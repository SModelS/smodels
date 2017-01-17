.. index:: What's New in version |version|

.. |invisible compression| replace:: :ref:`invisible compression <invComp>`
.. |mass compression| replace:: :ref:`mass compression <massComp>`
.. |element| replace:: :ref:`element <element>`
.. |elements| replace:: :ref:`elements <element>`
.. |topology| replace:: :ref:`topology <topology>`
.. |topologies| replace:: :ref:`topologies <topology>`
.. |decomposition| replace:: :doc:`decomposition <Decomposition>`
.. |constraint| replace:: :ref:`constraint <ULconstraint>`
.. |constraints| replace:: :ref:`constraints <ULconstraint>`
.. |runSModelS| replace:: :ref:`runSModelS.py <runSModelS>`
.. |database| replace:: :ref:`database <Database>`
.. |Fastlim| replace:: :ref:`Fastlim <addingFastlim>`
.. |output| replace:: :ref:`output <smodelsOutput>`
.. |results| replace:: :ref:`experimental results <ExpResult>`
.. |txnames| replace:: :ref:`txnames <TxName>`
.. |EM| replace:: :ref:`EM-type <EMtype>`
.. |UL| replace:: :ref:`UL-type <ULtype>`
.. |EMr| replace:: :ref:`EM-type result <EMtype>`
.. |ULr| replace:: :ref:`UL-type result <ULtype>`
.. |EMrs| replace:: :ref:`EM-type results <EMtype>`
.. |ULrs| replace:: :ref:`UL-type results <ULtype>`
.. |ExpRes| replace:: :ref:`Experimental Result<ExpResult>`
.. |ExpRess| replace:: :ref:`Experimental Results<ExpResult>`
.. |expres| replace:: :ref:`experimental result<ExpResult>`
.. |express| replace:: :ref:`experimental results<ExpResult>`
.. |Dataset| replace:: :ref:`Data Set<DataSet>`
.. |Datasets| replace:: :ref:`Data Sets<DataSet>`
.. |dataset| replace:: :ref:`data set<DataSet>`
.. |datasets| replace:: :ref:`data sets<DataSet>`
.. |parameters| replace:: :ref:`parameters file <parameterFile>`
.. |ssigBRe| replace:: :math:`\sum \sigma \times BR \times \epsilon`



What's New in Version 1.1
=========================
Since the publication of SModelS v1.0 in December 2014, the code base
has undergone significant structural changes. Version 1.1 comes
with many new features. The major novelties of this release are 
as follows:

* the inclusion of efficiency maps (see |EMrs|)
* a new and more flexible database format (see :ref:`Database structure <databaseStruct>`)
* inclusion of likelihood and :math:`\chi^2` calculation for |EMrs| 
  (see :ref:`likelihood calculation <likelihoodCalc>`)
* extended information on :ref:`topology coverage <topCoverage>` 
* inclusion of a database broswer tool for easy access to the information
  stored in the database (see :ref:`database browser <databaseBrowser>`)
* the database now supports also a more efficient :ref:`binary format <databasePickle>`
* performance improvement for the |decomposition| of the input model
* inclusion of new simplified results to the |database| (including a few 13 TeV results) 
* |Fastlim| efficiency maps can now also be used in SModelS


.. |decomposition| replace:: :doc:`decomposition <Decomposition>`
