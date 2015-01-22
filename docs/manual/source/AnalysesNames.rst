.. index:: Analyses Names
.. _Txname:

Analyses Names
==============

Experimental analyses containing cross-section upper limits (see :ref:`UL analysis <ULanalysis>`)
are labelled according to their publication reference. 
For ATLAS analyses, we use arXiv or confnote number, 
for CMS analyses the PAS (public analysis summary) number. 
However, as discussed in :ref:`UL analysis <ULanalysis>`,
a single publication may contain several :ref:`analyses <ULanalysis>`, corresponding to different
:ref:`constraints <ULconstraint>` (or sums of :ref:`elements <element>`).
Therefore the analysis label must also specify to which :ref:`constraint <ULconstraint>` it refers too.
Although the :ref:`constraint <ULconstraint>` can be specified using the :ref:`bracket notation <bracketnotation>`, this is usually too cumbersome for a simple label.
We therefore adopt a notation based on the CMS SMS conventions, where each specific :ref:`constraint <ULconstraint>` is
labeled as *T<constraint name>*.
The complete analysis name is then given by:

<publication label>:T<constraint name>

As an example, consider the analysis with the upper limits for slepton pair production:

.. _constraintplot:

.. image:: images/ATLASslepUL.png
   :height: 460px

The conference note corresponding to the analysis is ATLAS-SUSY-2013-11 and its :ref:`constraint <ULconstraint>` 
(:math:`[[[e^+]],[[e^-]]] + [[[\mu^+]],[[\mu^-]]]`) label is *TSlepSlep*. As a result, SModelS 
uses the analysis name:

*ATLAS-SUSY-2013-11:TSlepSlep*

A complete list of all constraint labels (or Tx names) can be found `here <http://smodels.hephy.at/wiki/SmsDictionary>`_.
