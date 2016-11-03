.. index:: Confronting Theory Predictions with Data

.. |EM| replace:: :ref:`EM-type <EMtype>`
.. |UL| replace:: :ref:`UL-type <ULtype>`
.. |EMr| replace:: :ref:`EM-type result <EMtype>`
.. |ULr| replace:: :ref:`UL-type result <ULtype>`
.. |EMrs| replace:: :ref:`EM-type results <EMtype>`
.. |ULrs| replace:: :ref:`UL-type results <ULtype>`
.. |ExpRes| replace:: :ref:`Experimental Result<ExpResult>`
.. |ExpRess| replace:: :ref:`Experimental Results<ExpResult>`
.. |Dataset| replace:: :ref:`Data Set<DataSet>`
.. |Datasets| replace:: :ref:`Data Sets<DataSet>`
.. |dataset| replace:: :ref:`data set<DataSet>`
.. |datasets| replace:: :ref:`data sets<DataSet>`
.. |element| replace:: :ref:`element <element>`
.. |elements| replace:: :ref:`elements <element>`
.. |topology| replace:: :ref:`topology <topology>`
.. |topologies| replace:: :ref:`topologies <topology>`
.. |sigBR| replace:: :math:`\sigma \times BR`
.. |sigBRe| replace:: :math:`\sigma \times BR \times \epsilon`
.. |ssigBRe| replace:: :math:`\sum \sigma \times BR \times \epsilon`

.. _confrontPredictions:


Confronting Predictions with Experimental Limits
================================================

Once the relevant signal cross-sections (or :ref:`theory predictions <theoryPredictions>`) have been computed
for the input model, these must be compared to the respective upper limits.
The upper limits for the signal are stored in the SModelS :ref:`Database <database>`
and depend on the type of |ExpRes|: |UL| or |EM|.

In the case of a |ULr|, the theory predictions typically consist of a list of signal
cross-sections (one for each cluster) for
the single |Dataset| (see :ref:`Theory  Predictions for Upper Limit Results <thePredUL>` for more details).
Each theory prediction must then be compared to its
corresponding upper limit.  This limit is simply the cross-section upper limit provided by
the experimental publication or conference note and is extracted from the corresponding UL map (see |ULrs|).

For |EMrs| there is a single cluster for each |Dataset|, and hence a single signal cross-section
value. This value must be compared to the upper limit for the corresponding signal region.
This upper limit is easily computed using the number of observed and expected events for the |Dataset|
and their uncertainties and is typically stored in the :ref:`Database <database>`.
Since most |EMrs| have several signal regions (|Datasets|), there will be one theory prediction/upper limit
for each |Dataset|. By default SModelS keeps only the best |Dataset|, which is the one which maximizes
the ratio :math:`\mbox{(expected signal)}/\mbox{(expected background)}`.
Thus each |EMr| will have a single theory prediction/upper limit, corresponding  to the best |Dataset|.
If the user wants to have access to all the |datasets|, the default
behavior can be disabled using the variable *useBestDataset*.


The procedure described above can be applied to all the |ExpRess| in the database, resulting
in a list of theory predictions and upper limits for each |ExpRes|. A model can then be considered
excluded by the experimental results if, for one or more predictions, we have *theory prediction* :math:`>` *upper limit* [*]_.

* **The upper limits for a given**  |ULr| **or** |EMr| **can be obtained using the** `getUpperLimitFor  method <../../../documentation/build/html/experiment.html#experiment.expResultObj.ExpResult.getUpperLimitFor>`_

For |EMrs| a likelihood can also be computed (in addition to the upper limits above). This is a simple Poisson
convoluted with a Gaussian using event counts and errors on the number of background events (from the :ref:`Database <database>`) and
signal efficiencies (by default 20%).

A :math:`\chi^2` is also computed from the likelihood using a test statistic 

-2*log(likelihood(:math:`s_{\mathrm{theo}}`)/max(likelihood (:math:`s_{\mathrm{max}}`))

where the likelihood in the denominator is maximized over the number of signal events `s`.
The maximum likelihood occurs by definition at :math:`n_s = n_{\mathrm{obs}} - n_b`, or the difference between the number
of observed and predicted background events.
The prediction for the number of signal events :math:`s_{\mathrm{theo}}` is computed from the given theory input.




.. [*] The statistical significance of the exclusion statement is difficult to quantify exactly, since the model
   is being tested by a large number of results simultaneously.

