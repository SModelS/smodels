.. index:: Confronting Theory Predictions with Data

.. |EM| replace:: :ref:`EM-type <EMtype>`
.. |UL| replace:: :ref:`UL-type <ULtype>`
.. |EMr| replace:: :ref:`EM-type result <EMtype>`
.. |ULr| replace:: :ref:`UL-type result <ULtype>`
.. |EMrs| replace:: :ref:`EM-type results <EMtype>`
.. |ULrs| replace:: :ref:`UL-type results <ULtype>`
.. |ExpRes| replace:: :ref:`Experimental Result<ExpResult>`
.. |ExpRess| replace:: :ref:`Experimental Results<ExpResult>`
.. |Dataset| replace:: :ref:`DataSet<DataSet>`
.. |Datasets| replace:: :ref:`DataSets<DataSet>`
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

Once the relevant signal cross sections (or :ref:`theory predictions <theoryPredictions>`) have been computed
for the input model, these must be compared to the respective upper limits.
The upper limits for the signal are stored in the SModelS :ref:`Database <database>`
and depend on the type of |ExpRes|: |UL| or |EM|.

In the case of a |ULr|, the theory predictions typically consist of a list of signal
cross sections (one for each cluster) for
the single |dataset| (see :ref:`Theory  Predictions for Upper Limit Results <thePredUL>` for more details).
Each theory prediction must then be compared to its
corresponding upper limit.  This limit is simply the cross section upper limit provided by
the experimental publication or conference note and is extracted from the corresponding UL map (see |ULrs|).

For |EMrs| there is a single cluster for each |dataset| (or signal region), and hence a single signal cross section
value. This value must be compared to the upper limit for the corresponding signal region.
This upper limit is easily computed using the number of observed and expected events for the |dataset|
and their uncertainties and is typically stored in the :ref:`Database <database>`.
Since most |EMrs| have several signal regions (|datasets|), there will be one theory prediction/upper limit
for each |dataset|. By default SModelS keeps only the best |dataset|, i.e. the one with the largest
ratio :math:`\mbox{(theory prediction)}/\mbox{(expected limit)}`.
Thus each |EMr| will have a single theory prediction/upper limit, corresponding to the best |dataset|
(based on the expected limit).
If the user wants to have access to all the |datasets|, the default
behavior can be disabled by setting *useBestDataset=False* in `theoryPredictionsFor <../../../documentation/build/html/theory.html#theory.theoryPrediction.theoryPredictionsFor>`_ (see :ref:`Example.py <exampleCode>`).


The procedure described above can be applied to all the |ExpRess| in the database, resulting
in a list of theory predictions and upper limits for each |ExpRes|. A model can then be considered
excluded by the experimental results if, for one or more predictions, we have *theory prediction* :math:`>` *upper limit* [*]_.

* **The upper limits for a given**  |ULr| **or** |EMr| **can be obtained using the** `getUpperLimitFor  method <../../../documentation/build/html/experiment.html#experiment.expResultObj.ExpResult.getUpperLimitFor>`_

.. _likelihoodCalc:

Likelihood Computation
----------------------


For |EMrs| a :math:`\chi^2`-value can also be computed (in addition to the upper limits above).
The :math:`\chi^2` is computed from the likelihood using:

.. math::
   \chi^2 = -2 \log\left(\frac{L(n_{\mathrm{signal}})}{L_{\mathrm{max}}}\right)
   
   
 
where :math:`L(n_{\mathrm{signal}})` is the likelihood for a given number of signal events. 
The likelihood is computed using a simple Poisson convoluted with a Gaussian (for the
background and signal uncertainties) as a function of 
the number of observed events (:math:`n_{\mathrm{obs}}`), the number of expected background events
(:math:`n_{b}`) and its error (:math:`\delta_{b}`)
and the number of signal events (:math:`n_{\mathrm{signal}}`) and its error (:math:`\delta_{s})`).
While :math:`n_{\mathrm{obs}}`, :math:`n_{b}` and :math:`\delta_{b}` are directly extracted from the |dataset|,
:math:`n_{\mathrm{signal}}` is obtained from the :ref:`theoryPredictions` calculation and 
:math:`\delta_{s} = 20\%~\cdot n_{\mathrm{signal}}` by default.

* **The** :math:`\chi^2` **for a given** |EMr| **is computed using the** `chi2  method <../../../documentation/build/html/tools.html#tools.statistics.chi2>`_


.. [*] The statistical significance of the exclusion statement is difficult to quantify exactly, since the model
   is being tested by a large number of results simultaneously.

