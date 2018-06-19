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
behavior can be disabled by setting *useBestDataset=False* in `theoryPredictionsFor <theory.html#theory.theoryPrediction.theoryPredictionsFor>`_ (see :ref:`Example.py <exampleCode>`).


The procedure described above can be applied to all the |ExpRess| in the database, resulting
in a list of theory predictions and upper limits for each |ExpRes|. A model can then be considered
excluded by the experimental results if, for one or more predictions, we have *theory prediction* :math:`>` *upper limit* [*]_.

* **The upper limits for a given**  |ULr| **or** |EMr| **can be obtained using the** `getUpperLimitFor  method <experiment.html#experiment.expResultObj.ExpResult.getUpperLimitFor>`_

.. _likelihoodCalc:

Likelihood Computation
----------------------


In the case of |EMrs|, additional statistical information
about the constrained model can be provided by the SModelS output.
Following the procedure detailed in `CMS NOTE 2017-001 <https://cds.cern.ch/record/2242860?ln=en>`_, we construct a simplified
likelihood which describes the plausibility of the data :math:`D`, given a signal strength :math:`\mu`:

.. math::
   \mathcal{L}(\mu,\theta|D) =  P\left(D|\mu + b + \theta \right) p(\theta)


Here, :math:`\theta` denotes the nuisance parameter that describes the
variations in the signal and background contribtions due to systematic
effects. We assume :math:`p(\theta)` to follow a Gaussian distribution centered
around zero and with a variance of :math:`\delta^2`,
whereas :math:`P(D)` corresponds to a counting variable and is thus
properly described by a Poissonian. The complete likelihood thus reads:

.. math::
   \mathcal{L}(\mu,\theta|D) = \frac{(\mu + b + \theta)^{n_{obs}} e^{\mu + b + \theta}}{n_{obs}!} exp \left( -\frac{\theta^2}{2\delta^2} \right)

where :math:`n_{obs}` is the number of observed events in the signal region.
A test statistic :math:`T` can now be constructed from a likelihood ratio test:

.. math::
   \begin{split}T = -2 \ln \frac{H_0}{H_1} = -2 \ln \left(\frac{\mathcal{L}(\mu=n_{\mathrm{signal}},\theta|D)}{sup\{\mathcal{L}(\mu,\theta|D) : \mu \in \mathbb{R}^+ \}}\right)\end{split}

As the signal hypothesis in the numerator presents a special case of the
likelihood in the denominator, the Neyman-Pearson lemma holds, and we
can assume :math:`T` to be distributed according to a :math:`\chi^2` distribution
with one degree of freedom. Because :math:`H_0` assumes the signal strength of
a particular model, :math:`T=0`  corresponds to a perfect match between that
model's prediction and the measured data. :math:`T \gtrsim 1.96` corresponds to
a 95\% confidence level upper limit.
While :math:`n_{\mathrm{obs}}`, :math:`b`  and :math:`\delta_{b}` are directly extracted from
the data set
(coined *observedN*, *expectedBG* and *bgError*, respectively),
:math:`n_{\mathrm{signal}}` is obtained from the calculation of the
theory predictions. A default 20\% systematical uncertainty is assumed for :math:`n_{\mathrm{signal}}`,
resulting in :math:`\delta^2 = \delta_{b}^2 + \left(0.2 n_{\mathrm{signal}}\right)^2`.

SModelS reports the :math:`\chi^2` (:math:`T` values) and likelihood *for each* |EMr|,
together with the observed and expected :math:`r` values.
We note that in the general case analyses may be correlated, so summing up the :math:`T`
values will no longer follow a :math:`\chi^2_{(n)}`  distribution.
Therefore, for a conservative interpretation, only the result with the best expected limit should be used.
Moreover, for a statistically rigorous usage in scans, it is recommended to check that the analysis giving the
best expected limit does not wildly jump within
continuous regions of parameter space that give roughly the same phenomenology.



* **The** :math:`\chi^2` **for a given** |EMr| **is computed using the** `chi2  method <tools.html#tools.simplifiedLikelihoods.LikelihoodComputer.chi2>`_
* **The likelihood for a given** |EMr| **is computed using the** `likelihood  method <tools.html#tools.simplifiedLikelihoods.LikelihoodComputer.likelihood>`_


Combination of Signal Regions
-----------------------------

In case, the experiment provides a covariance matrix, signal regions can be combined.
Just as before, we follow the procedure given in `CMS NOTE 2017-001 <https://cds.cern.ch/record/2242860?ln=en>`_. SModelS allows for a marginalization as well as a profiling of
the nuisances; profiling is the default. As CPU performance is a concern in SModelS, we take
the liberty of aggregating the official results to an acceptable number of aggregate regions, where *acceptable* is typically a number below 20.
In *runSModelS.py*, combining signal regions is turned off per default, only the result from the best expected signal region is reported.
from the best expected signal region is quoted. It can be turned of with the parameter **options:combineSRs**, see :ref:`parameter file <parameterFile>`.

+----------------------------------------+---------------------------------------+------------------------------------------------+
| .. image:: images/T2bbffff_bestSR.png  | .. image:: images/T2bbffff_17.png     | .. image:: images/T2bbffff_44.png              |
|            :width: 300px               |            :width: 300px              |            :width: 300px                       |
| Best signal region                     | 17 aggregate regions                  | 44 signal regions                              |
+----------------------------------------+---------------------------------------+------------------------------------------------+
  
Figure: Comparison of validation plots, for CMS-PAS-SUS-16-052: for best signal region (left), for combination of 17 aggregate signal regions (center), for combination of all 44 signal regions (right).

.. [*] The statistical significance of the exclusion statement is difficult to quantify exactly, since the model
   is being tested by a large number of results simultaneously.
