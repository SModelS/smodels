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
.. |covariace| replace:: :ref:`combination of signal regions <combineSRs>`

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
ratio :math:`\mbox{(theory prediction)}/\mbox{(expected limit)}`. (See below for |covariace|)
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
Following the procedure detailed in `CMS-NOTE-2017-001 <https://cds.cern.ch/record/2242860?ln=en>`_, we construct a simplified
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
   \mathcal{L}(\mu,\theta|D) = \frac{(\mu + b + \theta)^{n_{obs}} e^{-(\mu + b + \theta)}}{n_{obs}!} exp \left( -\frac{\theta^2}{2\delta^2} \right)

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


.. _combineSRs:

Combination of Signal Regions
-----------------------------

If the experiment provides a covariance matrix together with the efficiency maps, signal regions can be combined. 
This is implemented in SModelS v1.1.3 onwards, following as above the simplified likelihood approach described in `CMS-NOTE-2017-001 <https://cds.cern.ch/record/2242860?ln=en>`_. 

SModelS allows for a marginalization as well as a profiling of the nuisances, with profiling being the default (an example for using marginalisation can be found in :ref:`How To's <Examples>`).
Since CPU performance is a concern in SModelS, we try to aggregate the official results, which can comprise >100 signal regions, to an acceptable number of aggregate regions. Here *acceptable* means as few aggregate regions as possible without loosing in precision or constraining power. 
The CPU time scales roughly linearly with the number of signal regions, so aggregating e.g. from 80 to 20 signal regions means gaining a factor 4 in computing time.

Under the assumptions described in `CMS-NOTE-2017-001 <https://cds.cern.ch/record/2242860?ln=en>`_,
the likelihood for the signal hypothesis when combining signal regions is given by:

.. math::
   \mathcal{L}(\mu,\theta|D) = \prod_{i=1}^{N} \frac{(\mu s_i^r + b_i + \theta_i)^{n_{obs}^i} e^{-(\mu s_i^r + b_i + \theta_i)}}{n_{obs}^i!} exp \left( -\frac{1}{2} \vec{\theta}^T V^{-1} \vec{\theta} \right)

where the product is over all :math:`N` signal regions, :math:`\mu` is the overall signal strength, :math:`s_i^r` the relative signal strength
in each signal region and :math:`V` represents the covariance matrix. 
Note, however, that unlike the case of a single signal region, we do not include any signal uncertainties, since this
should correspond to a second order effect.


Using the above likelihood we compute a 95\% confidence level limit on :math:`\mu` using the :math:`CL_s` (:math:`CL_{sb}/CL_{b}`) limit from the 
test statistic :math:`q_\mu`, as described in Eq. 14 in G. Cowan et al., 
`Asymptotic formulae for likelihood-based tests <https://arxiv.org/abs/1007.1727>`_. 
We then search for the value :math:`CL_s = 0.95` using the Brent bracketing technique available through the `scipy optimize library <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html>`_.
Note that the limit computed through this procedure applies to the total signal yield summed over all signal regions and assumes
that the relative signal strengths in each signal region are fixed by the signal hypothesis. As a result, the above limit has to be computed
for each given input model (or each :ref:`theory prediction <theoryPredictions>`), thus considerably increasing CPU time.

When using *runSModelS.py*, the combination of signal regions is turned on or off with the parameter **options:combineSRs**, see :ref:`parameter file <parameterFile>`. Its default value is *False*, in which case only the result from the best expected signal region (best SR) is reported. 
If *combineSRs = True*, both the combined result and the one from the best SR are quoted. 

In the :ref:`figure below <combinedSRfig>` we show the constraints on the simplified model 
`T2bbWWoff <http://smodels.hephy.at/wiki/SmsDictionary#T2bbWWoff>`_ when using
the best signal region (left), all the 44 signal regions considered in `CMS-PAS-SUS-16-052 <http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/>`_ (center) and the aggregated signal regions included in the SModelS database (right).
As we can see, while the curve obtained from the combination of all 44 signal regions is much closer to the official exclusion than the one obtained using only the best SR. Finally, the aggregated result included in the SModelS database (total of 17 aggregate regions) comes with little loss in constraining power, although it considerable reduces the running time.

.. _combinedSRfig:

+-----------------------------------------+-----------------------------------------+-----------------------------------------+
| .. image:: images/T2bbWWoff_bestSR.png  | .. image:: images/T2bbWWoff_44.png      | .. image:: images/T2bbWWoff_17.png      |
|            :width: 300px                |            :width: 300px                |            :width: 300px                |
| Best signal region                      | 44 signal regions                       | 17 aggregate regions                    |
+-----------------------------------------+-----------------------------------------+-----------------------------------------+


Figure: Comparison of exclusion curves for `CMS-PAS-SUS-16-052 <http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-052/>`_ using only the best signal region (left), the combination of 17 aggregate signal regions (center), and the combination of all 44 signal regions (right).


.. [*] The statistical significance of the exclusion statement is difficult to quantify exactly, since the model
   is being tested by a large number of results simultaneously.
