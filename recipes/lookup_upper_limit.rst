
Recipe: look up the upper limit of a particular result, for a particular set of mother and daughter masses
==========================================================================================================

.. code:: python

    import sys ## ugly hack to make sure smodels is in your path
    sys.path.append("../")
.. code:: python

    from smodels.tools.physicsUnits import GeV
    from smodels.experiment.smsInterpolation import upperLimit
.. code:: python

    #specify analysis and topology as strings:
    analysis = "ATLAS_CONF_2013_048"
    topology = "T6bbWW"
.. code:: python

    masses = [500*GeV, 400*GeV, 100*GeV]
.. code:: python

    upperLimit(analysis, topology, masses)



.. parsed-literal::

    1.02E-01 [pb]



.. code:: python

    