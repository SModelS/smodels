
Recipe: load the database, filtering out only a few results. Will be obsolete with the database browser
=======================================================================================================

.. code:: python

    import sys ## ugly hack to make sure smodels is in your path
    sys.path.append("../")
.. code:: python

    from smodels.experiment import smsAnalysisFactory, smsHelpers
    from smodels.tools.physicsUnits import GeV
.. code:: python

    smsHelpers.base="../test/database/" ## define where the database resides
.. code:: python

    analyses=["SUS12011"]
    topologies=["T1"]
.. code:: python

    list_of_analyses=smsAnalysisFactory.load(analyses,topologies) ## load only the given analyses and topos
.. code:: python

    masses=[[300*GeV,100*GeV],[300*GeV,100*GeV]]
.. code:: python

    for result in list_of_analyses:
        print("result.label:     \t\t%s" % result.label )
        print("result.sqrts:     \t\t%s" % result.sqrts )
        print("result.conditions:\t\t%s" % result.conditions )
        print("result.constraint:\t\t%s" % result.constraint )
        print("upper limit for mass vector:\t\t%s" % result.getUpperLimitFor(masses))

.. parsed-literal::

    result.label:     		SUS12011:T1
    result.sqrts:     		7.00E+00 [TeV]
    result.conditions:		['None']
    result.constraint:		[[['jet','jet']],[['jet','jet']]]
    upper limit for mass vector:		5.12E+00 [pb]


.. code:: python

    