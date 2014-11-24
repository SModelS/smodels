
Recipe: compute NLL cross sections for a given slha file, recycling LO cross sections
=====================================================================================

Wolfgang Waltenberger, june 2014
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    import sys ## ugly hack to make sure smodels is in our path
    sys.path.append("../")
.. code:: python

    # now lets import those parts of smodels that are needed for this exercise
    from smodels.tools import xsecComputer
    from smodels.tools.physicsUnits import TeV, fb
    from smodels.tools.xsecComputer import LO, NLL
.. code:: python

    # fetch example1.slha via http, save as tmp.slha
    import urllib
    urllib.urlretrieve("http://smodels.hephy.at/slha/example2.slha","tmp.slha")



.. parsed-literal::

    ('tmp.slha', <httplib.HTTPMessage instance at 0x7fd2877834d0>)



.. code:: python

    # now lets compute the leading order (LO) cross sections for 8 TeV, simulating 1000
    # events with pythia.
    xsecsLO=xsecComputer.computeXSec ( 8*TeV, LO, 1000, "tmp.slha" )
.. code:: python

    print xsecsLO.getDictionary(groupBy="labels")

.. parsed-literal::

    {'8 TeV (LO)': {(1000001, 1000021): 1.22E-02 [pb], (-1000004, -1000001): 3.48E-04 [pb], (-1000024, 1000024): 2.82E-02 [pb], (1000021, 1000021): 1.39E-03 [pb], (-1000024, 1000023): 1.22E-02 [pb], (-2000002, 1000002): 3.48E-04 [pb], (-1000002, 1000023): 6.96E-04 [pb], (1000022, 1000037): 3.48E-04 [pb], (1000003, 1000021): 6.96E-04 [pb], (-1000004, 1000002): 6.96E-04 [pb], (1000001, 1000003): 3.48E-03 [pb], (-1000001, 2000002): 1.04E-03 [pb], (-1000003, 1000003): 2.44E-03 [pb], (-2000003, 1000002): 3.48E-04 [pb], (-1000037, 1000035): 3.48E-04 [pb], (1000035, 1000037): 2.44E-03 [pb], (1000001, 1000002): 5.53E-02 [pb], (1000002, 2000002): 4.18E-03 [pb], (-1000024, 1000002): 7.31E-03 [pb], (1000023, 1000024): 3.34E-02 [pb], (1000001, 2000001): 3.48E-04 [pb], (-1000024, 1000021): 3.48E-04 [pb], (-1000003, 1000021): 6.96E-04 [pb], (-1000002, 1000024): 6.96E-04 [pb], (1000002, 1000023): 6.96E-03 [pb], (-1000004, 1000004): 3.13E-03 [pb], (-1000003, 1000002): 1.39E-03 [pb], (1000002, 1000002): 5.46E-02 [pb], (-1000024, -1000003): 3.48E-04 [pb], (1000025, 1000035): 1.04E-03 [pb], (1000001, 1000001): 5.92E-03 [pb], (-1000001, 2000001): 6.96E-04 [pb], (2000001, 2000002): 3.48E-04 [pb], (1000001, 1000024): 2.44E-02 [pb], (-1000002, 2000002): 6.96E-04 [pb], (-1000002, -1000001): 6.96E-04 [pb], (-1000002, 1000021): 3.48E-04 [pb], (1000001, 1000023): 1.74E-03 [pb], (-1000024, 1000022): 6.96E-04 [pb], (-1000002, 1000002): 3.83E-03 [pb], (1000021, 1000023): 1.39E-03 [pb], (1000021, 2000002): 1.74E-03 [pb], (1000002, 2000001): 1.74E-03 [pb], (1000003, 1000024): 6.96E-04 [pb], (-1000003, 1000001): 3.48E-04 [pb], (1000022, 2000002): 1.04E-03 [pb], (1000002, 1000003): 5.22E-03 [pb], (1000021, 1000037): 3.48E-04 [pb], (-1000001, 1000002): 6.96E-04 [pb], (-1000037, 1000037): 3.48E-04 [pb], (-1000037, 1000002): 3.48E-04 [pb], (-1000024, 1000037): 6.96E-04 [pb], (-1000024, 1000025): 6.96E-04 [pb], (-1000024, 1000004): 3.48E-04 [pb], (1000022, 1000035): 3.48E-04 [pb], (1000023, 1000025): 3.48E-04 [pb], (-1000002, 2000001): 3.48E-04 [pb], (-1000016, 1000016): 3.48E-04 [pb], (1000001, 1000022): 6.96E-04 [pb], (-1000004, -1000002): 3.48E-04 [pb], (-2000001, 1000002): 6.96E-04 [pb], (1000021, 1000022): 1.04E-03 [pb], (1000021, 2000001): 1.39E-03 [pb], (1000002, 1000021): 3.31E-02 [pb], (-1000003, -1000001): 3.48E-04 [pb], (1000023, 1000023): 3.48E-04 [pb], (-1000024, -1000001): 6.96E-04 [pb], (1000002, 1000004): 2.44E-03 [pb], (1000024, 1000025): 3.48E-04 [pb], (-1000003, 1000023): 3.48E-04 [pb], (1000001, 2000002): 2.09E-03 [pb], (1000025, 1000037): 1.74E-03 [pb], (-1000001, 1000021): 1.04E-03 [pb], (-1000004, 1000024): 6.96E-04 [pb], (2000002, 2000002): 1.39E-03 [pb], (1000021, 1000024): 2.09E-03 [pb], (-1000001, 1000001): 3.83E-03 [pb]}}


.. code:: python

    # and write it back into the file
    xsecComputer.addXSecToFile ( xsecsLO, "tmp.slha" )



.. parsed-literal::

    True



.. code:: python

    # now compute the NLL cross sections, recycling the LO computations
    # so only the k-factors are applied, and pythia needs not be rerun
    xsecsNLL=xsecComputer.computeXSec ( 8*TeV, NLL, None, "tmp.slha", loFromSlha=True )

.. parsed-literal::

    15:36:35.195 INFO     smodels.tools.xsecComputer:78  Using LO cross-sections from tmp.slha


.. code:: python

    # ... which has a .getDictionary() method which contains the cross sections
    D=xsecsNLL.getDictionary(groupBy="labels")["8 TeV (NLO+NLL)"]
    print D

.. parsed-literal::

    {(1000001, 1000021): 2.29E+01 [fb], (1000001, 1000003): 3.83E+00 [fb], (1000021, 1000021): 5.64E+00 [fb], (-2000002, 1000002): 5.56E-01 [fb], (1000003, 1000021): 1.31E+00 [fb], (-1000004, 1000002): 1.01E+00 [fb], (-2000001, 1000002): 1.09E+00 [fb], (-1000003, 1000003): 3.53E+00 [fb], (-1000001, 2000002): 1.67E+00 [fb], (1000001, 1000002): 6.09E+01 [fb], (-1000003, 1000002): 2.02E+00 [fb], (2000002, 2000002): 1.78E+00 [fb], (1000002, 1000002): 6.01E+01 [fb], (1000001, 1000001): 6.51E+00 [fb], (-1000001, 2000001): 1.09E+00 [fb], (2000001, 2000002): 4.37E-01 [fb], (-1000002, 2000002): 1.11E+00 [fb], (-1000002, 1000002): 5.54E+00 [fb], (1000021, 2000002): 3.27E+00 [fb], (1000002, 2000001): 2.06E+00 [fb], (-1000003, 1000001): 5.04E-01 [fb], (1000002, 1000003): 5.75E+00 [fb], (-1000001, 1000002): 1.01E+00 [fb], (1000001, 2000002): 2.51E+00 [fb], (-1000002, 2000001): 5.45E-01 [fb], (-1000004, 1000004): 4.54E+00 [fb], (1000021, 2000001): 2.62E+00 [fb], (1000002, 1000021): 6.21E+01 [fb], (1000001, 2000001): 4.12E-01 [fb], (-2000003, 1000002): 5.45E-01 [fb], (1000002, 2000002): 4.98E+00 [fb], (1000002, 1000004): 2.68E+00 [fb], (-1000001, 1000001): 5.54E+00 [fb]}


.. code:: python

    # now lets make a simple bar chart of all cross sections, in fb 
    import pylab; import numpy; pylab.bar( range(len(D)), map ( lambda x: float(x/fb), D.values() ) )
    pylab.xticks( .5+ numpy.arange(len(D)), D.keys(), rotation="vertical" ); pylab.ylabel( "xsec [fb]");


.. image:: nll_xsecs_from_slha_files/nll_xsecs_from_slha_10_0.png


.. code:: python

    