
Recipe: compute LO cross sections for a given slha file
=======================================================

Wolfgang Waltenberger, june 2014
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    import sys ## ugly hack to make sure smodels is in your path
    sys.path.append("../")
.. code:: python

    # now lets import those parts of smodels that are needed for this exercise
    from smodels.tools import xsecComputer
    from smodels.tools.physicsUnits import TeV, fb
    from smodels.tools.xsecComputer import LO, NLL
.. code:: python

    # fetch example1.slha via http, save as tmp.slha
    import urllib
    urllib.urlretrieve("http://smodels.hephy.at/slha/example1.slha","tmp.slha")



.. parsed-literal::

    ('tmp.slha', <httplib.HTTPMessage instance at 0x7f502820a3b0>)



.. code:: python

    # now lets compute the leading order (LO) cross sections for 8 TeV, simulating 1000
    # events with pythia.
    xsecs=xsecComputer.computeXSec ( 8*TeV, LO, 1000, "tmp.slha" )
.. code:: python

    # the output is a XSectionList ...
    type(xsecs)



.. parsed-literal::

    smodels.theory.crossSection.XSectionList



.. code:: python

    # ... which has a .getDictionary() method which contains the cross sections
    D=xsecs.getDictionary(groupBy="labels")["8 TeV (LO)"]
    print D[(1000023,1000024)]

.. parsed-literal::

    3.50E-02 [pb]


.. code:: python

    # now lets make a simple bar chart of all cross sections, in fb 
    import pylab; import numpy; pylab.bar( range(len(D)), map ( lambda x: float(x/fb), D.values() ) )
    pylab.xticks( .5+ numpy.arange(len(D)), D.keys(), rotation="vertical" ); pylab.ylabel( "xsec [fb]");


.. image:: lo_xsecs_from_slha_files/lo_xsecs_from_slha_8_0.png


.. code:: python

    import os; os.unlink("tmp.slha")
.. code:: python

    