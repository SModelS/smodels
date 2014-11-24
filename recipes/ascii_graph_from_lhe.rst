
Recipe: create an element from an lhe event, with cross sections, and draw it as an ascii graph
===============================================================================================

.. code:: python

    import sys ## ugly hack to make sure smodels is in your path
    sys.path.append("../")
.. code:: python

    from smodels.theory import lheReader, lheDecomposer, crossSection
    from smodels.installation import installDirectory
    from smodels.tools import asciiGraph
.. code:: python

    filename="%s/inputFiles/lhe/T1_1.lhe" % installDirectory() 
.. code:: python

    reader = lheReader.LheReader ( filename )
.. code:: python

    event=reader.next()
.. code:: python

    xsecs=crossSection.getXsecFromLHEFile ( filename )
.. code:: python

    element=lheDecomposer.elementFromEvent ( event, xsecs )
.. code:: python

    print asciiGraph.asciidraw ( element )

.. parsed-literal::

       q  q 
       \ /  
    ----*----
    ----*----
       / \  
       q  q 
    


.. code:: python

    