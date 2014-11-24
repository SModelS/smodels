
Recipe: print out the theoretical decomposition of an slha file
===============================================================

.. code:: python

    import sys ## ugly hack to make sure smodels is in your path
    sys.path.append("../")
.. code:: python

    from smodels.theory import slhaDecomposer
    from smodels.tools.physicsUnits import fb, GeV
.. code:: python

    list_of_topos = slhaDecomposer.decompose ( "../inputFiles/slha/andrePT4.slha", 
                                              sigcut = 0.001 * fb, doCompress=True, doInvisible=True, 
                                              minmassgap = 5* GeV) ## perform the decomposition

.. parsed-literal::

    10:57:29.118 WARNING  smodels.theory.slhaDecomposer:132 Ignoring t+ decays
    10:57:29.126 WARNING  smodels.theory.slhaDecomposer:132 Ignoring higgs decays
    10:57:29.126 WARNING  smodels.theory.slhaDecomposer:132 Ignoring H0 decays
    10:57:29.127 WARNING  smodels.theory.slhaDecomposer:132 Ignoring A0 decays
    10:57:29.127 WARNING  smodels.theory.slhaDecomposer:132 Ignoring H+ decays


.. code:: python

    list_of_topos.printout() ## printout the whole topology table

.. parsed-literal::

       ======================================================= 
     || 	 						 || 
     || 	 	 Global topologies table 	 	 ||
     || 	 						 || 
       ======================================================= 
    ===================================================== 
    Number of vertices: [2, 3] 
    	 Number of vertex parts: [[1, 0], [1, 1, 0]]
    	 Total Global topology weight:
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.27E+01 [fb]
    ===================================================== 
    Number of vertices: [2, 2] 
    	 Number of vertex parts: [[1, 0], [1, 0]]
    	 Total Global topology weight:
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.52E+02 [fb]
    ===================================================== 
    Number of vertices: [2, 1] 
    	 Number of vertex parts: [[1, 0], [0]]
    	 Total Global topology weight:
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.11E+02 [fb]
    ===================================================== 
    Number of vertices: [3, 3] 
    	 Number of vertex parts: [[1, 1, 0], [1, 1, 0]]
    	 Total Global topology weight:
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.48E+01 [fb]
    ===================================================== 
    Number of vertices: [1, 1] 
    	 Number of vertex parts: [[0], [0]]
    	 Total Global topology weight:
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.37E+01 [fb]
    ===================================================== 
    Number of vertices: [1, 3] 
    	 Number of vertex parts: [[0], [1, 1, 0]]
    	 Total Global topology weight:
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.55E+01 [fb]
    


.. code:: python

    for topo in list_of_topos:  ## print out all elements
        for element in topo.elementList:
            element.printout()                   

.. parsed-literal::

    	 Particles in topology: [[['W-']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.04E-03 [fb]
    
    	 Particles in topology: [[['W-']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.03E-03 [fb]
    
    	 Particles in topology: [[['W-']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.03E-03 [fb]
    
    	 Particles in topology: [[['W-']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.02E-03 [fb]
    
    	 Particles in topology: [[['W-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.76E-03 [fb]
    
    	 Particles in topology: [[['W-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.76E-03 [fb]
    
    	 Particles in topology: [[['W+']], [['ta-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.04E-03 [fb]
    
    	 Particles in topology: [[['nu']], [['ta-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.27E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.27E-01 [fb]
    
    	 Particles in topology: [[['W+']], [['e-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.03E-03 [fb]
    
    	 Particles in topology: [[['nu']], [['e-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['W+']], [['mu-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.03E-03 [fb]
    
    	 Particles in topology: [[['nu']], [['mu-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['W+']], [['nu'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.02E-03 [fb]
    
    	 Particles in topology: [[['nu']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[['W+']], [['nu'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.76E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['nu'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.27E-01 [fb]
    
    	 Particles in topology: [[['e+']], [['nu'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['mu+']], [['nu'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['nu'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[['W+']], [['nu'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.76E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['nu'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.27E-01 [fb]
    
    	 Particles in topology: [[['e+']], [['nu'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['mu+']], [['nu'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['nu'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['W+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.92E-02 [fb]
    
    	 Particles in topology: [[['nu']], [['Z'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.89E-02 [fb]
    
    	 Particles in topology: [[['nu']], [['higgs'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.70E-02 [fb]
    
    	 Particles in topology: [[['ta-']], [['W+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.61E-02 [fb]
    
    	 Particles in topology: [[['ta-']], [['Z'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.26E-02 [fb]
    
    	 Particles in topology: [[['ta-']], [['higgs'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.13E-02 [fb]
    
    	 Particles in topology: [[['nu']], [['W-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.31E-02 [fb]
    
    	 Particles in topology: [[['nu']], [['Z'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.30E-03 [fb]
    
    	 Particles in topology: [[['nu']], [['higgs'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.66E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.08E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.06E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.06E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.04E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.65E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.65E-03 [fb]
    
    	 Particles in topology: [[['W+']], [['ta+'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.78E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.34E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.60E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.60E-01 [fb]
    
    	 Particles in topology: [[['W+']], [['ta-'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.78E-03 [fb]
    
    	 Particles in topology: [[['ta-']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.34E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.60E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.60E-01 [fb]
    
    	 Particles in topology: [[['W+']], [['nu'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.01E-02 [fb]
    
    	 Particles in topology: [[['nu']], [['nu'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.97E+00 [fb]
    
    	 Particles in topology: [[['W+']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.54E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.62E-01 [fb]
    
    	 Particles in topology: [[['e+']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['mu+']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[['W+']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.54E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.62E-01 [fb]
    
    	 Particles in topology: [[['e+']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['mu+']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[['W+']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.54E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.62E-01 [fb]
    
    	 Particles in topology: [[['e+']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['mu+']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[['W+']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.54E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.62E-01 [fb]
    
    	 Particles in topology: [[['e+']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['mu+']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['ta-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.73E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['ta-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.73E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.50E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.50E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.50E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.50E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['e-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['e-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['mu-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['mu-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['nu'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['nu'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['nu'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['nu'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['W-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.31E-02 [fb]
    
    	 Particles in topology: [[['ta+']], [['Z'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.30E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['higgs'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.66E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['W-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.62E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['Z'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.68E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['higgs'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.30E-03 [fb]
    
    	 Particles in topology: [[['ta-']], [['W+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.62E-03 [fb]
    
    	 Particles in topology: [[['Z']], [['W+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.00E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['W+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.70E-03 [fb]
    
    	 Particles in topology: [[['ta-']], [['Z'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.68E-03 [fb]
    
    	 Particles in topology: [[['Z']], [['W-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.00E-03 [fb]
    
    	 Particles in topology: [[['ta-']], [['higgs'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.30E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['W-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.70E-03 [fb]
    
    	 Particles in topology: [[['W-']], [['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.04E-03 [fb]
    
    	 Particles in topology: [[['W-']], [['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.03E-03 [fb]
    
    	 Particles in topology: [[['W-']], [['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.03E-03 [fb]
    
    	 Particles in topology: [[['W-']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.02E-03 [fb]
    
    	 Particles in topology: [[['ta-']], [['W+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.04E-03 [fb]
    
    	 Particles in topology: [[['ta-']], [['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.12E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.09E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.09E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['W+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.03E-03 [fb]
    
    	 Particles in topology: [[['e-']], [['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.09E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.05E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.05E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['W+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.03E-03 [fb]
    
    	 Particles in topology: [[['mu-']], [['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.09E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.05E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.05E-01 [fb]
    
    	 Particles in topology: [[['nu']], [['W+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.02E-03 [fb]
    
    	 Particles in topology: [[['nu']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.00E-01 [fb]
    
    	 Particles in topology: [[['mu+']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.97E+01 [fb]
    
    	 Particles in topology: [[['nu']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.62E+01 [fb]
    
    	 Particles in topology: [[['nu']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.02E+01 [fb]
    
    	 Particles in topology: [[['e+']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.07E+01 [fb]
    
    	 Particles in topology: [[['nu']], [['mu-']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 8.60E+00 [fb]
    
    	 Particles in topology: [[['e+']], [['e-']]]
    	 The element masses are 
    	 Branch 0: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.33E+00 [fb]
    
    	 Particles in topology: [[['mu+']], [['mu-']]]
    	 The element masses are 
    	 Branch 0: [4.20E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.20E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.23E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.80E-02 [fb]
    
    	 Particles in topology: [[['Z']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.52E-02 [fb]
    
    	 Particles in topology: [[['higgs']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.26E-02 [fb]
    
    	 Particles in topology: [[['ta+']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.20E-02 [fb]
    
    	 Particles in topology: [[['W+']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.61E-02 [fb]
    
    	 Particles in topology: [[['mu+']], [['mu-']]]
    	 The element masses are 
    	 Branch 0: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.44E+00 [fb]
    
    	 Particles in topology: [[['nu']], [['e-']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.92E+00 [fb]
    
    	 Particles in topology: [[['nu']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.60E-02 [fb]
    
    	 Particles in topology: [[['higgs']], [['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.08E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.06E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.06E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.04E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['W+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.78E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.34E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.99E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['W+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.78E-03 [fb]
    
    	 Particles in topology: [[['ta-']], [['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.34E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.99E-01 [fb]
    
    	 Particles in topology: [[['e+']], [['e-']]]
    	 The element masses are 
    	 Branch 0: [4.20E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.20E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.23E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.73E-01 [fb]
    
    	 Particles in topology: [[['ta-']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.73E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['e-']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['mu-']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['ta+']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.83E+00 [fb]
    
    	 Particles in topology: [[['ta+']], [['nu']]]
    	 The element masses are 
    	 Branch 0: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.28E+01 [fb]
    
    	 Particles in topology: [[['ta+']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.60E-02 [fb]
    
    	 Particles in topology: [[['ta+']], [['W-']]]
    	 The element masses are 
    	 Branch 0: [1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.31E-02 [fb]
    
    	 Particles in topology: [[['ta+']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.34E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['W-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.62E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['Z']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.68E-03 [fb]
    
    	 Particles in topology: [[['ta+']], [['higgs']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.30E-03 [fb]
    
    	 Particles in topology: [[['W+']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.62E-03 [fb]
    
    	 Particles in topology: [[['W+']], [['W-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.22E-03 [fb]
    
    	 Particles in topology: [[['W+']], [['Z']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.00E-03 [fb]
    
    	 Particles in topology: [[['W+']], [['higgs']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.70E-03 [fb]
    
    	 Particles in topology: [[['Z']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.68E-03 [fb]
    
    	 Particles in topology: [[['Z']], [['W-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.00E-03 [fb]
    
    	 Particles in topology: [[['Z']], [['Z']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.45E-03 [fb]
    
    	 Particles in topology: [[['Z']], [['higgs']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.60E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.30E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['W-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.70E-03 [fb]
    
    	 Particles in topology: [[['higgs']], [['higgs']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.17E-03 [fb]
    
    	 Particles in topology: [[['W-']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.02E-03 [fb]
    
    	 Particles in topology: [[['ta-']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[['e-']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['mu-']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['W+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.02E-03 [fb]
    
    	 Particles in topology: [[['ta+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[['e+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['mu+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['mu+']], []]
    	 The element masses are 
    	 Branch 0: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.97E+01 [fb]
    
    	 Particles in topology: [[['nu']], []]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.30E+01 [fb]
    
    	 Particles in topology: [[['e+']], []]
    	 The element masses are 
    	 Branch 0: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.07E+01 [fb]
    
    	 Particles in topology: [[['mu-']], []]
    	 The element masses are 
    	 Branch 0: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 8.60E+00 [fb]
    
    	 Particles in topology: [[['ta+']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.80E-02 [fb]
    
    	 Particles in topology: [[['W+']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.92E-02 [fb]
    
    	 Particles in topology: [[['Z']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.52E-02 [fb]
    
    	 Particles in topology: [[['higgs']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.26E-02 [fb]
    
    	 Particles in topology: [[['ta+']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.20E-02 [fb]
    
    	 Particles in topology: [[['W+']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.61E-02 [fb]
    
    	 Particles in topology: [[['Z']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.89E-02 [fb]
    
    	 Particles in topology: [[['higgs']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.70E-02 [fb]
    
    	 Particles in topology: [[['ta+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.21E-03 [fb]
    
    	 Particles in topology: [[['e+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.17E-03 [fb]
    
    	 Particles in topology: [[['mu+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.17E-03 [fb]
    
    	 Particles in topology: [[['nu']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.14E-03 [fb]
    
    	 Particles in topology: [[['e-']], []]
    	 The element masses are 
    	 Branch 0: [2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.92E+00 [fb]
    
    	 Particles in topology: [[['ta-']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.60E-02 [fb]
    
    	 Particles in topology: [[['W-']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.31E-02 [fb]
    
    	 Particles in topology: [[['higgs']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.04E-03 [fb]
    
    	 Particles in topology: [[['ta+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.99E-01 [fb]
    
    	 Particles in topology: [[['ta-']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.99E-01 [fb]
    
    	 Particles in topology: [[['W+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.01E-02 [fb]
    
    	 Particles in topology: [[['ta+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.03E+00 [fb]
    
    	 Particles in topology: [[['e+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.01E+00 [fb]
    
    	 Particles in topology: [[['mu+']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.01E+00 [fb]
    
    	 Particles in topology: [[['nu']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.97E+00 [fb]
    
    	 Particles in topology: [[['ta-']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.82E-01 [fb]
    
    	 Particles in topology: [[['e-']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.77E-01 [fb]
    
    	 Particles in topology: [[['mu-']], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.77E-01 [fb]
    
    	 Particles in topology: [[['ta-']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.60E-02 [fb]
    
    	 Particles in topology: [[['W-']], []]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.31E-02 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.12E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.09E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.09E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.27E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.27E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.09E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.05E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.05E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.09E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.05E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.05E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.03E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.00E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.27E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.54E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.54E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.27E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.24E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.54E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.54E-01 [fb]
    
    	 Particles in topology: [[['ta+'], ['ta-']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.34E-01 [fb]
    
    	 Particles in topology: [[['ta+'], ['ta-']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta+'], ['ta-']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta+'], ['ta-']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.28E-01 [fb]
    
    	 Particles in topology: [[['ta+'], ['ta-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.60E-01 [fb]
    
    	 Particles in topology: [[['ta+'], ['ta-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.60E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['ta+']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.34E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['ta+']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['ta+']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.30E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['ta+']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.28E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['ta+']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.60E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['ta+']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.60E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['nu']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.03E+00 [fb]
    
    	 Particles in topology: [[['nu'], ['nu']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.01E+00 [fb]
    
    	 Particles in topology: [[['nu'], ['nu']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.01E+00 [fb]
    
    	 Particles in topology: [[['nu'], ['nu']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.00E+00 [fb]
    
    	 Particles in topology: [[['nu'], ['nu']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.61E+00 [fb]
    
    	 Particles in topology: [[['nu'], ['nu']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.61E+00 [fb]
    
    	 Particles in topology: [[['e+'], ['e-']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.62E-01 [fb]
    
    	 Particles in topology: [[['e+'], ['e-']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['e+'], ['e-']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['e+'], ['e-']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.57E-01 [fb]
    
    	 Particles in topology: [[['e+'], ['e-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.98E-01 [fb]
    
    	 Particles in topology: [[['e+'], ['e-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.98E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['e+']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.62E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['e+']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['e+']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['e+']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.57E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['e+']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.98E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['e+']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.98E-01 [fb]
    
    	 Particles in topology: [[['mu+'], ['mu-']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.62E-01 [fb]
    
    	 Particles in topology: [[['mu+'], ['mu-']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['mu+'], ['mu-']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['mu+'], ['mu-']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.57E-01 [fb]
    
    	 Particles in topology: [[['mu+'], ['mu-']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.98E-01 [fb]
    
    	 Particles in topology: [[['mu+'], ['mu-']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.98E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['mu+']], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.62E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['mu+']], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['mu+']], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.59E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['mu+']], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 4.57E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['mu+']], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.98E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['mu+']], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.98E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['ta+'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.73E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['ta-'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.73E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['nu'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.82E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.50E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.50E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.50E-01 [fb]
    
    	 Particles in topology: [[['ta-'], ['nu']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.50E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['ta+'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['ta-'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['nu'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.77E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['e-'], ['nu']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['ta+'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['ta-'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.72E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['nu'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.77E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['mu-'], ['nu']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['ta+'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.71E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['ta-'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.71E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['nu'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.73E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.48E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.48E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.48E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['ta-']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.48E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['ta+'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['ta-'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['nu'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 8.46E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.29E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.29E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.29E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['e-']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.29E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['ta+'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['ta-'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.49E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['nu'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 8.46E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.29E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.29E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.29E-01 [fb]
    
    	 Particles in topology: [[['nu'], ['mu-']], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.29E-01 [fb]
    
    	 Particles in topology: [[['W+'], ['nu']], [['W-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.22E-03 [fb]
    
    	 Particles in topology: [[['W+'], ['nu']], [['Z'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.00E-03 [fb]
    
    	 Particles in topology: [[['W+'], ['nu']], [['higgs'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.70E-03 [fb]
    
    	 Particles in topology: [[['Z'], ['ta+']], [['W-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.00E-03 [fb]
    
    	 Particles in topology: [[['Z'], ['ta+']], [['Z'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.45E-03 [fb]
    
    	 Particles in topology: [[['Z'], ['ta+']], [['higgs'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.30E-03 [fb]
    
    	 Particles in topology: [[['higgs'], ['ta+']], [['W-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.70E-03 [fb]
    
    	 Particles in topology: [[['higgs'], ['ta+']], [['Z'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.30E-03 [fb]
    
    	 Particles in topology: [[['higgs'], ['ta+']], [['higgs'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.17E-03 [fb]
    
    	 Particles in topology: [[], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.00E-01 [fb]
    
    	 Particles in topology: [[], []]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV]]
    	 Branch 1: [1.87E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.62E+01 [fb]
    
    	 Particles in topology: [[], []]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.30E+01 [fb]
    
    	 Particles in topology: [[], []]
    	 The element masses are 
    	 Branch 0: [1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.14E-03 [fb]
    
    	 Particles in topology: [[], []]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.97E+00 [fb]
    
    	 Particles in topology: [[], []]
    	 The element masses are 
    	 Branch 0: [1.78E+02 [GeV]]
    	 Branch 1: [1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.83E+00 [fb]
    
    	 Particles in topology: [[], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[], [['nu'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[], [['nu'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.22E-01 [fb]
    
    	 Particles in topology: [[], [['Z'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.89E-02 [fb]
    
    	 Particles in topology: [[], [['higgs'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.70E-02 [fb]
    
    	 Particles in topology: [[], [['W+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.61E-02 [fb]
    
    	 Particles in topology: [[], [['ta+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.21E-03 [fb]
    
    	 Particles in topology: [[], [['e+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.17E-03 [fb]
    
    	 Particles in topology: [[], [['mu+'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.17E-03 [fb]
    
    	 Particles in topology: [[], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 7.14E-03 [fb]
    
    	 Particles in topology: [[], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.21E-03 [fb]
    
    	 Particles in topology: [[], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [1.78E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.21E-03 [fb]
    
    	 Particles in topology: [[], [['Z'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.30E-03 [fb]
    
    	 Particles in topology: [[], [['higgs'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [1.87E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 5.66E-03 [fb]
    
    	 Particles in topology: [[], [['nu'], ['ta+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 3.00E+00 [fb]
    
    	 Particles in topology: [[], [['nu'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.61E+00 [fb]
    
    	 Particles in topology: [[], [['nu'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 2.61E+00 [fb]
    
    	 Particles in topology: [[], [['e+'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[], [['e-'], ['e+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[], [['mu+'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[], [['mu-'], ['mu+']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 6.06E-01 [fb]
    
    	 Particles in topology: [[], [['nu'], ['ta-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 1.81E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 9.73E-01 [fb]
    
    	 Particles in topology: [[], [['nu'], ['e-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 8.46E-01 [fb]
    
    	 Particles in topology: [[], [['nu'], ['mu-']]]
    	 The element masses are 
    	 Branch 0: [3.72E+02 [GeV]]
    	 Branch 1: [3.72E+02 [GeV], 2.03E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 8.46E-01 [fb]
    
    	 Particles in topology: [[], [['W-'], ['nu']]]
    	 The element masses are 
    	 Branch 0: [1.78E+02 [GeV]]
    	 Branch 1: [4.30E+02 [GeV], 1.87E+02 [GeV], 1.78E+02 [GeV]]
    	 The element weights are: 
    	 Sqrts: 8.00E+00 [TeV]	 Weight: 1.31E-02 [fb]
    


.. code:: python

    