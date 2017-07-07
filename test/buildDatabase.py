#!/usr/bin/python

from __future__ import print_function
import sys, os, time
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.exceptions import SModelSExperimentError
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.colors import colors
from smodels.tools.physicsUnits import pb, fb, GeV
from smodels.theory import slhaDecomposer

colors.on = True
# setLogLevel ( "debug" )

smstoplist = smstoplist = slhaDecomposer.decompose( "T2tt.slha" )
# print ( "smstoplist=",len(smstoplist ) )
dir = "corrdb/"
d=Database( dir, discard_zeroes = True )
print(d)
results=d.getExpResults()

massvec = [[400.*GeV,75*GeV], [400*GeV,75*GeV]]

for e in results[:1]:
    print ( e.globalInfo.id )
    dsets = [ "SR1: MET > 200", "SR2: MET > 300" ]
    topo = "T2tt"
    effs = []
    for ds in dsets:
        eff = e.getEfficiencyFor ( topo, massvec, ds )
        if not eff: continue
        effs.append ( eff )
        ul = e.getUpperLimitFor ( ds )
        eul = e.getUpperLimitFor ( ds, expected=True )
        #print ( "observed upper limit",ds," (on sigma*eff)", ul )
        #print ( "expected upper limit",ds," (on sigma*eff)", eul )
        #print ( "observed upper limit",ds," (on sigma)", ul/eff )
        #print ( "expected upper limit",ds," (on sigma)", eul/eff )
        #print ( "eff", eff )
    try:
        t0=time.time()
        cul = None ## e.getCombinedUpperLimitFor ( effs )
        t1=time.time()
        dt=t1-t0
        print ( "observed upper limit combined (on sigma)", cul,"in",dt,"s" )
        ecul = None ## e.getCombinedUpperLimitFor ( effs, expected=True )
        t2=time.time()
        dt2=t2-t1
        print ( "expected upper limit combined (on sigma)", ecul )
    except SModelSExperimentError as e:
        print ( "exception: %s" % e )
        
    predictions = theoryPredictionsFor ( e, smstoplist, useBestDataset=False, combinedResults=True )
    print ( "theorypredfor=%s" % predictions )
    for theoryPrediction in predictions:
        dataset = theoryPrediction.dataset
        datasetID = dataset.dataInfo.dataId            
        mass = theoryPrediction.mass
        txnames = [str(txname) for txname in theoryPrediction.txnames]
        PIDs =  theoryPrediction.PIDs         
        print ( "------------------------" )
        print ( "Dataset = ",datasetID )  #Analysis name
        print ( "TxNames = ",txnames )  
        print ( "Prediction Mass = ",mass )   #Value for average cluster mass (average mass of the elements in cluster)
        print ( "Prediction PIDs = ",PIDs )   #Value for average cluster mass (average mass of the elements in cluster)
        print ( "Theory Prediction = ",theoryPrediction.xsection )  #Signal cross section
        print ( "Effective efficiency= ",theoryPrediction.effectiveEff ) # effective efficiency
        print ( "Condition Violation = ",theoryPrediction.conditions ) #Condition violation values
          
        ul = theoryPrediction.getUpperLimitFor(txname=txnames[0],mass=mass,dataID=datasetID)                     
        print ( "UL for theory prediction = ",ul )

        # Compute the r-value
        r = theoryPrediction.xsection.value/ul
        print ( "r = ",r )
        #Compute likelihhod and chi^2 for EM-type results:
        if dataset.dataInfo.dataType == 'efficiencyMap':
            theoryPrediction.computeStatistics()
            print ( 'Chi2, likelihood=', theoryPrediction.chi2, theoryPrediction.likelihood )
        if r > 1.0:
            rmax = r
            bestResult = theoryPrediction.expResult.globalInfo.id
        
