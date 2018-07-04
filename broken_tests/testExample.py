#!/usr/bin/env python

"""
.. module:: testExample
   :synopsis: Tests a simple full smodels run

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
from smodels.share.models.MSSMparticles import BSMList
from smodels.share.models.SMparticles import SMList
import unittest
from smodels.tools.physicsUnits import GeV, fb, pb
from smodels.theory import decomposer
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObj import Database
from smodels.theory.model import Model
import pickle

class ExampleTest(unittest.TestCase):        
    def testExample(self):
        """
        Main program. Displays basic use case.
    
        """
        slhafile = '../inputFiles/slha/lightEWinos.slha'
        model = Model(BSMList,SMList,slhafile)
        model.updateParticles()
    
        #Set main options for decomposition:
        sigmacut = 0.3 * fb
        mingap = 5. * GeV    
        """ Decompose model (use slhaDecomposer for SLHA input or lheDecomposer for LHE input) """
        smstoplist = decomposer.decompose(model, sigmacut, doCompress=True, 
                                          doInvisible=True, minmassgap=mingap)
        
        
        topweights = [0.45069449699999997, 0.8060860593513246, 303.2711905718159, 
                      0.3637955190752928, 15.339097274505018, 3358.5703119038644, 
                      2.160128693731249, 4.5031659201699235, 9.486866215694839, 
                      17.980093334558738, 2.164988614601289, 147.823524958894, 
                      0.8102285793970079, 39.61814354617472, 16.436088764209956, 
                      9.291691094022376, 299.60965243794095, 46.34557917490008, 
                      15.658947303819541, 7.991672275037888, 158.71034751774343, 
                      56.78304329128415, 15.597472309381796]
       
        self.assertEqual(len(smstoplist), 23)        
        #self.assertEqual(len(smstoplist.getElements()), 669)
        for itop,top in enumerate(smstoplist):
            self.assertAlmostEqual(top.getTotalWeight()[0].value.asNumber(fb), 
                                   topweights[itop],4)

        
        # Load all analyses from database
        database = Database("./database",force_load = "txt")
        listOfExpRes = database.getExpResults()
        self.assertEqual(len(listOfExpRes), 4)
                
        # Compute the theory predictions for each analysis
        for expResult in listOfExpRes[1:2]:
            predictions = theoryPredictionsFor(expResult, smstoplist, useBestDataset=False)            
            if not predictions: continue
            print('\n',expResult.getValuesFor('id')[0])
            for theoryPrediction in predictions:
                dataset = theoryPrediction.dataset
                datasetID = dataset.getValuesFor('dataId')[0]            
                mass = theoryPrediction.mass
                txnames = [str(txname) for txname in theoryPrediction.txnames]
#                 PIDs =  theoryPrediction.PIDs         
                print("------------------------")
                print("Dataset = ",datasetID)   #Analysis name
                print("TxNames = ",txnames)   
                print("Prediction Mass = ",mass)    #Value for average cluster mass (average mass of the elements in cluster)
#                 print("Prediction PIDs = ",PIDs)    #Value for average cluster mass (average mass of the elements in cluster)
                print("Theory Prediction = ",theoryPrediction.value[0].value)   #Value for the cluster signal cross-section
                print("Condition Violation = ",theoryPrediction.conditions)  #Condition violation values
                  
                #Get upper limit for the respective prediction:
                if expResult.getValuesFor('dataType')[0] == 'upperLimit':
                    ul = expResult.getUpperLimitFor(txname=theoryPrediction.txnames[0],mass=mass)                     
                elif expResult.getValuesFor('dataType')[0] == 'efficiencyMap':
                    ul = expResult.getUpperLimitFor(dataID=datasetID)
                else: print('weird:',expResult.getValuesFor('dataType'))
                print("Theory Prediction UL = ",ul)
                print("R = ",theoryPrediction.value[0].value/ul)
        

        
if __name__ == "__main__":
    unittest.main()

