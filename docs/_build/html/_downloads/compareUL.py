
# coding: utf-8

## How To: Compare theory predictions with experimental limits

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys
sys.path.append("../")


# In[2]:

#Import those parts of smodels that are needed for this exercise
#(We will assume the input is a SLHA file. For LHE files, use the lheDecomposer instead)
from smodels.theory import slhaDecomposer
from smodels.installation import installDirectory
from smodels.tools.physicsUnits import fb, GeV
from smodels.theory.theoryPrediction import theoryPredictionFor
from smodels.experiment import smsAnalysisFactory, smsHelpers


# In[3]:

#Define the SLHA input file name
filename="%s/inputFiles/slha/gluino_squarks.slha" % installDirectory()


# In[4]:

#Load the database, do the decomposition and compute theory predictions:
#(Look at the theory predictions HowTo to learn how to compute theory predictions)
smsHelpers.base="../test/database/"
listofanalyses = smsAnalysisFactory.load()
listOfTopologies = slhaDecomposer.decompose (filename, sigcut = 0.03 * fb, doCompress=True, doInvisible=True,minmassgap = 5* GeV)
analysesPredictions = [theoryPredictionFor(analysis, listOfTopologies) for analysis in listofanalyses]


# In[5]:

#Print the value of each theory prediction (cluster) for each analysis and the corresponding analysis upper limit:
for anaPrediction in analysesPredictions:
    if not anaPrediction: continue #skip analyses without results
    for theoryPred in anaPrediction:
        print "Analysis name = ",theoryPred.analysis.label
        print "Theory prediction = ",theoryPred.value[0].value
        print "Upper limit = ",theoryPred.analysis.getUpperLimitFor(theoryPred.mass)

