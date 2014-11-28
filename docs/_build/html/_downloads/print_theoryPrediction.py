
# coding: utf-8

## How To: Print out the theory predictions

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

# define where the database resides
smsHelpers.base="../test/database/"
#and load the analyses:
listofanalyses = smsAnalysisFactory.load()


# In[4]:

#Define the SLHA input file name
filename="%s/inputFiles/slha/gluino_squarks.slha" % installDirectory()


# In[5]:

#Perform the decomposition:
listOfTopologies = slhaDecomposer.decompose (filename, sigcut = 0.03 * fb, doCompress=True, doInvisible=True,minmassgap = 5* GeV)


# In[6]:

#Compute the theory prediction for each analysis using the results from the decomposition:
analysesPredictions = [theoryPredictionFor(analysis, listOfTopologies) for analysis in listofanalyses]


# In[7]:

#Print information about each theory prediction (cluster) for each analysis:
#(Since this is a test database, there are very few results)
for anaPrediction in analysesPredictions:
    if not anaPrediction: continue #skip analyses without results
    for theoryPred in anaPrediction:
        print "Analysis name = ",theoryPred.analysis.label
        print "Theory prediction = ",theoryPred.value
        print "Conditions violation (if any) = ",theoryPred.conditions
        print "Mass of cluster = ",theoryPred.mass

