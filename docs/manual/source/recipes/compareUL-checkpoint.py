
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
filename="%s/inputFiles/slha/lightSquarks.slha" % installDirectory()


# In[4]:

#Load the database, do the decomposition and compute theory predictions:
#(Look at the theory predictions HowTo to learn how to compute theory predictions)
smsHelpers.base="../test/database/"
listofanalyses = smsAnalysisFactory.load()
listOfTopologies = slhaDecomposer.decompose (filename, sigcut = 0.03 * fb, doCompress=True, doInvisible=True,minmassgap = 5* GeV)
analysesPredictions = [theoryPredictionFor(analysis, listOfTopologies) for analysis in listofanalyses]


# Out[4]:

#     12:02:03.189 INFO     smodels.experiment.smsAnalysisFactory:61  SUS12022 has been superseded by SUS13006, skipping SUS12022
#     12:02:03.254 INFO     smodels.theory.slhaDecomposer:132 Ignoring t+ decays
#     12:02:03.269 INFO     smodels.theory.slhaDecomposer:132 Ignoring higgs decays
#     12:02:03.270 INFO     smodels.theory.slhaDecomposer:132 Ignoring H0 decays
#     12:02:03.270 INFO     smodels.theory.slhaDecomposer:132 Ignoring A0 decays
#     12:02:03.271 INFO     smodels.theory.slhaDecomposer:132 Ignoring H+ decays
#     12:02:03.335 INFO     smodels.theory.crossSection:512 Ignoring 76 lower order cross-sections
#     12:02:21.672 INFO     smodels.experiment.smsInterpolation:174 Masses out of range for ATLAS_CONF_2013_048/T6bbWW (no extrapolation)
# 

# In[5]:

#Print the value of each theory prediction (cluster) for each analysis and the corresponding analysis upper limit:
for anaPrediction in analysesPredictions:
    if not anaPrediction: continue #skip analyses without results
    for theoryPred in anaPrediction:
        print "Analysis name = ",theoryPred.analysis.label
        print "Theory prediction = ",theoryPred.value[0].value
        print "Upper limit = ",theoryPred.analysis.getUpperLimitFor(theoryPred.mass)


# Out[5]:

#     Analysis name =  SUS13006:TChiWZ
#     Theory prediction =  1.85E-02 [pb]
#     Upper limit =  4.44E-01 [pb]
#     Analysis name =  SUS12028:T2
#     Theory prediction =  1.77E-03 [pb]
#     Upper limit =  1.64E-02 [pb]
#     Analysis name =  SUS12028:T1
#     Theory prediction =  3.92E-04 [pb]
#     Upper limit =  3.14E-02 [pb]
# 

# In[5]:



