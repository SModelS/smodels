
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


# Out[3]:

#     12:12:57.073 INFO     smodels.experiment.smsAnalysisFactory:61  SUS12022 has been superseded by SUS13006, skipping SUS12022
# 

# In[4]:

#Define the SLHA input file name
filename="%s/inputFiles/slha/lightSquarks.slha" % installDirectory()


# In[5]:

#Perform the decomposition:
listOfTopologies = slhaDecomposer.decompose (filename, sigcut = 0.03 * fb, doCompress=True, doInvisible=True,minmassgap = 5* GeV)


# Out[5]:

#     12:12:57.168 INFO     smodels.theory.slhaDecomposer:132 Ignoring t+ decays
#     12:12:57.182 INFO     smodels.theory.slhaDecomposer:132 Ignoring higgs decays
#     12:12:57.185 INFO     smodels.theory.slhaDecomposer:132 Ignoring H0 decays
#     12:12:57.186 INFO     smodels.theory.slhaDecomposer:132 Ignoring A0 decays
#     12:12:57.187 INFO     smodels.theory.slhaDecomposer:132 Ignoring H+ decays
#     12:12:57.249 INFO     smodels.theory.crossSection:512 Ignoring 76 lower order cross-sections
# 

# In[6]:

#Compute the theory prediction for each analysis using the results from the decomposition:
analysesPredictions = [theoryPredictionFor(analysis, listOfTopologies) for analysis in listofanalyses]


# Out[6]:

#     12:13:15.689 INFO     smodels.experiment.smsInterpolation:174 Masses out of range for ATLAS_CONF_2013_048/T6bbWW (no extrapolation)
# 

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


# Out[7]:

#     Analysis name =  SUS13006:TChiWZ
#     Theory prediction =  ['8.00E+00 [TeV]:1.85E-02 [pb]']
#     Conditions violation (if any) =  {'None': None}
#     Mass of cluster =  [[2.69E+02 [GeV], 1.29E+02 [GeV]], [2.69E+02 [GeV], 1.29E+02 [GeV]]]
#     Analysis name =  SUS12028:T2
#     Theory prediction =  ['8.00E+00 [TeV]:1.77E-03 [pb]']
#     Conditions violation (if any) =  {'None': None}
#     Mass of cluster =  [[9.91E+02 [GeV], 1.29E+02 [GeV]], [9.91E+02 [GeV], 1.29E+02 [GeV]]]
#     Analysis name =  SUS12028:T1
#     Theory prediction =  ['8.00E+00 [TeV]:3.92E-04 [pb]']
#     Conditions violation (if any) =  {'None': None}
#     Mass of cluster =  [[8.65E+02 [GeV], 1.29E+02 [GeV]], [8.65E+02 [GeV], 1.29E+02 [GeV]]]
# 

# In[7]:



