
## How To: Compare theory predictions with experimental limits

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


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
smsHelpers.base=os.path.join(os.getenv("HOME"),"smodels-database/")
listofanalyses = smsAnalysisFactory.load()
listOfTopologies = slhaDecomposer.decompose (filename, sigcut = 0.03 * fb, doCompress=True, doInvisible=True,minmassgap = 5* GeV)
analysesPredictions = [theoryPredictionFor(analysis, listOfTopologies) for analysis in listofanalyses]


# Out[4]:

#     12:23:58.229 INFO     smodels.experiment.smsAnalysisFactory:39  useSuperseded is not set, skipping superseded results
#     12:23:58.279 INFO     smodels.theory.slhaDecomposer:132 Ignoring t+ decays
#     12:23:58.286 INFO     smodels.theory.slhaDecomposer:132 Ignoring higgs decays
#     12:23:58.287 INFO     smodels.theory.slhaDecomposer:132 Ignoring H0 decays
#     12:23:58.287 INFO     smodels.theory.slhaDecomposer:132 Ignoring A0 decays
#     12:23:58.287 INFO     smodels.theory.slhaDecomposer:132 Ignoring H+ decays
#     12:23:58.320 INFO     smodels.theory.crossSection:512 Ignoring 76 lower order cross-sections
#     12:24:12.437 INFO     smodels.experiment.smsInterpolation:174 Masses out of range for ATLAS-SUSY-2013-19/T6bbWW (no extrapolation)
#     12:24:12.502 INFO     smodels.experiment.smsInterpolation:174 Masses out of range for CMS-SUS-13-011/T6bbWW (no extrapolation)
#     12:24:12.942 INFO     smodels.experiment.smsInterpolation:174 Masses out of range for CMS-PAS-SUS-13-008/T6ttWW (no extrapolation)
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

#     Analysis name =  CMS-SUS-12-028:T2
#     Theory prediction =  1.77E-03 [pb]
#     Upper limit =  1.64E-02 [pb]
#     Analysis name =  CMS-SUS-12-028:T1
#     Theory prediction =  3.92E-04 [pb]
#     Upper limit =  3.14E-02 [pb]
#     Analysis name =  ATLAS-SUSY-2013-12:TChiWZ
#     Theory prediction =  1.85E-02 [pb]
#     Upper limit =  3.30E-01 [pb]
#     Analysis name =  ATLAS-SUSY-2013-11:TChiWZ
#     Theory prediction =  1.85E-02 [pb]
#     Upper limit =  9.22E-01 [pb]
#     Analysis name =  CMS-SUS-13-012:T2
#     Theory prediction =  1.77E-03 [pb]
#     Upper limit =  6.10E-03 [pb]
#     Analysis name =  CMS-SUS-13-012:T1
#     Theory prediction =  3.92E-04 [pb]
#     Upper limit =  8.03E-03 [pb]
#     Analysis name =  CMS-PAS-SUS-13-019:T2
#     Theory prediction =  1.77E-03 [pb]
#     Upper limit =  3.76E-03 [pb]
#     Analysis name =  CMS-PAS-SUS-13-019:T1
#     Theory prediction =  3.92E-04 [pb]
#     Upper limit =  9.69E-03 [pb]
#     Analysis name =  ATLAS-SUSY-2013-02:T2
#     Theory prediction =  1.77E-03 [pb]
#     Upper limit =  6.10E-03 [pb]
#     Analysis name =  ATLAS-SUSY-2013-02:T1
#     Theory prediction =  3.92E-04 [pb]
#     Upper limit =  1.08E-02 [pb]
#     Analysis name =  CMS-SUS-13-006:TChiWZ
#     Theory prediction =  1.85E-02 [pb]
#     Upper limit =  4.64E-01 [pb]
# 
