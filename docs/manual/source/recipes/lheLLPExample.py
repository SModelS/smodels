#!/usr/bin/env python
# coding: utf-8

# # How To: Run SModelS with LHE input with additional width information

# In[1]:


# Set up the path to SModelS installation folder
import sys; sys.path.append("."); import smodels_paths


# In[2]:


from smodels.tools import runtime
#Define your model
runtime.modelFile = 'smodels.share.models.mssm' 

from smodels.theory import decomposer
from smodels.tools.physicsUnits import fb, GeV, TeV
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObj import Database
from smodels.tools import coverage
from smodels.tools.smodelsLogging import setLogLevel
from smodels.particlesLoader import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
setLogLevel("info")


# ### Read LHE input:

# In[8]:


# Path to input file (either a SLHA or LHE file)
lhefile = 'inputFiles/lhe/simplyGluino.lhe'
model = Model(BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=lhefile)


# ### Replace widths

# In[9]:


# At this point all BSM particles in model have either zero (stable) or infinite (unstable) widths.
# However the widths can be overwritten for the desired (long-lived) particles.
# For instance, the gluino width can be set to 1e-15 GeV as below:
gluino = model.getParticlesWith(pdg = 1000021)[0]
print('old width=',gluino.totalwidth) #Check that the width was indeed infinite
#Assign new width value:
gluino.totalwidth = 1e-15*GeV
print('new width=',gluino.totalwidth)


# In[10]:


#The same has to be done for the anti-particles:
gluino = model.getParticlesWith(pdg = -1000021)[0]
gluino.totalwidth = 1e-15*GeV
#The remaining steps are as for any input:


# ### Decompose the input model:

# In[12]:


# Set main options for decomposition
sigmacut = 0.01 * fb
mingap = 5. * GeV
# Decompose model (use slhaDecomposer for SLHA input or lheDecomposer for LHE input)
toplist = decomposer.decompose(model, sigmacut, doCompress=True, doInvisible=True, minmassgap=mingap)

# Access basic information from decomposition, using the topology list and topology objects:
print( "\n Decomposition Results: " )
print( "\t  Total number of topologies: %i " %len(toplist) )
nel = sum([len(top.elementList) for top in toplist])
print( "\t  Total number of elements = %i " %nel )


# ### Load the Database of experimental results:

# In[13]:


# Set the path to the database
database = Database("official")
# Load the experimental results to be used.
# In this case, all results are employed.
listOfExpRes = database.getExpResults()

# Print basic information about the results loaded.
# Count the number of loaded UL and EM experimental results:
nUL, nEM = 0, 0
for exp in listOfExpRes:
    expType = exp.getValuesFor('dataType')[0]
    if expType == 'upperLimit':
        nUL += 1
    elif  expType == 'efficiencyMap':
        nEM += 1
print("\n Loaded Database with %i UL results and %i EM results " %(nUL,nEM))


# ### Match the decomposed simplified models with the experimental database of constraints:

# In[14]:


# Compute the theory predictions for each experimental result and print them:
print("\n Theory Predictions and Constraints:")
rmax = 0.
bestResult = None
for expResult in listOfExpRes:
    predictions = theoryPredictionsFor(expResult, toplist)
    if not predictions: continue # Skip if there are no constraints from this result
    print('\n %s (%i TeV)' %(expResult.globalInfo.id,expResult.globalInfo.sqrts.asNumber(TeV)))
    for theoryPrediction in predictions:
        dataset = theoryPrediction.dataset
        datasetID = theoryPrediction.dataId()
        mass = theoryPrediction.mass
        txnames = [str(txname) for txname in theoryPrediction.txnames]
        PIDs =  theoryPrediction.PIDs         
        print( "------------------------" )
        print( "TxNames = ",txnames )  
        print( "Theory Prediction = ",theoryPrediction.xsection.value )  #Signal cross section
        # Get the corresponding upper limit:
        print( "UL for theory prediction = ",theoryPrediction.upperLimit )
        # Compute the r-value
        r = theoryPrediction.getRValue()
        print( "r = ",r )
        #Compute likelihhod and chi^2 for EM-type results:
        if dataset.dataInfo.dataType == 'efficiencyMap':
            theoryPrediction.computeStatistics()
            print( 'Chi2, likelihood=', theoryPrediction.chi2, theoryPrediction.likelihood )
        if r > rmax:
            rmax = r
            bestResult = expResult.globalInfo.id

# Print the most constraining experimental result
print( "\nThe largest r-value (theory/upper limit ratio) is ",rmax )
if rmax > 1.:
    print( "(The input model is likely excluded by %s)" %bestResult )
else:
    print( "(The input model is not excluded by the simplified model results)" )


# ### Check for simplified models in the input model which were not tested by the Database
# ### (result with displaced gluino decays will now show up)

# In[16]:


#Find out missing topologies for sqrts=8*TeV:
uncovered = coverage.Uncovered(toplist,sqrts=8.*TeV)
#First sort coverage groups by label
groups = sorted(uncovered.groups[:], key = lambda g: g.label)
#Print uncovered cross-sections:
for group in groups:
    print("\nTotal cross-section for %s (fb): %10.3E\n" %(group.description,group.getTotalXSec()))


# In[17]:


missingTopos = uncovered.getGroup('missing (prompt)')
#Print some of the missing topologies:
if missingTopos.generalElements:
    print('Missing topologies (up to 3):' )
    for genEl in missingTopos.generalElements[:3]:
        print('Element:', genEl)
        print('\tcross-section (fb):', genEl.missingX)
else:
    print("No missing topologies found\n")

missingDisplaced = uncovered.getGroup('missing (displaced)')
#Print elements with displaced decays:
if missingDisplaced.generalElements:
    print('\nElements with displaced vertices (up to 2):' )    
    for genEl in missingDisplaced.generalElements[:2]:
        print('Element:', genEl)
        print('\tcross-section (fb):', genEl.missingX)
else:
    print("\nNo displaced decays")


# In[ ]:




