
## How To: Look up the upper limit of a particular result, for a particular set masses

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


# In[2]:

#Import those parts of smodels that are needed for this exercise
from smodels.tools.physicsUnits import GeV
from smodels.experiment import smsAnalysisFactory, smsHelpers


# In[3]:

## define where the database resides
smsHelpers.base=os.path.join(os.getenv("HOME"),"smodels-database/")
#specify analysis and topology as strings:
#and load analysis
ananame = ["CMS-SUS-12-028"]
txname = ["T2"]
analysis = smsAnalysisFactory.load(ananame, topologies=txname)[0]


# Out[3]:

#     12:14:47.608 INFO     smodels.experiment.smsAnalysisFactory:39  useSuperseded is not set, skipping superseded results
# 

# In[4]:

masses = [[500*GeV, 150*GeV],[500*GeV, 150*GeV]]


# In[5]:

analysis.getUpperLimitFor(masses)


# Out[5]:

#     1.28E-01 [pb]

# In[5]:



