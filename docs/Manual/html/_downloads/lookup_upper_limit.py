
# coding: utf-8

## How To: Look up the upper limit of a particular result, for a particular set masses

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys
sys.path.append("../")


# In[2]:

#Import those parts of smodels that are needed for this exercise
from smodels.tools.physicsUnits import GeV
from smodels.experiment import smsAnalysisFactory, smsHelpers


# In[3]:

## define where the database resides
smsHelpers.base="../test/database/"
#specify analysis and topology as strings:
#and load analysis
ananame = ["SUS12028"]
txname = ["T2"]
analysis = smsAnalysisFactory.load(ananame, topologies=txname)[0]


# In[4]:

masses = [[500*GeV, 150*GeV],[500*GeV, 150*GeV]]


# In[5]:

analysis.getUpperLimitFor(masses)


# In[5]:



