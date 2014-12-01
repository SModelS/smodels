
# coding: utf-8

## How To: Find missing topologies that are not covered by the database

# In[1]:

# Set up the path to SModelS installation folder if running on a different folder
import sys
sys.path.append("../")


# In[2]:

# Import those parts of smodels that are needed for this exercise
from smodels.tools.physicsUnits import TeV, GeV, fb


# In[3]:

from smodels.installation import installDirectory


# In[4]:

from smodels.theory import slhaDecomposer


# In[5]:

from smodels.experiment import smsAnalysisFactory, smsHelpers


# In[6]:

from smodels.tools import missingTopologies


# In[7]:

# define where the database resides
smsHelpers.base = "../test/database/"


# In[8]:

# load list of analyses from database
listOfAnalyses = smsAnalysisFactory.load()


# In[9]:

# Define the SLHA file name
filename = "%s/inputFiles/slha/gluino_squarks.slha" % installDirectory()


# In[10]:

# Perform the decomposition:
listOfTopologies = slhaDecomposer.decompose (filename, sigcut=0.5*fb, doCompress=True, doInvisible=True, minmassgap=5*GeV)


# In[11]:

# Initiate missing Topologies for 8 TeV
missingtopos = missingTopologies.MissingTopoList(8*TeV)


# In[12]:

# Check listOfTopologies against listOfAnalyses to find missing topologies
missingtopos.findMissingTopos(listOfTopologies, listOfAnalyses, minmassgap=5*GeV,doCompress=True, doInvisible=True)


# In[13]:

# to print a sorted list of missing topologies with high weights use
missingtopos.printout()


# In[14]:

# To access the missing topologies direcly
# For the i-th entry, where the entries are not sorted, do
i = 3
topology = missingtopos.topos[i]
print topology.topo
print topology.weights
print topology.value

