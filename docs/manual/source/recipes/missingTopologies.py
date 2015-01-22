
## How To: Find missing topologies that are not covered by the database

# In[1]:

# Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


# In[2]:

# Import those parts of smodels that are needed for this exercise
from smodels.tools.physicsUnits import TeV, GeV, fb
from smodels.installation import installDirectory
from smodels.theory import slhaDecomposer
from smodels.experiment import smsAnalysisFactory, smsHelpers
from smodels.tools import missingTopologies


# In[3]:

# define where the database resides
smsHelpers.base=os.path.join(os.getenv("HOME"),"smodels-database/")


# In[4]:

# load list of analyses from database
listOfAnalyses = smsAnalysisFactory.load()


# Out[4]:

#     12:25:29.818 INFO     smodels.experiment.smsAnalysisFactory:39  useSuperseded is not set, skipping superseded results
# 

# In[5]:

# Define the SLHA file name
filename = "%s/inputFiles/slha/gluino_squarks.slha" % installDirectory()


# In[6]:

# Perform the decomposition:
listOfTopologies = slhaDecomposer.decompose (filename, sigcut=0.5*fb, doCompress=True, doInvisible=True, minmassgap=5*GeV)


# Out[6]:

#     12:25:33.310 INFO     smodels.theory.slhaDecomposer:132 Ignoring t+ decays
#     12:25:33.315 INFO     smodels.theory.slhaDecomposer:132 Ignoring higgs decays
#     12:25:33.316 INFO     smodels.theory.slhaDecomposer:132 Ignoring H0 decays
#     12:25:33.316 INFO     smodels.theory.slhaDecomposer:132 Ignoring A0 decays
#     12:25:33.316 INFO     smodels.theory.slhaDecomposer:132 Ignoring H+ decays
#     12:25:33.348 INFO     smodels.theory.crossSection:512 Ignoring 76 lower order cross-sections
# 

# In[7]:

# Initiate missing Topologies for 8 TeV
missingtopos = missingTopologies.MissingTopoList(8*TeV)


# In[8]:

# Check listOfTopologies against listOfAnalyses to find missing topologies
missingtopos.findMissingTopos(listOfTopologies, listOfAnalyses, minmassgap=5*GeV,doCompress=True, doInvisible=True)


# In[9]:

# to print a sorted list of missing topologies with high weights use
missingtopos.printout()


# Out[9]:

#     
#     ================================================================================
#     Missing topologies with the highest cross-sections (up to 10):
#     Sqrts (TeV)   Weight (fb)        Element description
#     8.00E+00   5.958E+01    #                                 [[[W]],[[W]]]
#     8.00E+00   1.511E+01    #                 [[[jet],[W]],[[jet,jet],[W]]]
#     8.00E+00   9.069E+00    #             [[[jet,jet],[W]],[[jet,jet],[W]]]
#     8.00E+00   7.822E+00    #                     [[[jet]],[[jet,jet],[W]]]
#     8.00E+00   5.500E+00    #                     [[[jet],[W]],[[jet],[W]]]
#     8.00E+00   4.964E+00    #             [[[jet],[W]],[[jet,jet],[higgs]]]
#     8.00E+00   4.932E+00    #             [[[jet],[higgs]],[[jet,jet],[W]]]
#     8.00E+00   4.603E+00    #                 [[[b,t],[W]],[[jet,jet],[W]]]
#     8.00E+00   4.603E+00    #                 [[[jet,jet],[W]],[[t,b],[W]]]
#     8.00E+00   4.596E+00    #                 [[[jet],[W]],[[jet],[higgs]]]
#     
# 

# In[10]:

# To access the missing topologies direcly
# For the i-th entry, where the entries are not sorted, do
i = 3
topology = missingtopos.topos[i]
print topology.topo
print topology.weights
print topology.value


# Out[10]:

#     [[[jet]],[[jet],[W]]]
#     ['8.00E+00 [TeV]:2.11E-03 [pb]']
#     2.11E+00 
# 
