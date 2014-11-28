
# coding: utf-8

## How To: Load the database, selecting only a few results.

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys
sys.path.append("../")


# In[2]:

from smodels.experiment import smsAnalysisFactory, smsHelpers
from smodels.tools.physicsUnits import GeV


# In[3]:

## define where the database resides
smsHelpers.base="../test/database/"


### How to load results from one publication (or conference note)

# In[4]:

#Select only the CMS SUS-12-028 conference note
analyses=["SUS12028"]


# In[5]:

#Loads the selected analyses
list_of_analyses=smsAnalysisFactory.load(analyses)


# In[6]:

#Print the analyses that were loaded:
for analysis in list_of_analyses: analysis.printout(outputLevel=1)


# In[7]:

#To see which elements are constrained by the analyses (in bracket notation), set outputLevel=2
for analysis in list_of_analyses: analysis.printout(outputLevel=2)


# In[8]:

#To print basic information about one analysis:
analysisT1, analysisT2 = list_of_analyses
print "Name = ",analysisT1.label,", Sqrts = ",analysisT1.sqrts, ", Luminosity =",analysisT1.lum


# In[9]:

#To obtain the upper limit for a given analysis and a given mass vector.
#Note that the number of masses in the mass vector must be consitent with the analysis. For the T1 analysis, for instance:
massesT1 = [[300*GeV,100*GeV],[300*GeV,100*GeV]]
analysisT1 = list_of_analyses[0]
print analysisT1.getUpperLimitFor(massesT1)


# In[10]:

#For the T2 analysis:
massesT2 = [[300*GeV,50*GeV],[300*GeV,50*GeV]]
print analysisT2.getUpperLimitFor(massesT2)


# In[11]:

#If you try with the wrong mass format, an error will be printed:
masses = [[300*GeV],[300*GeV,50*GeV]]
print analysisT2.getUpperLimitFor(masses)


### How to load results for one constraint (Txname)

# In[12]:

#It is also possible to load all the results for a single constraint (using the Txname convention)
#(The warning are just because we are using a test databse)
Txnames = ["T1"]
new_list = smsAnalysisFactory.load(topologies=Txnames)


# In[13]:

#In this "test" example there is only one analysis:
for analysis in new_list: print analysis.label


### How to load all experimental analyses, including the superseded publications

# In[14]:

#By default only non-supersed analyses are loaded:
analysisList = smsAnalysisFactory.load()
for analysis in analysisList: print analysis.label


# In[15]:

#To load all analyses (included the superseded ones), set useSuperseded=True
full_list = smsAnalysisFactory.load(useSuperseded=True)
for analysis in full_list: print analysis.label

