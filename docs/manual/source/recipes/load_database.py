
## How To: Load the database, selecting only a few results.

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


# In[2]:

from smodels.experiment import smsAnalysisFactory, smsHelpers
from smodels.tools.physicsUnits import GeV


# In[3]:

## define where the database resides
smsHelpers.base=os.path.join(os.getenv("HOME"),"smodels-database/")


### How to load results from one publication (or conference note)

# In[4]:

#Select only the CMS SUS-12-028 conference note
analyses=["CMS-SUS-12-028"]


# In[5]:

#Loads the selected analyses
#(The INFO tells you that superseded analyses are not loaded, see below)
list_of_analyses=smsAnalysisFactory.load(analyses)


# Out[5]:

#     12:12:28.879 INFO     smodels.experiment.smsAnalysisFactory:39  useSuperseded is not set, skipping superseded results
# 

# In[6]:

#Print the analyses that were loaded:
for analysis in list_of_analyses: analysis.printout(outputLevel=1)


# Out[6]:

#     ========================================================
#     Analysis Name: CMS-SUS-12-028
#     Tx Label: T1bbbb
#     Analysis Sqrts: 8.00E+00 [TeV]
#     
#     ========================================================
#     Analysis Name: CMS-SUS-12-028
#     Tx Label: T1tttt
#     Analysis Sqrts: 8.00E+00 [TeV]
#     
#     ========================================================
#     Analysis Name: CMS-SUS-12-028
#     Tx Label: T2
#     Analysis Sqrts: 8.00E+00 [TeV]
#     
#     ========================================================
#     Analysis Name: CMS-SUS-12-028
#     Tx Label: T2bb
#     Analysis Sqrts: 8.00E+00 [TeV]
#     
#     ========================================================
#     Analysis Name: CMS-SUS-12-028
#     Tx Label: T1
#     Analysis Sqrts: 8.00E+00 [TeV]
#     
# 

# In[7]:

#To see which elements are constrained by the analyses (in bracket notation), set outputLevel=2
for analysis in list_of_analyses: analysis.printout(outputLevel=2)


# Out[7]:

#     ========================================================
#     Analysis Name: CMS-SUS-12-028
#     Tx Label: T1bbbb
#     Analysis Sqrts: 8.00E+00 [TeV]
#     	 -----------------------------
#     	 Elements tested by analysis:
#     	    [[['b', 'b']], [['b', 'b']]]
#     
#     ========================================================
#     Analysis Name: CMS-SUS-12-028
#     Tx Label: T1tttt
#     Analysis Sqrts: 8.00E+00 [TeV]
#     	 -----------------------------
#     	 Elements tested by analysis:
#     	    [[['t', 't']], [['t', 't']]]
#     
#     ========================================================
#     Analysis Name: CMS-SUS-12-028
#     Tx Label: T2
#     Analysis Sqrts: 8.00E+00 [TeV]
#     	 -----------------------------
#     	 Elements tested by analysis:
#     	    [[['jet']], [['jet']]]
#     
#     ========================================================
#     Analysis Name: CMS-SUS-12-028
#     Tx Label: T2bb
#     Analysis Sqrts: 8.00E+00 [TeV]
#     	 -----------------------------
#     	 Elements tested by analysis:
#     	    [[['b']], [['b']]]
#     
#     ========================================================
#     Analysis Name: CMS-SUS-12-028
#     Tx Label: T1
#     Analysis Sqrts: 8.00E+00 [TeV]
#     	 -----------------------------
#     	 Elements tested by analysis:
#     	    [[['jet', 'jet']], [['jet', 'jet']]]
#     
# 

# In[8]:

#To print basic information about one analysis:
analysisT1 = list_of_analyses[4]
analysisT2 = list_of_analyses[2]
print "Name = ",analysisT1.label,", Sqrts = ",analysisT1.sqrts, ", Luminosity =",analysisT1.lum


# Out[8]:

#     Name =  CMS-SUS-12-028:T1 , Sqrts =  8.00E+00 [TeV] , Luminosity = 1.17E+01 [1/fb]
# 

# In[9]:

#To obtain the upper limit for a given analysis and a given mass vector.
#Note that the number of masses in the mass vector must be consitent with the analysis. For the T1 analysis, for instance:
massesT1 = [[300*GeV,100*GeV],[300*GeV,100*GeV]]
analysisT1 = list_of_analyses[0]
print analysisT1.getUpperLimitFor(massesT1)


# Out[9]:

#     6.27E-01 [pb]
# 

# In[10]:

#For the T2 analysis:
massesT2 = [[300*GeV,50*GeV],[300*GeV,50*GeV]]
print analysisT2.getUpperLimitFor(massesT2)


# Out[10]:

#     1.07E+00 [pb]
# 

# In[11]:

#If you try with the wrong mass format, an error will be printed:
masses = [[300*GeV],[300*GeV,50*GeV]]
print analysisT2.getUpperLimitFor(masses)


# Out[11]:

#     12:12:29.789 ERROR    smodels.experiment.limitGetter:76  Masses differ between branches.
#     False
# 

### How to load results for one constraint (Txname)

# In[12]:

#It is also possible to load all the results for a single constraint (using the Txname convention)
Txnames = ["T1"]
new_list = smsAnalysisFactory.load(topologies=Txnames)


# Out[12]:

#     12:12:29.923 INFO     smodels.experiment.smsAnalysisFactory:39  useSuperseded is not set, skipping superseded results
# 

# In[13]:

#Print all the analyses containing the required Txname:
for analysis in new_list: print analysis.label


# Out[13]:

#     CMS-SUS-12-028:T1
#     CMS-SUS-13-012:T1
#     CMS-PAS-SUS-13-019:T1
#     ATLAS-SUSY-2013-02:T1
# 

### How to load all experimental analyses, including the superseded publications

# In[14]:

#By default only non-supersed analyses are loaded:
analysisList = smsAnalysisFactory.load()
for analysis in analysisList: print analysis.label


# Out[14]:

#     12:12:30.224 INFO     smodels.experiment.smsAnalysisFactory:39  useSuperseded is not set, skipping superseded results
#     CMS-SUS-12-024:T1bbbb
#     CMS-SUS-12-024:T1ttttoff
#     CMS-SUS-12-024:T1tttt
#     CMS-SUS-12-024:T5tttt
#     ATLAS-CONF-2013-065:T2tt
#     CMS-SUS-12-028:T1bbbb
#     CMS-SUS-12-028:T1tttt
#     CMS-SUS-12-028:T2
#     CMS-SUS-12-028:T2bb
#     CMS-SUS-12-028:T1
#     ATLAS-CONF-2012-105:T1tttt
#     ATLAS-CONF-2013-007:T1tttt
#     ATLAS-SUSY-2013-12:TChiWH
#     ATLAS-SUSY-2013-12:TChiWZ
#     ATLAS-SUSY-2013-12:TChiWZoff
#     ATLAS-SUSY-2013-11:TChiWZ
#     ATLAS-SUSY-2013-11:TSlepSlep
#     ATLAS-SUSY-2013-15:T2tt
#     ATLAS-SUSY-2013-15:T2bbWW
#     ATLAS-SUSY-2013-14:TStauStau
#     CMS-SUS-13-012:T1tttt
#     CMS-SUS-13-012:T2
#     CMS-SUS-13-012:T1ttttoff
#     CMS-SUS-13-012:T1
#     CMS-SUS-13-013:T1tttt
#     CMS-SUS-13-013:T1ttttoff
#     ATLAS-SUSY-2013-19:T2tt
#     ATLAS-SUSY-2013-19:T2bbWW
#     ATLAS-SUSY-2013-19:T6bbWW
#     CMS-SUS-13-011:T2tt
#     CMS-SUS-13-011:T6bbWW
#     ATLAS-CONF-2013-024:T2tt
#     CMS-PAS-SUS-13-019:T1tttt
#     CMS-PAS-SUS-13-019:T2tt
#     CMS-PAS-SUS-13-019:T1bbbb
#     CMS-PAS-SUS-13-019:T2
#     CMS-PAS-SUS-13-019:T1
#     CMS-PAS-SUS-13-019:T1ttttoff
#     CMS-PAS-SUS-13-019:T2bb
#     CMS-SUS-13-002:T1tttt
#     CMS-PAS-SUS-13-018:T2bb
#     ATLAS-SUSY-2013-04:T1tttt
#     ATLAS-SUSY-2013-05:T2bb
#     ATLAS-SUSY-2013-02:T2
#     ATLAS-SUSY-2013-02:T1
#     CMS-PAS-SUS-13-008:T6ttWW
#     CMS-PAS-SUS-13-008:T1tttt
#     CMS-SUS-13-007:T1tttt
#     CMS-SUS-13-007:T1ttttoff
#     CMS-SUS-13-006:TChiChipmSlepL
#     CMS-SUS-13-006:TChiWZ
#     CMS-SUS-13-006:TChiChipmSlepStau
#     CMS-SUS-13-006:TSlepSlep
#     CMS-SUS-13-006:TChiWZoff
#     CMS-PAS-SUS-14-011:T1bbbb
#     CMS-PAS-SUS-14-011:T2tt
#     CMS-PAS-SUS-14-011:T1tttt
#     CMS-PAS-SUS-14-011:T1ttttoff
#     CMS-PAS-SUS-13-016:T1tttt
#     CMS-PAS-SUS-13-016:T1ttttoff
#     ATLAS-CONF-2013-061:T1bbbb
#     ATLAS-CONF-2013-061:T1tttt
# 

# In[15]:

#To load all analyses (included the superseded ones), set useSuperseded=True
full_list = smsAnalysisFactory.load(useSuperseded=True)
for analysis in full_list: print analysis.label


# Out[15]:

#     ATLAS-CONF-2013-049:TSlepSlep
#     ATLAS-CONF-2013-048:T2bbWW
#     ATLAS-CONF-2013-048:T6bbWW
#     CMS-SUS-12-024:T1bbbb
#     CMS-SUS-12-024:T1ttttoff
#     CMS-SUS-12-024:T1tttt
#     CMS-SUS-12-024:T5tttt
#     ATLAS-CONF-2013-065:T2tt
#     CMS-SUS-12-028:T1bbbb
#     CMS-SUS-12-028:T1tttt
#     CMS-SUS-12-028:T2
#     CMS-SUS-12-028:T2bb
#     CMS-SUS-12-028:T1
#     ATLAS-CONF-2013-047:T2
#     ATLAS-CONF-2013-047:T1
#     ATLAS-CONF-2012-105:T1tttt
#     ATLAS-CONF-2013-007:T1tttt
#     ATLAS-SUSY-2013-12:TChiWH
#     ATLAS-SUSY-2013-12:TChiWZ
#     ATLAS-SUSY-2013-12:TChiWZoff
#     ATLAS-SUSY-2013-11:TChiWZ
#     ATLAS-SUSY-2013-11:TSlepSlep
#     ATLAS-SUSY-2013-15:T2tt
#     ATLAS-SUSY-2013-15:T2bbWW
#     ATLAS-SUSY-2013-14:TStauStau
#     CMS-SUS-13-012:T1tttt
#     CMS-SUS-13-012:T2
#     CMS-SUS-13-012:T1ttttoff
#     CMS-SUS-13-012:T1
#     CMS-SUS-13-013:T1tttt
#     CMS-SUS-13-013:T1ttttoff
#     ATLAS-SUSY-2013-19:T2tt
#     ATLAS-SUSY-2013-19:T2bbWW
#     ATLAS-SUSY-2013-19:T6bbWW
#     CMS-SUS-13-011:T2tt
#     CMS-SUS-13-011:T6bbWW
#     ATLAS-CONF-2013-024:T2tt
#     CMS-PAS-SUS-13-019:T1tttt
#     CMS-PAS-SUS-13-019:T2tt
#     CMS-PAS-SUS-13-019:T1bbbb
#     CMS-PAS-SUS-13-019:T2
#     CMS-PAS-SUS-13-019:T1
#     CMS-PAS-SUS-13-019:T1ttttoff
#     CMS-PAS-SUS-13-019:T2bb
#     ATLAS-CONF-2012-166:T2tt
#     CMS-SUS-13-002:T1tttt
#     CMS-PAS-SUS-12-026:T1tttt
#     CMS-PAS-SUS-13-018:T2bb
#     CMS-PAS-SUS-12-022:TChiChipmSlepL
#     CMS-PAS-SUS-12-022:TChiWZ
#     CMS-PAS-SUS-12-022:TSlepSlep
#     CMS-PAS-SUS-12-022:TChiChipmSlepStau
#     ATLAS-SUSY-2013-04:T1tttt
#     ATLAS-SUSY-2013-05:T2bb
#     ATLAS-CONF-2013-053:T2bb
#     ATLAS-SUSY-2013-02:T2
#     ATLAS-SUSY-2013-02:T1
#     CMS-PAS-SUS-13-004:T2tt
#     CMS-PAS-SUS-13-004:T1bbbb
#     CMS-PAS-SUS-13-004:T1ttttoff
#     CMS-PAS-SUS-13-004:T1tttt
#     CMS-PAS-SUS-13-008:T6ttWW
#     CMS-PAS-SUS-13-008:T1tttt
#     ATLAS-CONF-2013-037:T2tt
#     ATLAS-CONF-2013-037:T6bbWW
#     ATLAS-CONF-2013-035:TChiWZ
#     ATLAS-CONF-2013-035:TChiWZoff
#     CMS-SUS-13-007:T1tttt
#     CMS-SUS-13-007:T1ttttoff
#     CMS-SUS-13-006:TChiChipmSlepL
#     CMS-SUS-13-006:TChiWZ
#     CMS-SUS-13-006:TChiChipmSlepStau
#     CMS-SUS-13-006:TSlepSlep
#     CMS-SUS-13-006:TChiWZoff
#     CMS-PAS-SUS-14-011:T1bbbb
#     CMS-PAS-SUS-14-011:T2tt
#     CMS-PAS-SUS-14-011:T1tttt
#     CMS-PAS-SUS-14-011:T1ttttoff
#     CMS-PAS-SUS-13-016:T1tttt
#     CMS-PAS-SUS-13-016:T1ttttoff
#     ATLAS-CONF-2013-061:T1bbbb
#     ATLAS-CONF-2013-061:T1tttt
# 

# In[15]:



