
## How To: Compute NLL cross sections for a given SLHA file

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys
sys.path.append("../")


# In[2]:

#Import those parts of smodels that are needed for this exercise
from smodels.tools import xsecComputer
from smodels.tools.physicsUnits import TeV, fb
from smodels.installation import installDirectory
from smodels.tools.xsecComputer import LO, NLL


# In[3]:

#Define the SLHA file name
filename="%s/inputFiles/slha/lightSquarks.slha" % installDirectory()


# In[4]:

#Lets compute the NLL cross-sections for 8 TeV. The xsecComputer will first use Pythia to compute
#the LO cross-sections and then NLLfast to compute the k-factors.
#The output will contain only the processes contained in NLLfast (gluino and squark production)
#For the Pythia step we have to define the number of MC events (1k)
#(To see how to also include LO cross-sections for all process check the LO HowTo)
xsecsNLL=xsecComputer.computeXSec( 8*TeV, NLL, 1000, filename)


# In[5]:

# the output is a XSectionList ...
type(xsecsNLL)


# Out[5]:

#     smodels.theory.crossSection.XSectionList

# In[6]:

#Each entry in the list contains the cross-section value:
print(xsecsNLL[0].value)
#The PDGs of the particles produced:
print(xsecsNLL[0].pid)
#And some additional info
print("label =",xsecsNLL[0].info.label,"Sqrts =",xsecsNLL[0].info.sqrts, "QCD order =",xsecsNLL[0].info.order)


# Out[6]:

#     2.22E-02 [pb]
#     (1000001, 1000021)
#     ('label =', '8 TeV (NLO+NLL)', 'Sqrts =', 8.00E+00 [TeV], 'QCD order =', 2)
# 

# In[7]:

#Finally, lets write the cross-sections back to the file (will write only the cross-sections not overlapping the existing ones):
xsecComputer.addXSecToFile(xsecsNLL,filename)


# Out[7]:

#     12:10:10.253 WARNING  smodels.tools.xsecComputer:150 SLHA file already contains XSECTION blocks. Adding only missing cross-sections.
# 

#     True

# In[7]:



