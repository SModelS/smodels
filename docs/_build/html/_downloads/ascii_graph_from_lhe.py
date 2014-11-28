
# coding: utf-8

## How To: Create ASCII graphs for elements in a LHE event file

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys
sys.path.append("../")


# In[2]:

#Import those parts of smodels that are needed for this exercise
from smodels.theory import lheReader, lheDecomposer, crossSection
from smodels.installation import installDirectory
from smodels.tools import asciiGraph


# In[3]:

#Load an input file containing LHE events and start the LHE reader
filename="%s/inputFiles/lhe/gluino_squarks.lhe" % installDirectory()
reader = lheReader.LheReader ( filename )


# In[4]:

#Read the next event and generate the corresponding element
event=reader.next()
element=lheDecomposer.elementFromEvent (event)


# In[5]:

#Print the corresponding ASCII graph
print asciiGraph.asciidraw ( element )


# In[6]:

#Do the same for the next event:
event=reader.next()
element=lheDecomposer.elementFromEvent ( event )
print asciiGraph.asciidraw ( element )

