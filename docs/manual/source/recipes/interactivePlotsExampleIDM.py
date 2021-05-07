#!/usr/bin/env python
# coding: utf-8

# ### Import the plot creator

# In[1]:


# Set up the path to SModelS installation folder
import sys,os

sys.path.append(".")
import smodels_paths
from smodels.tools.interactivePlots import PlotMaster
from smodels.installation import installDirectory
import smodels


# ### Define the input folders, the parameters file and where to store the plots

# In[2]:


smodelsOutput = f'{installDirectory()}/inputFiles/scanExample/smodels-output.tar.gz'
slhaFiles = f'{installDirectory()}/inputFiles/scanExample/slhas.tar.gz'
parameters = './iplots_parameters_IDM.py'
indexfile = 'main.html'
modelfile=f'{installDirectory()}/smodels/share/models/idm.py'
outputFolder = './'
npoints = -1 #If negative will use all the points in the folders


# ### Create plots

# In[3]:


plotMaster = PlotMaster(smodelsOutput,slhaFiles,parameters,indexfile)
loadData = plotMaster.loadData(npoints)
plotMaster.makePlots(outputFolder)  


# ### Display the output

# In[4]:


from IPython.core.display import display,HTML


# In[5]:


#Main page (to open the plots go to the plots output folder and open them with the browser)
mainFile = open(outputFolder + indexfile,'r')
display(HTML(mainFile.read()))
mainFile.close()

