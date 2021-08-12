#!/usr/bin/env python
# coding: utf-8

# ### Import the plot creator

# In[1]:


# Set up the path to SModelS installation folder
import sys,os
sys.path.append(".")

import smodels_paths
from smodels.tools.interactivePlots import Plotter
from smodels.installation import installDirectory


# ### Define the input folders, the parameters file and where to store the plots

# In[2]:


smodelsOutput = './inputFiles/scanExample/smodels-output.tar.gz'
slhaFiles = './inputFiles/scanExample/slhas.tar.gz'
parameters = './iplots_parameters.py'
indexfile = 'main.html'
modelfile=f'{installDirectory()}/smodels/share/models/mssm.py'
outputFolder = './'
npoints = -1 #If negative will use all the points in the folders


# ### Create plots

# In[3]:


plotter = Plotter(smodelsOutput,slhaFiles,parameters,modelfile)
plotter.loadData(npoints)
plotter.plot(outputFolder,indexfile)


# ### Display the output

# In[4]:


plotter.display()

