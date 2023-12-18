#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os,sys
# Set up the path to SModelS installation folder
import sys; sys.path.append("."); import smodels_paths
import numpy as np
import matplotlib.pyplot as plt

from smodels.tools.modelTester import getCombiner


# ### Define input (SLHA) file and the parameters file

# In[2]:


slhafile = os.path.expanduser("./inputFiles/slha/ewino_example.slha")
# Define parameters file with combineAnas = ATLAS-SUSY-2019-08,ATLAS-SUSY-2019-09:
parfile = os.path.expanduser("./parameters_comb.ini")


# ### Define some basic parameters for plotting the likelihoods

# In[3]:


expected = False # whether to plot the observed or expected likelihood
normalize = True # whether to normalize the likelihoods
muvals = np.linspace(0.,5.,200) # Signal strength values for which to evaluate the likelihoods


# ### Run SModelS and get the analysis combination

# In[4]:


combiner = getCombiner(slhafile, parfile)


# ### Use the combination to evaluate the likelihoods

# In[5]:


llhdDict = combiner.getLlhds(muvals,expected,normalize)


# ### Compute L_SM, L_BSM, L_max and the UL on mu

# In[6]:


muhat = combiner.muhat()
lmax = combiner.lmax()
lsm = combiner.lsm()
lbsm = combiner.likelihood(mu=1.0)
muULDict = {'combined' : combiner.getUpperLimitOnMu()}
for theoryPred in combiner.theoryPredictions:
    anaID = theoryPred.analysisId()
    muULDict[anaID] = theoryPred.getUpperLimitOnMu()


# ### Plot the results

# In[7]:


fig = plt.figure(figsize=(8,4))
ymin = 0.
ymax = 0.
for anaID,l in llhdDict.items():
    if anaID == 'combined':
        zorder = 100
        linestyle = '--'
    else:
        zorder = None
        linestyle = '-'
        
    p = plt.plot(muvals,l,label=anaID,linewidth=3,zorder=zorder,linestyle=linestyle)
    ymin = min(ymin,min(l))
    ymax = max(ymax,max(l))
    
    plt.vlines(muULDict[anaID],ymin=ymin,ymax=ymax,linestyle='--',linewidth=2,
           label=r'$\mu_{UL}$ (%s)' %anaID,color=p[0].get_color(),alpha=0.7)
    
    
plt.vlines(muhat,ymin=ymin,ymax=ymax,linestyle='--',linewidth=2,
           label=r'$\hat{\mu}$',color='black',alpha=0.7)


# plt.yscale('log')
# plt.ylim(1e-1,1e1)
plt.xlim(muvals.min(),muvals.max())
plt.xlabel(r'$\mu$')
if normalize:
    plt.ylabel('Normalized Likelihood')
else:
    plt.ylabel('Likelihood')
plt.legend(framealpha=1)
plt.title(r'$\hat{\mu} = $ %1.2f, $L_{max} =$ %1.2e, $L_{SM} =$ %1.2e, $L_{BSM} =$ %1.2e' %(muhat,lmax,lsm,lbsm))
plt.tight_layout()
plt.show()


# In[ ]:




