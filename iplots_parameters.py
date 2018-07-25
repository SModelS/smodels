
#Plot title
plot_title= 'My interactive plot example'

#Label, block and code number for the variables you want to plot, e.g. 'm_gluino':['MASS',1000021]. These will be your x and y axes.
#Note that spaces are not allowed in variable's names. 

variable_x= {'m_gluino[GeV]':['MASS',1000021]}
variable_y= {'m_suR[GeV]':['MASS',2000002]}

#SLHA hover information: In a dictionary form, give the name* of your variable, the block and code number to find it in the SLHA file.
slha_hover_information={'m_gluino': ['MASS',1000021],'m_suR':['MASS',2000002],'m_LSP': ['MASS',1000022]} 

#For which particles you want to get the mean decay length.
ctau_hover_information={'ctau_gluino':1000021}
#For which particles you want to display the decay channels and branching ratios. 
BR_hover_information={'BR_gluino':1000021}
#WARNING: Lists of branching ratios lists can be very long, so the may not fit in the hover box. 
#You can tell how many entries you want to print with BR_get_top, e.g. BR_get_top=5 (default 'all').
BR_get_top='all'


#SModelS hover information; options are:
#
#SmodelS_excluded -> prints whether the point is excluded or not by SModelS
#r_max -> shows the highest r-value for each parameter point
#chi2 -> shows the chi^2 value, if available (if not, the output is 'none')
#Tx -> shows the topology/ies which give r_max
#Analysis -> shows the experimental analysis from which the strongest constraint (r_max) comes from
#MT_max -> shows the missing topology with the largest cross section (in SModelS bracket notation)
#MT_max_xsec -> shows the cross section of MT_max
#MT_total_xsec -> shows the total missing cross section (i.e. the sum of all missing topologies cross sections)  
#MT_long_xsec -> shows the total missing cross section in long cascade decays  
#MT_asym_xsec -> shows the total missing cross section in decays with asymmetric branches 
#MT_outgrid_xsec -> shows the total missing cross section outside the mass grids of the experimental results
#file -> shows the name of the input spectrum file 

smodels_hover_information=['SmodelS_excluded','r_max','chi2','Tx','Analysis','MT_max','MT_max_xsec','file','MT_total_xsec']

#Set which plots do you want. 
#Same options as for SModels hover information are available.
plot_data=['all','non-excluded','excluded'] # choice: all, non-excluded, excluded points
plot_list=['SmodelS_excluded','r_max','chi2','Tx','Analysis','MT_max','MT_max_xsec','MT_long_xsec','MT_asym_xsec','MT_total_xsec']
