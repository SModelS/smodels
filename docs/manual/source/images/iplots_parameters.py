#Note: html notation like <sub>NAME</sub> for subscript and <sup>NAME</sup> for superscript 
#can be used in axes labels and spectrum hover information (slha blocks, BRs, ctau,...), 
#but not in the SModelS hover information.


#Plot title
plot_title = 'interactive plot example from a small pMSSM scan'

#Label*, block and code number for the variables you want to plot, e.g. 'm(gluino)': ['MASS', 1000021]. These will be your x and y axes.
variable_x = {'m<sub>gluino</sub>': ['MASS', 1000021]}
variable_y = {'m<sub>&#967;<sub>1</sub><sup>0</sup></sub>': ['MASS', 1000022]}


#SLHA hover information: In a dictionary form, give the name* of your variable, the block and code number to find it in the SLHA file.
slha_hover_information = {'m(gluino)</sub>': ['MASS', 1000021], 'm(stop1)': ['MASS', 1000006], 'm(chi10)': ['MASS', 1000022]} 

#For which particles you want to get the mean decay length.
ctau_hover_information = {'ctau(chi1+)': 1000024}
#For which particles you want to display the decay channels and branching ratios. 
BR_hover_information = {'BR(gluino)': 1000021, 'BR(chi1+)': 1000024}
#The output is written in the form '.25[1000022,1,-1]',  where the first number (0.25) is the branching ratio, and the numbers in [,] are the PDG codes of the decay products.
#WARNING: Lists of branching ratios lists can be very long, so the may not fit in the hover box. 
#You can tell how many entries you want to print with BR_get_top, e.g. BR_get_top = 5 (default 'all').
BR_get_top = 5


#SModelS hover information; options are:
#
#SModelS_status -> prints whether the point is excluded or not by SModelS.
#r_max -> shows the highest r-value for each parameter point.
#chi2 -> shows the chi^2 value which corresponds to r_max, if available (if not, the output is 'False').
#Tx -> shows the topology/ies which give r_max
#Analysis -> shows the experimental analysis from which the strongest constraint (r_max) comes from.
#MT_max -> shows the missing topology with the largest cross section (in SModelS bracket notation).
#MT_max_xsec -> shows the cross section of MT_max.
#MT_total_xsec -> shows the total missing cross section (i.e. the sum of all missing topologies cross sections) . 
#MT_long_xsec -> shows the total missing cross section in long cascade decays . 
#MT_asym_xsec -> shows the total missing cross section in decays with asymmetric branches. 
#MT_outgrid_xsec -> shows the total missing cross section outside the mass grids of the experimental results.
#file -> shows the name of the input spectrum file. 

SModelS_hover_information = ['SModelS_status', 'r_max', 'Tx', 'Analysis','chi2', 'MT_max', 'MT_max_xsec', 'MT_total_xsec', 'MT_long_xsec', 'MT_asym_xsec', 'MT_outgrid_xsec', 'file']

#Set which plots do you want
# choice for plot data: all, non-excluded, excluded points.
plot_data = ['all', 'non-excluded', 'excluded'] 
# choice for plot list: same options as for SModels hover information.
plot_list = ['SModelS_status','r_max', 'chi2', 'Tx', 'Analysis', 'MT_max', 'MT_max_xsec', 'MT_total_xsec', 'MT_long_xsec', 'MT_asym_xsec','MT_outgrid_xsec']
