
#Plot title
plot_title= 'My plot'

#Name* block and code number of the variable you want to plot, e. g. 'm_gluino':['MASS',1000021]. These will be your x and y axis. 

variable_x= {'m_gluino[GeV]':['MASS',1000021]}
variable_y= {'m_suR[GeV]':['MASS',2000002]}


#SLHA Hover information: In a dictionary form, give the name* of your variable, the block and code number to find it in the SLHA file.

slha_hover_information={'m_gluino': ['MASS',1000021],'m_suR':['MASS',2000002],'m_LSP': ['MASS',1000022]} 

#For which particles you want to get its decay width and mean decay length.
ctau_hover_information={'ctau_suR':2000002}
#For which particles you want to get its branching ratio. 
BR_hover_information={'BR_gluino':1000021}
#WARNING: Branching ratio's lists can be very long, so the may not fit in the hover,but you can tell how many branches you want to print with BR_get_top (default 'all').
BR_get_top='all'


#*Espaces are not allowed in variable's names.

#Smodels hover information. Options are:
#SModelS hover information. 
#Options are:
#SmodelS_excluded -> Prints whether the point is excluded or not by SModelS
#r_max -> Shows the highest r-value for each parameter point
#chi2 -> Shows the chi^2 value,if available (if not the output is -.1).
#Tx -> Shows the topologies which give r_max
#Analysis -> Shows the experimental analysis from which the strongest constraint (r_max) comes from
#MT_max -> Shows the missing topology with the largest cross section (in SModelS bracket notation)
#MT_max_xsec -> Shows the cross section of MT_max
#MT_total_mxsec -> Shows the total missing cross section (i.e. the sum of all missing topologies cross sections)  
#MT_long_xsec -> Shows the total missing cross section in long cascade decays  
#MT_asym_xsec -> Shows the total missing cross section in decays with asymmetric branches 
#MT_outgrid_xsec -> Shows the total missing cross section where we are outside the mass grid 
#file -> Shows the name of the input spectrum file 
smodels_hover_information=['SmodelS_excluded','chi2','r_max','Tx','Analysis','MT_asym_xsec','MT_total_mxsec','MT_outgrid_xsec','file','MT_long_xsec']

#Set which plots do you want. 
#Same options from Smodels hover information are available.
plot_data=['all','non-excluded','excluded'] # choice: all, non-excluded, excluded points
plot_list=['SmodelS_excluded','chi2','r_max','Tx','Analysis','MT_asym_xsec','MT_total_mxsec','MT_outgrid_xsec','MT_max','MT_long_xsec']
