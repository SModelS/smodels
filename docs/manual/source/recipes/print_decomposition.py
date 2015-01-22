
## How To: Print out the theoretical decomposition

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


# In[2]:

#Import those parts of smodels that are needed for this exercise
#(We will assume the input is a SLHA file. For LHE files, use the lheDecomposer instead)
from smodels.theory import slhaDecomposer
from smodels.installation import installDirectory
from smodels.tools.physicsUnits import fb, GeV


# In[3]:

#Define the SLHA file name
filename="%s/inputFiles/slha/gluino_squarks.slha" % installDirectory()


# In[4]:

#Perform the decomposition:
listOfTopologies = slhaDecomposer.decompose (filename, sigcut = 0.5 * fb, doCompress=True, doInvisible=True,minmassgap = 5* GeV)


# Out[4]:

#     12:21:45.668 INFO     smodels.theory.slhaDecomposer:132 Ignoring t+ decays
#     12:21:45.674 INFO     smodels.theory.slhaDecomposer:132 Ignoring higgs decays
#     12:21:45.674 INFO     smodels.theory.slhaDecomposer:132 Ignoring H0 decays
#     12:21:45.674 INFO     smodels.theory.slhaDecomposer:132 Ignoring A0 decays
#     12:21:45.674 INFO     smodels.theory.slhaDecomposer:132 Ignoring H+ decays
#     12:21:45.708 INFO     smodels.theory.crossSection:512 Ignoring 76 lower order cross-sections
# 

# In[5]:

#Print a summary of all the topologies generated:
listOfTopologies.printout()


# Out[5]:

#        ======================================================= 
#      || 	 						 || 
#      || 	 	 Global topologies table 	 	 ||
#      || 	 						 || 
#        ======================================================= 
#     ===================================================== 
#     Topology:
#     Number of vertices: [3, 2] 
#     Number of vertex parts: [[1, 1, 0], [1, 0]]
#     Total Global topology weight :
#     Sqrts: 8.00E+00 [TeV], Weight:6.72E-03 [pb]
#     
#     Total Number of Elements: 7
#     ===================================================== 
#     Topology:
#     Number of vertices: [3, 3] 
#     Number of vertex parts: [[1, 1, 0], [1, 1, 0]]
#     Total Global topology weight :
#     Sqrts: 8.00E+00 [TeV], Weight:1.16E-02 [pb]
#     
#     Total Number of Elements: 7
#     ===================================================== 
#     Topology:
#     Number of vertices: [2, 2] 
#     Number of vertex parts: [[2, 0], [1, 0]]
#     Total Global topology weight :
#     Sqrts: 8.00E+00 [TeV], Weight:1.63E-03 [pb]
#     
#     Total Number of Elements: 1
#     ===================================================== 
#     Topology:
#     Number of vertices: [3, 2] 
#     Number of vertex parts: [[2, 1, 0], [1, 0]]
#     Total Global topology weight :
#     Sqrts: 8.00E+00 [TeV], Weight:1.66E-02 [pb]
#     
#     Total Number of Elements: 7
#     ===================================================== 
#     Topology:
#     Number of vertices: [1, 2] 
#     Number of vertex parts: [[0], [1, 0]]
#     Total Global topology weight :
#     Sqrts: 8.00E+00 [TeV], Weight:6.15E-04 [pb]
#     
#     Total Number of Elements: 1
#     ===================================================== 
#     Topology:
#     Number of vertices: [2, 2] 
#     Number of vertex parts: [[1, 0], [1, 0]]
#     Total Global topology weight :
#     Sqrts: 8.00E+00 [TeV], Weight:1.76E-01 [pb]
#     
#     Total Number of Elements: 7
#     ===================================================== 
#     Topology:
#     Number of vertices: [3, 2] 
#     Number of vertex parts: [[1, 1, 0], [2, 0]]
#     Total Global topology weight :
#     Sqrts: 8.00E+00 [TeV], Weight:3.49E-03 [pb]
#     
#     Total Number of Elements: 2
#     ===================================================== 
#     Topology:
#     Number of vertices: [3, 3] 
#     Number of vertex parts: [[1, 1, 0], [2, 1, 0]]
#     Total Global topology weight :
#     Sqrts: 8.00E+00 [TeV], Weight:4.03E-02 [pb]
#     
#     Total Number of Elements: 17
#     	 .................................................. 
#     Number of vertex parts: [[2, 1, 0], [2, 1, 0]]
#     Total Global topology weight :
#     Sqrts: 8.00E+00 [TeV], Weight:2.06E-02 [pb]
#     
#     Total Number of Elements: 10
#     
# 

# In[6]:

#To print specific information about othe i-th topology:
i = 3
top = listOfTopologies[i]
print "Number of vertices = ",top.vertnumb
print "Number of final states = ",top.vertparts
print "Number of elements = ",len(top.elementList)


# Out[6]:

#     Number of vertices =  [3, 2]
#     Number of final states =  [[2, 1, 0], [1, 0]]
#     Number of elements =  7
# 

# In[7]:

#We can also print information for each element in the list:
for element in top.elementList:
    element.printout()    


# Out[7]:

#     		 Particles in element: [[['b', 't-'], ['W+']], [['jet']]]
#     		 The element masses are 
#     		 Branch 0: [8.65E+02 [GeV], 2.69E+02 [GeV], 1.29E+02 [GeV]]
#     		 Branch 1: [9.91E+02 [GeV], 1.29E+02 [GeV]]
#     		 The element weights are: 
#      		 Sqrts: 8.00E+00 [TeV], Weight:1.98E-03 [pb]
#     
#     		 Particles in element: [[['t+', 'b'], ['W-']], [['jet']]]
#     		 The element masses are 
#     		 Branch 0: [8.65E+02 [GeV], 2.69E+02 [GeV], 1.29E+02 [GeV]]
#     		 Branch 1: [9.91E+02 [GeV], 1.29E+02 [GeV]]
#     		 The element weights are: 
#      		 Sqrts: 8.00E+00 [TeV], Weight:1.98E-03 [pb]
#     
#     		 Particles in element: [[['jet', 'jet'], ['W+']], [['jet']]]
#     		 The element masses are 
#     		 Branch 0: [8.65E+02 [GeV], 2.69E+02 [GeV], 1.29E+02 [GeV]]
#     		 Branch 1: [9.91E+02 [GeV], 1.29E+02 [GeV]]
#     		 The element weights are: 
#      		 Sqrts: 8.00E+00 [TeV], Weight:3.91E-03 [pb]
#     
#     		 Particles in element: [[['jet', 'jet'], ['W-']], [['jet']]]
#     		 The element masses are 
#     		 Branch 0: [8.65E+02 [GeV], 2.69E+02 [GeV], 1.29E+02 [GeV]]
#     		 Branch 1: [9.91E+02 [GeV], 1.29E+02 [GeV]]
#     		 The element weights are: 
#      		 Sqrts: 8.00E+00 [TeV], Weight:3.91E-03 [pb]
#     
#     		 Particles in element: [[['b', 'b'], ['higgs']], [['jet']]]
#     		 The element masses are 
#     		 Branch 0: [8.65E+02 [GeV], 2.69E+02 [GeV], 1.29E+02 [GeV]]
#     		 Branch 1: [9.91E+02 [GeV], 1.29E+02 [GeV]]
#     		 The element weights are: 
#      		 Sqrts: 8.00E+00 [TeV], Weight:1.06E-03 [pb]
#     
#     		 Particles in element: [[['jet', 'jet'], ['higgs']], [['jet']]]
#     		 The element masses are 
#     		 Branch 0: [8.65E+02 [GeV], 2.69E+02 [GeV], 1.29E+02 [GeV]]
#     		 Branch 1: [9.91E+02 [GeV], 1.29E+02 [GeV]]
#     		 The element weights are: 
#      		 Sqrts: 8.00E+00 [TeV], Weight:3.28E-03 [pb]
#     
#     		 Particles in element: [[['t+', 't-'], ['higgs']], [['jet']]]
#     		 The element masses are 
#     		 Branch 0: [8.65E+02 [GeV], 2.69E+02 [GeV], 1.29E+02 [GeV]]
#     		 Branch 1: [9.91E+02 [GeV], 1.29E+02 [GeV]]
#     		 The element weights are: 
#      		 Sqrts: 8.00E+00 [TeV], Weight:5.15E-04 [pb]
#     
# 

# In[8]:

#The element information can be also accessed directly:
el = top.elementList[0]
print "Final states = ",el.getParticles()
print "Mass array = ",el.getMasses()
print "Weight = ",el.weight.niceStr()


# Out[8]:

#     Final states =  [[['b', 't-'], ['W+']], [['jet']]]
#     Mass array =  [[8.65E+02 [GeV], 2.69E+02 [GeV], 1.29E+02 [GeV]], [9.91E+02 [GeV], 1.29E+02 [GeV]]]
#     Weight =  Sqrts: 8.00E+00 [TeV], Weight:1.98E-03 [pb]
#     
# 
