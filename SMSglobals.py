#Define some global variables
from SMSHelpers import addunit

ListOfAnalyses=[]
DoCompress = False     #Flag to turn on/off compression
DoInvisible = False    #Flag to turn on/off invisible compression

def initglob():
    global ncomp         #Overall counter for compressed topologies
    global min_lep_pt, min_tau_pt, min_jet_pt, min_b_pt, min_pt
    global evcount,nComp,nInvis
    global minmassgap
    global DistAnalyses

    evcount = 0
    nComp = 0              #Counter for mass compressed topologies
    nInvis = 0             #Counter for invisible compressed topologies
    minmassgap = addunit(10.,'GeV')        #Minimum mass gap (for mass compression)
    DistAnalyses = ""      #List of analyses used to compute mass distances (used by MassDist function)
    
    ListOfAnalyses = []  #Initialize list of analyses
