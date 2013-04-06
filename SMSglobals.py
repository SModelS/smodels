#Define some global variables

def initglob():
    global ncomp         #Overall counter for compressed topologies
    global min_lep_pt, min_tau_pt, min_jet_pt, min_b_pt, min_pt
    global ListOfAnalyses
    global evcount,nComp,nInvis
    global minmassgap
    global DoCompress, DoInvisible

    DoCompress = False     #Flag to turn on/off compression
    DoInvisible = False    #Flag to turn on/off invisible compression
    evcount = 0
    nComp = 0              #Counter for mass compressed topologies
    nInvis = 0             #Counter for invisible compressed topologies
    minmassgap = 1.         #Minimum mass gap (for compression)
    
    ListOfAnalyses = []  #Initialize list of analyses

    
