#Define some global variables

def initglob():
    global ncomp         #Overall counter for compressed topologies
    global min_lep_pt, min_tau_pt, min_jet_pt, min_b_pt, min_pt
    global ListOfAnalyses
    global AnalysesRes
    global evcount
    global massequiv,maxgap
    evcount = 0
    ncomp = [0,0]
    min_pt = 0.         #Minimum pT for any particle (= 0 to turn off compression)
    min_lep_pt = 0.     #Minimum lepton pT (= 0 to turn off compression)
    min_tau_pt = 0.     #Minimum tau pT (= 0 to turn off compression)
    min_jet_pt = 0.     #Minimum jet pT (= 0 to turn off compression)
    min_b_pt = 0.       #Minimum b pT (= 0 to turn off compression)
    massequiv = 0.05     #Maximum difference (in %) for the definition of equal masses
    maxgap = 100.         #Max absolute gap between masses for the def. of equal masses
    
    ListOfAnalyses = []  #Initialize list of analyses
    AnalysesRes = []   #Initialize dictionary for all analyses
    
