def load():
    from SMSmethods import EAnalysis
    
    ListOfAnalyses = []
    #List analyses and classify in topologies and elements
        
     #Global Topology:
    Analysis = EAnalysis()
    Analysis.label = "T1-simple" #it can be whatever you want
    Analysis.Top.vertnumb = [2,2]    #Number of vertices 
    Analysis.Top.vertparts = [[2,0],[2,0]] #Number of particle insertions 
    Analysis.results = {
#T1    
    "[[[jet,jet]],[[jet,jet]]]" : "None",
#T1tttt    
    "[[[t,t]],[[t,t]]]" : "None",
#T1bbbb    
    "[[[b,b]],[[b,b]]]" : "None",
#T1tbtb    
    "[[[t,b]],[[t,b]]]" : "None",
#T1ttttpm    
    "[[[t+,t-]],[[t+,t-]]]" : "None"     
    }
    Analysis.plots  = {
    "[[[jet,jet]],[[jet,jet]]]" : ["T1",["alphaT","alphaT8TeV"]],
    "[[[t,t]],[[t,t]]]" : ["T1tttt",["alphaT","RA48TeV","RA2b8TeV","MultiLepton8TeV","alphaT8TeV","SUS13008"]],
    "[[[b,b]],[[b,b]]]" : ["T1bbbb",["alphaT","RA2b8TeV","alphaT8TeV"]],
    "[[[t,b]],[[t,b]]]" : ["T1tbtb",["ATLAS_CONF_2013_007"]],
    "[[[t+,t-]],[[t+,t-]]]" : ["T1tttt",["ATLAS_CONF_2012_105","ATLAS_CONF_2013_007"]]
    }    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)
    


#     #Global Topology:
    Analysis = EAnalysis()
    Analysis.label = "T2-simple" #it can be whatever you want
    Analysis.Top.vertnumb = [2,2]    #Number of vertices 
    Analysis.Top.vertparts = [[1,0],[1,0]] #Number of particle insertions 
    Analysis.results = {
#T2    
    "[[[jet]],[[jet]]]" : "None",
#T2tt    
    "[[[t]],[[t]]]" : "None",
#T2bb    
    "[[[b]],[[b]]]" : "None",
#TChiWZ    
    "[[[W]],[[Z]]]" : "None",
#TChiZZ    
#    "[[Z],[Z]]" : "None",
#TSlepSlep    
    "[[[l+]],[[l-]]]" : "None"
    

    }
    Analysis.plots  = {
    "[[[jet]],[[jet]]]" : ["T2",["alphaT","alphaT8TeV"]],
    "[[[t]],[[t]]]" : ["T2tt",["alphaT","LeptonicStop8TeV","alphaT8TeV","ATLAS_CONF_2013_024","ATLAS_CONF_2013_037","ATLAS_CONF_2012_166"]],
    "[[[b]],[[b]]]" : ["T2bb",["alphaT","alphaT8TeV","ATLAS_CONF_2013_001"]],
    "[[[W]],[[Z]]]" : ["TChiwz",["Weakinos8TeV","ATLAS_CONF_2013_035"]],
    "[[[l+]],[[l-]]]" : ["TSlepSlep",["Weakinos8TeV"]]
    }    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)
    


    
#         #Global Topology:
    Analysis = EAnalysis()
    Analysis.label = "T6-type" #it can be whatever you want
    Analysis.Top.vertnumb = [3,3]    #Number of vertices 
    Analysis.Top.vertparts = [[1,1,0],[1,1,0]] #Number of particle insertions 
    Analysis.results = {
#TChiChiSlep (Democratic chargino and neutralino decays):    
    "2.*([[[L],[L]],[[L],[nu]]] + [[[L],[L]],[[nu],[L]]])" : "Csim([[[L],[L]],[[L],[nu]]],[[[L],[L]],[[nu],[L]]]); Cgtr([[[L],[L]],[[L],[nu]]],3.*[[[ta],[ta]],[[L],[nu]]]); Cgtr([[[L],[L]],[[L],[nu]]], 3.*[[[L],[L]],[[ta],[nu]]]); Cgtr([[[L],[L]],[[nu],[L]]],3.*[[[ta],[ta]],[[nu],[L]]]); Cgtr([[[L],[L]],[[nu],[L]]],3.*[[[L],[L]],[[nu],[ta]]])",
#TChiChiSlep (tau-enriched chargino decay):    
    "[[[L],[L]],[[nu],[ta]]]" : "Cgtr([[[L],[L]],[[nu],[ta]]], 3.*[[[ta],[ta]],[[nu],[ta]]])",
#TChiChiSlep (tau-dominated  chargino and neutralino decays):
    "[[[ta],[ta]],[[nu],[ta]]]" : "None",
#Chargino-Chargino:
    "[[[L+],[nu]],[[L-],[nu]]] + [[[L+],[nu]],[[nu],[L-]]] + [[[L-],[nu]],[[nu],[L+]]] + [[[nu],[L+]],[[nu],[L-]]]" : "Csim([[[L+],[nu]],[[L-],[nu]]],[[[L+],[nu]],[[nu],[L-]]], [[[L-],[nu]],[[nu],[L+]]], [[[nu],[L+]],[[nu],[L-]]]); Cgtr([[[L+],[nu]],[[L-],[nu]]], 3.*[[[ta+],[nu]],[[L-],[nu]]]); Cgtr([[[L+],[nu]],[[L-],[nu]]], 3.*[[[L+],[nu]],[[ta-],[nu]]]); Cgtr([[[L+],[nu]],[[nu],[L-]]], 3.*[[[ta+],[nu]],[[nu],[L-]]]); Cgtr([[[L+],[nu]],[[nu],[L-]]], 3.*[[[L+],[nu]],[[nu],[ta-]]]); Cgtr([[[L-],[nu]],[[nu],[L+]]], 3.*[[[ta-],[nu]],[[nu],[L+]]]); Cgtr([[[L-],[nu]],[[nu],[L+]]], 3.*[[[L-],[nu]],[[nu],[ta+]]]); Cgtr([[[nu],[L+]],[[nu],[L-]]], 3.*[[[nu],[ta+]],[[nu],[L-]]]); Cgtr([[[nu],[L+]],[[nu],[L-]]], 3.*[[[nu],[L+]],[[nu],[ta-]]])",
#T6bbww (stop pair production)
    "[[[b],[W]],[[b],[W]]]" : "None",
#T5tttt
    "[[[t+],[t-]],[[t+],[t-]]] + [[[t-],[t+]],[[t+],[t-]]] + [[[t-],[t+]],[[t-],[t+]]]" : "None",
#T5tttt
    "[[[t],[t]],[[t],[t]]]]" : "None",    
#T6ttWW
    "[[[t+],[W-]],[[t+],[W-]]] + [[[t-],[W+]],[[t+],[W-]]] + [[[t-],[W+]],[[t-],[W+]]]" : "None",
#T6ttWW
    "[[[t],[W]],[[t],[W]]]" : "None",    
#T6ttzz    
    "[[[Z],[t]],[[Z],[t]]]" : "None",
#T6bbzz    
    "[[[b],[Z]],[[b],[Z]]]" : "None",    
#TChiChiSlepSlep
    "[[[l],[l]],[[l],[l]]]" : "Cgtr([[[l],[l]],[[l],[l]]],2.*[[[e],[l]],[[l],[l]]]); Cgtr([[[l],[l]],[[l],[l]]],2.*[[[l],[e]],[[l],[l]]])",
#TChiChipmStauL    
    "2.*([[[ta+],[ta-]],[[nu],[ta]]] + [[[ta+],[ta-]],[[ta],[nu]]] + [[[ta-],[ta+]],[[nu],[ta]]] + [[[ta-],[ta+]],[[ta],[nu]]])" : "Csim([[[ta+],[ta-]],[[nu],[ta]]],[[[ta+],[ta-]],[[ta],[nu]]],[[[ta-],[ta+]],[[nu],[ta]]],[[[ta-],[ta+]],[[ta],[nu]]])",
#TChipChimStauSnu
    "[[[ta+],[nu]],[[ta-],[nu]]] + [[[ta+],[nu]],[[nu],[ta-]]] + [[[ta-],[nu]],[[nu],[ta+]]] + [[[nu],[ta+]],[[nu],[ta-]]]" : "Csim([[[ta+],[nu]],[[ta-],[nu]]],[[[ta+],[nu]],[[nu],[ta-]]],[[[ta-],[nu]],[[nu],[ta+]]],[[[nu],[ta+]],[[nu],[ta-]]])"
    }
    Analysis.plots  = {   
    "2.*([[[L],[L]],[[L],[nu]]] + [[[L],[L]],[[nu],[L]]])" : ["TChiChipmSlepL",["Weakinos8TeV","ATLAS_CONF_2013_035"]],
    "[[[L],[L]],[[nu],[ta]]]" : ["TChiChipmSlepStau",["Weakinos8TeV"]],
    "[[[ta],[ta]],[[nu],[ta]]]" : ["TChiChipmStauStau",["Weakinos8TeV"]],
    "[[[L+],[nu]],[[L-],[nu]]] + [[[L+],[nu]],[[nu],[L-]]] + [[[L-],[nu]],[[nu],[L+]]] + [[[nu],[L+]],[[nu],[L-]]]" : ["TChipChimSlepSnu",["Weakinos8TeV"]],
    "[[[b],[W]],[[b],[W]]]" : ["T6bbWW",["LeptonicStop8TeV","ATLAS_CONF_2013_037","ATLAS_CONF_2012_166"]],
    "[[[t+],[t-]],[[t+],[t-]]] + [[[t-],[t+]],[[t+],[t-]]] + [[[t-],[t+]],[[t-],[t+]]]" : ["T5tttt",["ATLAS_CONF_2013_007"]],
    "[[[t],[t]],[[t],[t]]]]" : ["T5tttt",["SUS13008"]],
    "[[[t+],[W-]],[[t+],[W-]]] + [[[t-],[W+]],[[t+],[W-]]] + [[[t-],[W+]],[[t-],[W+]]]" : ["T6ttWW",["ATLAS_CONF_2013_007"]],
    "[[[t],[W]],[[t],[W]]]" : ["T6ttWW",["SUS13008"]],
    "[[[Z],[t]],[[Z],[t]]]" : ["T6ttZZ",["ATLAS_CONF_2013_025"]],
    "[[[b],[Z]],[[b],[Z]]]" : ["T6bbZZ",["SUS13008"]],
    "[[[l],[l]],[[l],[l]]]" : ["TChiChiSlepSlep",["ATLAS_CONF_2013_036"]],
    "2.*([[[ta+],[ta-]],[[nu],[ta]]] + [[[ta+],[ta-]],[[ta],[nu]]] + [[[ta-],[ta+]],[[nu],[ta]]] + [[[ta-],[ta+]],[[ta],[nu]]])" : ["TChiChipmStauL",["ATLAS_CONF_2013_028"]],
    "[[[ta+],[nu]],[[ta-],[nu]]] + [[[ta+],[nu]],[[nu],[ta-]]] + [[[ta-],[nu]],[[nu],[ta+]]] + [[[nu],[ta+]],[[nu],[ta-]]]" : ["TChipChimStauSnu",["ATLAS_CONF_2013_028"]]
    }    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)
    
    
#         #Global Topology:
    Analysis = EAnalysis()
    Analysis.label = "T5-type" #it can be whatever you want
    Analysis.Top.vertnumb = [3,3]    #Number of vertices 
    Analysis.Top.vertparts = [[2,1,0],[2,1,0]] #Number of particle insertions 
    Analysis.results = {
#T5WW
    "[[[jet,jet],[W]],[[jet,jet],[W]]]" : "None"        
    }
    Analysis.plots  = {   
    "[[[jet,jet],[W]],[[jet,jet],[W]]]" : ["T5WW",["ATLAS_CONF_2013_007"]]  
    }    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)
    
    
#         #Global Topology:
    Analysis = EAnalysis()
    Analysis.label = "T8-type" #it can be whatever you want
    Analysis.Top.vertnumb = [4,4]    #Number of vertices 
    Analysis.Top.vertparts = [[1,1,1,0],[1,1,1,0]] #Number of particle insertions 
    Analysis.results = {
#T8ChiSlep
    "1.454545([[[jet],[nu],[L]],[[jet],[nu],[L]]] + [[[jet],[nu],[L]],[[jet],[L],[nu]]] + [[[jet],[nu],[L]],[[jet],[L],[L]]] + [[[jet],[L],[nu]],[[jet],[L],[nu]]] + [[[jet],[L],[nu]],[[jet],[L],[L]]] + [[[jet],[L],[L]],[[jet],[L],[L]]] + [[[jet],[L],[L]],[[jet],[nu],[nu]]])" : "Csim([[[jet],[L],[nu]],[[jet],[L],[nu]]], [[[jet],[nu],[L]],[[jet],[nu],[L]]], [[[jet],[nu],[nu]],[[jet],[L],[L]]]); Csim([[[jet],[L],[nu]],[[jet],[nu],[L]]], 2.*([[[jet],[nu],[L]],[[jet],[nu],[L]]])); Cgtr(2.*([[[jet],[L],[L]],[[jet],[L],[L]]]), [[[jet],[nu],[L]],[[jet],[L],[L]]]); Cgtr(2.*([[[jet],[L],[L]],[[jet],[L],[L]]]), [[[jet],[L],[nu]],[[jet],[L],[L]]]); Cgtr([[[jet],[L],[nu]],[[jet],[L],[L]]], [[[jet],[nu],[L]],[[jet],[L],[nu]]]); Cgtr([[[jet],[L],[nu]],[[jet],[L],[L]]], [[[jet],[nu],[L]],[[jet],[nu],[L]]]); Cgtr([[[jet],[L],[nu]],[[jet],[L],[L]]], [[[jet],[L],[nu]],[[jet],[L],[nu]]]); Cgtr([[[jet],[L],[nu]],[[jet],[L],[L]]], [[[jet],[L],[L]],[[jet],[nu],[nu]]])",
    "[[[b],[t],[W]],[[b],[t],[W]]]" : "None"
    }
    Analysis.plots  = {   
    "1.454545([[[jet],[nu],[L]],[[jet],[nu],[L]]] + [[[jet],[nu],[L]],[[jet],[L],[nu]]] + [[[jet],[nu],[L]],[[jet],[L],[L]]] + [[[jet],[L],[nu]],[[jet],[L],[nu]]] + [[[jet],[L],[nu]],[[jet],[L],[L]]] + [[[jet],[L],[L]],[[jet],[L],[L]]] + [[[jet],[L],[L]],[[jet],[nu],[nu]]])" : ["T8ChiSlep",["ATLAS_CONF_2013_007"]],
    "[[[b],[t],[W]],[[b],[t],[W]]]" : ["T7btbtWW",["SUS13008"]]
    }    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)   
    
    
#         #Global Topology:
    Analysis = EAnalysis()
    Analysis.label = "T7-type" #it can be whatever you want
    Analysis.Top.vertnumb = [4,4]    #Number of vertices 
    Analysis.Top.vertparts = [[2,1,1,0],[2,1,1,0]] #Number of particle insertions 
    Analysis.results = {
#T7ChiSlep
     "1.454545([[[jet,jet],[nu],[L]],[[jet,jet],[nu],[L]]] + [[[jet,jet],[nu],[L]],[[jet,jet],[L],[nu]]] + [[[jet,jet],[nu],[L]],[[jet,jet],[L],[L]]] + [[[jet,jet],[L],[nu]],[[jet,jet],[L],[nu]]] + [[[jet,jet],[L],[nu]],[[jet,jet],[L],[L]]] + [[[jet,jet],[L],[L]],[[jet,jet],[L],[L]]] + [[[jet,jet],[L],[L]],[[jet,jet],[nu],[nu]]])" : "Csim([[[jet,jet],[L],[nu]],[[jet,jet],[L],[nu]]], [[[jet,jet],[nu],[L]],[[jet,jet],[nu],[L]]], [[[jet,jet],[nu],[nu]],[[jet,jet],[L],[L]]]); Csim([[[jet,jet],[L],[nu]],[[jet,jet],[nu],[L]]], 2.*([[[jet,jet],[nu],[L]],[[jet,jet],[nu],[L]]])); Cgtr(2.*([[[jet,jet],[L],[L]],[[jet,jet],[L],[L]]]), [[[jet,jet],[nu],[L]],[[jet,jet],[L],[L]]]); Cgtr(2.*([[[jet,jet],[L],[L]],[[jet,jet],[L],[L]]]), [[[jet,jet],[L],[nu]],[[jet,jet],[L],[L]]]); Cgtr([[[jet,jet],[L],[nu]],[[jet,jet],[L],[L]]], [[[jet,jet],[nu],[L]],[[jet,jet],[L],[nu]]]); Cgtr([[[jet,jet],[L],[nu]],[[jet,jet],[L],[L]]], [[[jet,jet],[nu],[L]],[[jet,jet],[nu],[L]]]); Cgtr([[[jet,jet],[L],[nu]],[[jet,jet],[L],[L]]], [[[jet,jet],[L],[nu]],[[jet,jet],[L],[nu]]]); Cgtr([[[jet,jet],[L],[nu]],[[jet,jet],[L],[L]]], [[[jet,jet],[L],[L]],[[jet,jet],[nu],[nu]]])"       
    }
    Analysis.plots  = {   
    "1.454545([[[jet,jet],[nu],[L]],[[jet,jet],[nu],[L]]] + [[[jet,jet],[nu],[L]],[[jet,jet],[L],[nu]]] + [[[jet,jet],[nu],[L]],[[jet,jet],[L],[L]]] + [[[jet,jet],[L],[nu]],[[jet,jet],[L],[nu]]] + [[[jet,jet],[L],[nu]],[[jet,jet],[L],[L]]] + [[[jet,jet],[L],[L]],[[jet,jet],[L],[L]]] + [[[jet,jet],[L],[L]],[[jet,jet],[nu],[nu]]])" : ["T7ChiSlep",["ATLAS_CONF_2013_007"]]  
    }    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)     
    

# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights. Check if listed histograms
# exist and store binning information in masscomp
    for Analy in ListOfAnalyses:
        Analy.GenerateElements()
        Analy.GetPlots()


    return ListOfAnalyses        
    
    
    

