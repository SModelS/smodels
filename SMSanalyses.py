def load():
    from SMSmethods import EAnalysis
    from SMSglobals import ListOfAnalyses
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
    "[[[t,t]],[[t,t]]]" : ["T1tttt",["alphaT","RA48TeV","RA2b8TeV","MultiLepton8TeV","alphaT8TeV","ATLAS_CONF_2013_007"]],
    "[[[b,b]],[[b,b]]]" : ["T1bbbb",["alphaT","RA2b8TeV","alphaT8TeV"]],
    "[[[t,b]],[[t,b]]]" : ["T1tbtb",["ATLAS_CONF_2013_007"]],
    "[[[t+,t-]],[[t+,t-]]]" : ["T1tttt",["ATLAS_CONF_2012_015"]]
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
    "2.*([[[L],[L]],[[L],[nu]]] + [[[L],[L]],[[nu],[L]]])" : "[[[L],[L]],[[L],[nu]]] ~ [[[L],[L]],[[nu],[L]]], [[[L],[L]],[[L],[nu]]] > 2.7*[[[ta],[ta]],[[L],[nu]]], [[[L],[L]],[[L],[nu]]] > 2.7*[[[L],[L]],[[ta],[nu]]], [[[L],[L]],[[nu],[L]]] > 2.7*[[[ta],[ta]],[[nu],[L]]], [[[L],[L]],[[nu],[L]]] > 2.7*[[[L],[L]],[[nu],[ta]]]",
#TChiChiSlep (tau-enriched chargino decay):    
    "[[[L],[L]],[[nu],[ta]]]" : "[[[L],[L]],[[nu],[ta]]] > 2.7*[[[ta],[ta]],[[nu],[ta]]]",
#TChiChiSlep (tau-dominated  chargino and neutralino decays):
    "[[[ta],[ta]],[[nu],[ta]]]" : "None",
#Chargino-Chargino:
    "[[[L+],[nu]],[[L-],[nu]]] + [[[L+],[nu]],[[nu],[L-]]] + [[[L-],[nu]],[[nu],[L+]]] + [[[nu],[L+]],[[nu],[L-]]]" : "[[[L+],[nu]],[[L-],[nu]]] ~ [[[L+],[nu]],[[nu],[L-]]], [[[L+],[nu]],[[nu],[L-]]] ~ [[[L-],[nu]],[[nu],[L+]]], [[[L-],[nu]],[[nu],[L+]]] ~ [[[nu],[L+]],[[nu],[L-]]], [[[L+],[nu]],[[L-],[nu]]] > 2.7*[[[ta+],[nu]],[[L-],[nu]]], [[[L+],[nu]],[[L-],[nu]]] > 2.7*[[[L+],[nu]],[[ta-],[nu]]], [[[L+],[nu]],[[nu],[L-]]] > 2.7*[[[ta+],[nu]],[[nu],[L-]]], [[[L+],[nu]],[[nu],[L-]]] > 2.7*[[[L+],[nu]],[[nu],[ta-]]], [[[L-],[nu]],[[nu],[L+]]] > 2.7*[[[ta-],[nu]],[[nu],[L+]]], [[[L-],[nu]],[[nu],[L+]]] > 2.7*[[[L-],[nu]],[[nu],[ta+]]], [[[nu],[L+]],[[nu],[L-]]] > 2.7*[[[nu],[ta+]],[[nu],[L-]]], [[[nu],[L+]],[[nu],[L-]]] > 2.7*[[[nu],[L+]],[[nu],[ta-]]]",
#T6bbww (stop pair production)
    "[[[b],[W]],[[b],[W]]]" : "None",
#T6tttt
    "[[[t],[t]],[[t],[t]]]" : "None",
#T6ttww
    "[[[t],[W]],[[t],[W]]]" : "None",
#T6ttzz    
    "[[[Z],[t]],[[Z],[t]]]" : "None",
#TChiChiSlepSlep
    "[[[l],[l]],[[l],[l]]]" : "[[[e],[l]],[[l],[l]]] < 0.63*[[[l],[l]],[[l],[l]]], [[[l],[e]],[[l],[l]]] < 0.63*[[[l],[l]],[[l],[l]]]",
#TChiChipmStauL    
    "2.*([[[ta+],[ta-]],[[nu],[ta]]] + [[[ta+],[ta-]],[[ta],[nu]]] + [[[ta-],[ta+]],[[nu],[ta]]] + [[[ta-],[ta+]],[[ta],[nu]]])" : "[[[ta+],[ta-]],[[nu],[ta]]] ~ [[[ta+],[ta-]],[[ta],[nu]]], [[[ta+],[ta-]],[[ta],[nu]]] ~ [[[ta-],[ta+]],[[nu],[ta]]], [[[ta-],[ta+]],[[nu],[ta]]] ~  [[[ta-],[ta+]],[[ta],[nu]]]",
#TChipChimStauSnu
    "[[[ta+],[nu]],[[ta-],[nu]]] + [[[ta+],[nu]],[[nu],[ta-]]] + [[[ta-],[nu]],[[nu],[ta+]]] + [[[nu],[ta+]],[[nu],[ta-]]]" : "[[[ta+],[nu]],[[ta-],[nu]]] ~ [[[ta+],[nu]],[[nu],[ta-]]], [[[ta+],[nu]],[[nu],[ta-]]] ~ [[[ta-],[nu]],[[nu],[ta+]]], [[[ta-],[nu]],[[nu],[ta+]]] ~ [[[nu],[ta+]],[[nu],[ta-]]]"    
    }
    Analysis.plots  = {   
    "2.*([[[L],[L]],[[L],[nu]]] + [[[L],[L]],[[nu],[L]]])" : ["TChiChipmSlepL",["Weakinos8TeV","ATLAS_CONF_2013_035"]],
    "[[[L],[L]],[[nu],[ta]]]" : ["TChiChipmSlepStau",["Weakinos8TeV"]],
    "[[[ta],[ta]],[[nu],[ta]]]" : ["TChiChipmStauStau",["Weakinos8TeV"]],
    "[[[L+],[nu]],[[L-],[nu]]] + [[[L+],[nu]],[[nu],[L-]]] + [[[L-],[nu]],[[nu],[L+]]] + [[[nu],[L+]],[[nu],[L-]]]" : ["TChipChimSlepSnu",["Weakinos8TeV"]],
    "[[[b],[W]],[[b],[W]]]" : ["T6bbww",["LeptonicStop8TeV","ATLAS_CONF_2013_037","ATLAS_CONF_2012_166"]],
    "[[[t],[t]],[[t],[t]]]" : ["T6tttt",["ATLAS_CONF_2013_007"]],
    "[[[t],[W]],[[t],[W]]]" : ["T6ttww",["ATLAS_CONF_2013_007"]],
    "[[[Z],[t]],[[Z],[t]]]" : ["T6ttZZ",["ATLAS_CONF_2013_025"]],
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
    "[[[jet],[nu],[L]],[[jet],[nu],[L]]] + [[[jet],[nu],[L]],[[jet],[L],[nu]]] + [[[jet],[nu],[L]],[[jet],[L],[L]]] + [[[jet],[nu],[L]],[[jet],[nu],[nu]]] + [[[jet],[L],[nu]],[[jet],[L],[nu]]] + [[[jet],[L],[nu]],[[jet],[L],[L]]] + [[[jet],[L],[nu]],[[jet],[nu],[nu]]] + [[[jet],[L],[L]],[[jet],[L],[L]]] + [[[jet],[L],[L]],[[jet],[nu],[nu]]] + [[[jet],[nu],[nu]],[[jet],[nu],[nu]]]" : "None?"        
    }
    Analysis.plots  = {   
    "[[[jet],[nu],[L]],[[jet],[nu],[L]]] + [[[jet],[nu],[L]],[[jet],[L],[nu]]] + [[[jet],[nu],[L]],[[jet],[L],[L]]] + [[[jet],[nu],[L]],[[jet],[nu],[nu]]] + [[[jet],[L],[nu]],[[jet],[L],[nu]]] + [[[jet],[L],[nu]],[[jet],[L],[L]]] + [[[jet],[L],[nu]],[[jet],[nu],[nu]]] + [[[jet],[L],[L]],[[jet],[L],[L]]] + [[[jet],[L],[L]],[[jet],[nu],[nu]]] + [[[jet],[nu],[nu]],[[jet],[nu],[nu]]]" : ["T8ChiSlep",["ATLAS_CONF_2013_007"]]  
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
    "[[[jet,jet],[nu],[L]],[[jet,jet],[nu],[L]]] + [[[jet,jet],[nu],[L]],[[jet,jet],[L],[nu]]] + [[[jet,jet],[nu],[L]],[[jet,jet],[L],[L]]] + [[[jet,jet],[nu],[L]],[[jet,jet],[nu],[nu]]] + [[[jet,jet],[L],[nu]],[[jet,jet],[L],[nu]]] + [[[jet,jet],[L],[nu]],[[jet,jet],[L],[L]]] + [[[jet,jet],[L],[nu]],[[jet,jet],[nu],[nu]]] + [[[jet,jet],[L],[L]],[[jet,jet],[L],[L]]] + [[[jet,jet],[L],[L]],[[jet,jet],[nu],[nu]]] + [[[jet,jet],[nu],[nu]],[[jet,jet],[nu],[nu]]]" : "None?"        
    }
    Analysis.plots  = {   
    "[[[jet,jet],[nu],[L]],[[jet,jet],[nu],[L]]] + [[[jet,jet],[nu],[L]],[[jet,jet],[L],[nu]]] + [[[jet,jet],[nu],[L]],[[jet,jet],[L],[L]]] + [[[jet,jet],[nu],[L]],[[jet,jet],[nu],[nu]]] + [[[jet,jet],[L],[nu]],[[jet,jet],[L],[nu]]] + [[[jet,jet],[L],[nu]],[[jet,jet],[L],[L]]] + [[[jet,jet],[L],[nu]],[[jet,jet],[nu],[nu]]] + [[[jet,jet],[L],[L]],[[jet,jet],[L],[L]]] + [[[jet,jet],[L],[L]],[[jet,jet],[nu],[nu]]] + [[[jet,jet],[nu],[nu]],[[jet,jet],[nu],[nu]]]" : ["T7ChiSlep",["ATLAS_CONF_2013_007"]]  
    }    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)     
    
    
    
    
  
    

# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights. Check if listed histograms
# exist and store binning information in masscomp
    for Analy in ListOfAnalyses:
        Analy.GenerateElements()
        Analy.GetPlots()
    
    
    
    

