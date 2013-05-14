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
    "[[[b,b]],[[b,b]]]" : "None"
    }
    Analysis.plots  = {
    "[[[jet,jet]],[[jet,jet]]]" : ["T1",["alphaT"]],
    "[[[t,t]],[[t,t]]]" : ["T1tttt",["alphaT","RA48TeV","RA2b8TeV"]],
    "[[[b,b]],[[b,b]]]" : ["T1bbbb",["alphaT","RA2b8TeV"]]
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
    "[[[jet]],[[jet]]]" : ["T2",["alphaT"]],
    "[[[t]],[[t]]]" : ["T2tt",["alphaT","LeptonicStop8TeV"]],
    "[[[b]],[[b]]]" : ["T2bb",["alphaT"]],
    "[[[W]],[[Z]]]" : ["TChiwz",["Weakinos8TeV"]],
#    "[[Z],[Z]]" : ["TChiZZ",["Weakinos8TeV"]],
    "[[[l+]],[[l-]]]" : ["TSlepSlep",["Weakinos8TeV"]]
    }    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)
    


    
#         #Global Topology:
    Analysis = EAnalysis()
    Analysis.run = "8TeV"  #Run label
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
    "[[[b],[W]],[[b],[W]]]" : "None"    
    }
    Analysis.plots  = {   
    "2.*([[[L],[L]],[[L],[nu]]] + [[[L],[L]],[[nu],[L]]])" : ["TChiChipmSlepL",["Weakinos8TeV"]],
    "[[[L],[L]],[[nu],[ta]]]" : ["TChiChipmSlepStau",["Weakinos8TeV"]],
    "[[[ta],[ta]],[[nu],[ta]]]" : ["TChiChipmStauStau",["Weakinos8TeV"]],
    "[[[L+],[nu]],[[L-],[nu]]] + [[[L+],[nu]],[[nu],[L-]]] + [[[L-],[nu]],[[nu],[L+]]] + [[[nu],[L+]],[[nu],[L-]]]" : ["TChipChimSlepSnu",["Weakinos8TeV"]],
    "[[[b],[W]],[[b],[W]]]" : ["T6bbww",["LeptonicStop8TeV"]]
    }    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)
    
    
    
    
    
   
    
    
    
    
    
    
    
    

# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights. Check if listed histograms
# exist and store binning information in masscomp
    for Analy in ListOfAnalyses:
        Analy.GenerateElements()
        Analy.GetPlots()
    
    
    
    

