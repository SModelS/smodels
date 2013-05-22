""" Upon 'loading', this unit creates a list of EAnalyses """

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



# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights. Check if listed histograms
# exist and store binning information in masscomp
    for Analy in ListOfAnalyses:
        Analy.GenerateElements()
        Analy.GetPlots()


    return ListOfAnalyses        
    
    
    

