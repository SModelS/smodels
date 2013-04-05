#!/usr/bin/env python

def load():
    """ This method creates the analysis objects from the info given in the SMS database.
        These EAnalysis objects are registered at SMSglobals.ListOfAnalyses """
    from SMSmethods import EAnalysis
    from SMSglobals import ListOfAnalyses
    import SMSResults

    analysis="RAZOR"
    topo="T1tttt"

    print SMSResults.getLumi( analysis )

    #List analyses and classify in topologies and elements
        
    #---Analyses example----
    #common info:
    sqrts = 7.    #Each Analyses "group" (EAnalysis object) should only combine results with a common sqrts (for simplicity)
    masscomp = 0.2  #Percentual difference to define equal masses (to be changed/improved later)
    #Global Topology:
    Analysis = EAnalysis()
    Analysis.sqrts = sqrts
    Analysis.masscomp = masscomp
    Analysis.plotpath = "2012"  #Run label
    Analysis.label = "T2-type" #it can be whatever you want
    Analysis.Top.B[0].vertnumb = 2    #Number of vertices of 1st branch
    Analysis.Top.B[0].vertparts = [1,0] #Number of particle insertions of 1st branch
    Analysis.Top.B[1].vertnumb = 2    #Number of vertices of 2nd branch
    Analysis.Top.B[1].vertparts = [1,0] #Number of particle insertions of 2nd branch    
    
    #Now list the elements constrained by CMS analyses which match this global topology and the corresponding conditions (assumptions):
    Analysis.results = {
#T2:
    "[[jet],[jet]]" : "None",
#T2bb:
    "[[b],[b]]" : "None",
#T2tt:
    "[[t],[t]]" : "None",
#TChizz:    (nothing forbids us to include the elements below in a different EAnalysis object, if we want to keep leptonic and hadronic analysis separated)
    "[[Z],[Z]]" : "None",
#TChiww:
    "[[W],[W]]" : "None",
#TChiwz:
    "[[W],[Z]]" : "None",
#TSlepSlep:
    "[[l],[l]]" : "None"        
}
#Now list the information necessary to find the experimental results for each element (CMS label + plus list of analyses)
    Analysis.plots  = {
    "[[jet],[jet]]" : ["T2",["RAZOR","RA2","alphaTb","alphaT"]],
    "[[b],[b]]" : ["T2bb",["RAZORb","alphaTb","alphaT"]],
    "[[t],[t]]" : ["T2tt",["RAZOR","HadronicStop","alphaTb","alphaT","RAZORb","RAZORjets"]],
    "[[Z],[Z]]" : ["TChizz",["MULTILEPTON"]],
    "[[W],[W]]" : ["TChiww",[]],
    "[[W],[Z]]" : ["TChiwz",["MULTILEPTON","LeptonsMT"]],
    "[[l],[l]]" : ["TSlepSlep",[]]}
    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)

# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights
    for Analy in ListOfAnalyses:
        Analy.GenerateElements() 
    
if __name__ == "__main__":
  from SMSglobals import ListOfAnalyses
  load()
