def load():
    from SMSmethods import EAnalysis
    from SMSglobals import ListOfAnalyses
    #List analyses and classify in topologies and elements
        
     #Global Topology:
    Analysis = EAnalysis()
    Analysis.label = "T1-simple" #it can be whatever you want
    Analysis.Top.B[0].vertnumb = 2    #Number of vertices of 1st branch
    Analysis.Top.B[0].vertparts = [2,0] #Number of particle insertions of 1st branch
    Analysis.Top.B[1].vertnumb = 2    #Number of vertices of 2nd branch
    Analysis.Top.B[1].vertparts = [2,0] #Number of particle insertions of 2nd branch    
    Analysis.results = { "[[[jet,jet]],[[jet,jet]]]" : "None", } # T1
    Analysis.plots  = { "[[[jet,jet]],[[jet,jet]]]" : ["T1",["alphaT"]] }    
#Add analysis to list of analyses:
    ListOfAnalyses.append(Analysis)
    
# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights. Check if listed histograms
# exist and store binning information in masscomp
    for Analy in ListOfAnalyses:
        Analy.GenerateElements()
        Analy.GetPlots()
    
    
    
    

