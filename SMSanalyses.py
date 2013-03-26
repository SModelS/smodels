def load():
    from SMSmethods import EAnalysis
    from SMSglobals import ListOfAnalyses
    #List analyses and classify in topologies and elements
        

    
    #---LeptonsMT(1209.6620)----
    #common info:
    lum = 4.98
    sqrts = 7.
    histname="/home/lessa/SMS_pythia/2012/LeptonsMT/hist.root"
    
    
    
    #Slepton Topologies:
    Analysis = EAnalysis()
    Analysis.lum = lum
    Analysis.sqrts = sqrts
    Analysis.histname = histname     
    Analysis.label = "LeptonsMT - C1N2_slepton"
    Analysis.Top.B[0].vertnumb = 3
    Analysis.Top.B[0].vertparts = [1,1,0]
    Analysis.Top.B[1].vertnumb = 3
    Analysis.Top.B[1].vertparts = [1,1,0]
    #Constraints
    Analysis.results = {
#Left-handed sleptons
"[[L,L],[L,nu]] + [[L,L],[nu,L]]" : "[[L,L],[L,nu]] ~ [[L,L],[nu,L]], [[L,L],[L,nu]] ~ 9.*[[ta,ta],[nu,ta]] ",
#Right-handed sleptons
"[[L,L],[nu,ta]]" : "[[L,L],[nu,ta]] ~ (3./2.)*[[l,l],[nu,ta]]"
}
    Analysis.plotpath = "/home/lessa/SMS_pythia/2012/LeptonsMT/"
    Analysis.plots  = {
"[[L,L],[L,nu]] + [[L,L],[nu,L]]" : "plot1",
"[[L,L],[nu,ta]]" : "plot2"
}
    ListOfAnalyses.append(Analysis)

#    
#    
#    #VV Topologies:
    Analysis = EAnalysis()
    Analysis.lum = lum
    Analysis.sqrts = sqrts
    Analysis.histname = histname  
    Analysis.label = "LeptonsMT - C1N2_VV"
    Analysis.Top.B[0].vertnumb = 2
    Analysis.Top.B[0].vertparts = [1,0]
    Analysis.Top.B[1].vertnumb = 2
    Analysis.Top.B[1].vertparts = [1,0]
    #Constraints    
    Analysis.results = {
#C1N2_WZ
"[[W],[Z]]" : "None",
#C1N2_ZZ
"[[Z],[Z]]" : "None"
}    
    ListOfAnalyses.append(Analysis)
    
    
    
    
    
    #---LeptonsMT8(SUS-12-022)----
    #general info:
    lum = 9.2
    sqrts = 8.
    histname="/home/lessa/SMS_pythia/2012/LeptonsMT8/hist.root"
    
    
    
    #Slepton Topologies:
    Analysis = EAnalysis()
    Analysis.lum = lum
    Analysis.sqrts = sqrts
    Analysis.histname = histname  
    Analysis.label = "LeptonsMT8 - C1N2_slepton + C1C1_slepton"
    Analysis.Top.B[0].vertnumb = 3
    Analysis.Top.B[0].vertparts = [1,1,0]
    Analysis.Top.B[1].vertnumb = 3
    Analysis.Top.B[1].vertparts = [1,1,0]
    #Contraints
    Analysis.results = {
#Left-handed sleptons
#"[[L,L],[L,nu]] + [[L,L],[nu,L]]" : "[[L,L],[L,nu]] ~ [[L,L],[nu,L]]",
"[[L,nu],[L,L]] + [[nu,L],[L,L]]" : "[[L,nu],[L,L]] ~ [[nu,L],[L,L]]",
#Right-handed sleptons
"[[L,L],[nu,ta]]" : "None",
#Full right-handed sleptons
"[[ta,ta],[nu,ta]]" : "None",
#C1C1 production
"[[L,nu],[nu,L]] + [[L,nu],[L,nu]] + [[nu,L],[nu,L]]" : "[[L,nu],[nu,L]] ~ 2*[[L,nu],[L,nu]], [[L,nu],[nu,L]] ~ 2*[[nu,L],[nu,L]]"
}
    ListOfAnalyses.append(Analysis)
    
    
    
    #VV Topologies:
    Analysis = EAnalysis()
    Analysis.lum = lum
    Analysis.sqrts = sqrts
    Analysis.histname = histname    
    Analysis.label = "LeptonsMT8 - C1N2_VV + SlepSlep"
    Analysis.Top.B[0].vertnumb = 2
    Analysis.Top.B[0].vertparts = [1,0]
    Analysis.Top.B[1].vertnumb = 2
    Analysis.Top.B[1].vertparts = [1,0]
    #Constraints    
    Analysis.results = {
#C1N2_WZ
"[[W],[Z]]" : "None",
#C1N2_ZZ
"[[Z],[Z]]" : "None",
#SlepSlep
"[[l+],[l-]]" : "None"
}
    ListOfAnalyses.append(Analysis)
    
    
    #Test Topology:
    Analysis = EAnalysis()
    Analysis.lum = lum
    Analysis.sqrts = sqrts
    Analysis.histname = histname    
    Analysis.label = "Test"
    Analysis.Top.B[0].vertnumb = 3
    Analysis.Top.B[0].vertparts = [1,1,0]
    Analysis.Top.B[1].vertnumb = 2
    Analysis.Top.B[1].vertparts = [1,0]
    #Constraints    
    Analysis.results = {
"[[W,nu],[ta]]" : "None"
}
    ListOfAnalyses.append(Analysis)

# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights
    for Analy in ListOfAnalyses:
        Analy.GenerateElements() 
    
    
    
    

