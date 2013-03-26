import SMSglobals, SMSanalyses, sys, SMSmethods

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Initialize global variables:
SMSglobals.initglob()
#Creat analyses list:
SMSanalyses.load()

     


nevts = 10000
Sqrts = 7
pytres = {"xsecfb" : 1.}
slhafile = "AndreSLHA/andrePT4.slha"
pytres = SMSmethods.runpythia(slhafile,nevts,Sqrts)

filename = open("./data/fort.68","r")

SMSTopList = []     
#Read events and get topologies (fills SMSTopList)
for iev in range(nevts):
    weight = {"7 TeV" : pytres["xsecfb"]/float(nevts), "8 TeV" : pytres["xsecfb"]/float(nevts)}

#Read event    
    PList = SMSmethods.getNextEvent(filename)  
#Get event topology    
    SMSTopListEv = SMSmethods.getEventTop(PList, weight)
#Add event topology to topology list:  
    for TopEv in SMSTopListEv:  
        SMSmethods.AddToList(TopEv,SMSTopList)
    SMSglobals.evcount +=1


eltot = 0
#Print Results:
print "Number of Global topologies = ",len(SMSTopList)
for i in range(len(SMSTopList)):
    print "Topology ",i    
    sumw = weight
    for w in sumw.keys():
        sumw[w] = 0.
    print "Number of vertices = ",SMSTopList[i].B[0].vertnumb,SMSTopList[i].B[1].vertnumb
    print "Vertice insertions = ",SMSTopList[i].B[0].vertparts,SMSTopList[i].B[1].vertparts
    print "Number of Elements = ",len(SMSTopList[i].B[0].ElList)
    for j in range(len(SMSTopList[i].WeightList)):
        for w in SMSTopList[i].WeightList[j].keys():
            sumw[w] += SMSTopList[i].WeightList[j][w]            
    print "Sum of weights=", sumw,"\n"
    eltot += len(SMSTopList[i].B[0].ElList)
 
            
#Print element list for Topology[i]:    
    if i < 0:
        for j in range(len(SMSTopList[i].B[0].ElList)):
            print SMSTopList[i].B[0].ElList[j].masses,SMSTopList[i].B[1].ElList[j].masses,SMSTopList[i].B[0].ElList[j].particles,SMSTopList[i].B[1].ElList[j].particles,SMSTopList[i].WeightList[j]


print "Total Number of Elements = ",eltot   
print '\n \n \n'


#Add event topologies to analyses:
for Analysis in SMSglobals.ListOfAnalyses:
    SMSmethods.AddToAnalysis(SMSTopList,Analysis)

#Print analyses output:
for Anal in SMSglobals.ListOfAnalyses:
    print Anal.label
    for i in range(len(Anal.Top.B[0].ElList)):
        print Anal.Top.B[0].ElList[i].particles,Anal.Top.B[1].ElList[i].particles
        for j in range(len(Anal.Top.B[0].ElList[i].masses)):
            if Anal.Top.WeightList[i][j]["7 TeV"] >= 0.:
                print Anal.Top.B[0].ElList[i].masses[j],Anal.Top.B[1].ElList[i].masses[j],Anal.Top.WeightList[i][j]
    print '\n'

print '\n \n \n'

#Compute theoretical predictions to analyses results:
for Analysis in SMSglobals.ListOfAnalyses:
    print Analysis.label
    for i in range(len(Analysis.results)):
        const = Analysis.results.items()[i][0]
        cond = Analysis.results.items()[i][1]
        constRes = SMSmethods.Aeval(Analysis,const)
        condRes = SMSmethods.Aeval(Analysis,cond)
        

        try:
            plot = Analysis.plots[const]
        except KeyError:    
            plot = "bla"
                
        print const,":     ",cond,":"
        for j in range(len(constRes)):            
            sigmalimit = SMSmethods.GetPlotLimit(constRes[j][0],plot)
            print constRes[j][0],constRes[j][1],condRes[j][1]
        print '\n'        



sys.exit()

