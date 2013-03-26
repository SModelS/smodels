import SMSglobals, SMSanalyses, sys, SMSmethods

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Initialize global variables:
SMSglobals.initglob()
#Creat analyses list:
SMSanalyses.load()

#test
nevts = 10000
Sqrts = 7
pytres = {"xsecfb" : 1.}
slhafile = "/home/lessa/SMS_pythia/SLHApts/andrePT4.slha"
#pytres = SMSmethods.runpythia(slhafile,nevts,Sqrts)

filename = open("./data/fort.68","r")

SMSTopList = []

#Read events and get topologies (fills SMSTopList)
for iev in range(nevts):
    weight = [pytres["xsecfb"]/float(nevts),pytres["xsecfb"]/float(nevts)]
    
    PList = SMSmethods.getNextEvent(filename)  
    SMSTop = SMSmethods.GTop()
    SMSTop = SMSmethods.getEventTop(PList, weight)

#Collect topologies:    
    SMSmethods.AddToList(SMSTop,SMSTopList)
    SMSglobals.evcount +=1
        


eltot = 0
#Print Results:
print "Number of Global topologies = ",len(SMSTopList)
for i in range(len(SMSTopList)):
    print "Topology ",i
    sumw = [0.]*len(weight)
    print "Number of vertices = ",SMSTopList[i].B[0].vertnumb,SMSTopList[i].B[1].vertnumb
    print "Vertice insertions = ",SMSTopList[i].B[0].vertins,SMSTopList[i].B[1].vertins
    print "Number of Elements = ",len(SMSTopList[i].B[0].ElList)
    for j in range(len(SMSTopList[i].WeightList)):
        for iw in range(len(SMSTopList[i].WeightList[j])):
            sumw[iw] += SMSTopList[i].WeightList[j][iw]            
    print "Sum of weights=", sumw,"\n"
    eltot += len(SMSTopList[i].B[0].ElList)
 
            
#Print element list for Topology[i]:            
    if i < 0:
        for j in range(len(SMSTopList[i].B[0].ElList)):
            print SMSTopList[i].B[0].ElList[j].masses,SMSTopList[i].B[1].ElList[j].masses,SMSTopList[i].B[0].ElList[j].particles,SMSTopList[i].B[1].ElList[j].particles,SMSTopList[i].WeightList[j]

            
            
            

print "Total Number of Elements = ",eltot   
print '\n \n \n'


#sys.exit()
#Add Elements to analyses:
SMSmethods.AddToAnalyses(SMSTopList)


for itop in range(len(SMSglobals.AnalysesRes)):
#    print "Topology ",itop
#    print "Number of vertices = ",SMSglobals.AnalysesRes[itop].Top.B[0].vertnumb,SMSglobals.AnalysesRes[itop].Top.B[1].vertnumb
#    print "Vertice insertions = ",SMSglobals.AnalysesRes[itop].Top.B[0].vertins,SMSglobals.AnalysesRes[itop].Top.B[1].vertins


    SMSglobals.AnalysesRes[itop].CreateAll()

#    for i in range(len(SMSglobals.AnalysesRes[itop].StrList)):        
#        print SMSglobals.AnalysesRes[itop].StrList[i]
#        for j in range(len(SMSglobals.AnalysesRes[itop].ElList[i].MassList)):
#            if sum(SMSglobals.AnalysesRes[itop].ElList[i].WeightList[j]) >= 0.:
#                print SMSglobals.AnalysesRes[itop].ElList[i].MassList[j],SMSglobals.AnalysesRes[itop].ElList[i].WeightList[j]
#        print '\n'
#sys.exit()

print '\n \n \n'

for Analysis in SMSglobals.ListOfAnalyses:
    for i in range(len(Analysis.constraints)):
        const = Analysis.constraints.items()[i][0]
        cond = Analysis.constraints.items()[i][1]
        constRes = SMSmethods.Aeval(Analysis,const)
        condRes = SMSmethods.Aeval(Analysis,cond)
        

        try:
            plot = Analysis.plots[const]
        except KeyError:    
            plot = "bla"
                
        print const,":     ",cond,":"
        for j in range(len(constRes)):            
            if sum(constRes[j][1][0]) != 0.:
                sigmalimit = SMSmethods.GetPlotLimit(constRes[j][0],plot)
                print constRes[j][0],constRes[j][1],condRes[j][1]
        print '\n'        



#If cuts on final state particles have been included, check for number of
# compressed topologies:
#if SMSglobals.min_pt + SMSglobals.min_lep_pt + SMSglobals.min_tau_pt + SMSglobals.min_jet_pt + SMSglobals.min_b_pt != 0:
#    compfr = float(sum(SMSglobals.ncomp))/float(2*nevts)
#    print "Fraction of compressed topologies = ",100.*compfr,"%"
#    if compfr > 0.25:
#        print "***Compressed Spectrum Detected*** \n"
#    
#    
#    
sys.exit()

