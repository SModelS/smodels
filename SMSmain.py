#!/usr/bin/env python

import SMSglobals, SMSanalyses, sys, SMSmethods, SMSxsec, SMSgetlimit
from prettytable import PrettyTable
#import SMSxsec

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Initialize global variables:
SMSglobals.initglob()
#Creat analyses list:
SMSanalyses.load()
#Generate events and compute cross-sections:
nevts = 5
slhafile = "AndreSLHA/andrePT4.slha"
Wv = SMSxsec.pytinit(nevts,slhafile)
W = Wv["Wdic"]
lhefile = Wv["lhefile"]
LHEfile = open(lhefile,"r")


SMSTopList = []     
#Read events and get topologies (fills SMSTopList)
for iev in range(nevts):

##Read event    
    PList = SMSmethods.getNextEvent(LHEfile) 
##Get mother PDGs:
    momPDG = tuple(SMSmethods.getMom(PList))
#Get event weight list:
    weight = {}
    
    for k in W.keys(): 
        if W[k].has_key(momPDG):
            weight.update({k : W[k][momPDG]})
        else:
            print "Error getting weight"
            sys.exit()    

#Get event topology    
    SMSTopListEv = SMSmethods.getEventTop(PList, weight)
#Add event topology to topology list:  
    for TopEv in SMSTopListEv:  
        SMSmethods.AddToList(TopEv,SMSTopList)
    SMSglobals.evcount +=1

element_table = PrettyTable(["#Vertices B[0]","#Vertices B[1]", "#Insertions B[0]","#Insertions B[1]", "#Elements", "Sum of weights"])          

eltot = 0
#Print Results:
print "Number of Global topologies = ",len(SMSTopList)
for i in range(len(SMSTopList)):
    print "Topology ",i    
    sumw = weight
    for w in sumw.keys():
        sumw[w] = 0.
#    print "Number of vertices = ",SMSTopList[i].B[0].vertnumb,SMSTopList[i].B[1].vertnumb
#    print "Vertice insertions = ",SMSTopList[i].B[0].vertparts,SMSTopList[i].B[1].vertparts
#    print "Number of Elements = ",len(SMSTopList[i].B[0].ElList)
    for j in range(len(SMSTopList[i].WeightList)):
        for w in SMSTopList[i].WeightList[j].keys():
            sumw[w] += SMSTopList[i].WeightList[j][w]            
#    print "Sum of weights=", sumw,"\n"
    element_table.add_row([SMSTopList[i].B[0].vertnumb,SMSTopList[i].B[1].vertnumb,SMSTopList[i].B[0].vertparts,SMSTopList[i].B[1].vertparts,len(SMSTopList[i].B[0].ElList),sumw])
    eltot += len(SMSTopList[i].B[0].ElList)

 
            
#Print element list for Topology[i]:    
    if i < 0:
        for j in range(len(SMSTopList[i].B[0].ElList)):
            print SMSTopList[i].B[0].ElList[j].masses,SMSTopList[i].B[1].ElList[j].masses,SMSTopList[i].B[0].ElList[j].particles,SMSTopList[i].B[1].ElList[j].particles,SMSTopList[i].WeightList[j]
print(element_table)

print "Total Number of Elements = ",eltot   
print '\n \n \n'


#Add event topologies to analyses:
for Analysis in SMSglobals.ListOfAnalyses:
    SMSmethods.AddToAnalysis(SMSTopList,Analysis)

#Print analyses output:
#for Anal in SMSglobals.ListOfAnalyses:
#    print Anal.label
#    for i in range(len(Anal.Top.B[0].ElList)):
#        print Anal.Top.B[0].ElList[i].particles,Anal.Top.B[1].ElList[i].particles
#        for j in range(len(Anal.Top.B[0].ElList[i].masses)):
#            if Anal.Top.WeightList[i][j]["7 TeV"] >= 0.:
#                print Anal.Top.B[0].ElList[i].masses[j],Anal.Top.B[1].ElList[i].masses[j],Anal.Top.WeightList[i][j]
#    print '\n'

print '\n \n \n'

y = PrettyTable(["#Vertices B[0]","#Vertices B[1]", "#Insertions B[0]","#Insertions B[1]", "#Elements", "Sum of weights"])          

#Compute theoretical predictions to analyses results:
for Analysis in SMSglobals.ListOfAnalyses:
    print "---------------Analysis Label = ",Analysis.label
    for i in range(len(Analysis.results)):
        const = Analysis.results.items()[i][0]
        cond = Analysis.results.items()[i][1]
        constRes = SMSmethods.Aeval(Analysis,const)
        condRes = SMSmethods.Aeval(Analysis,cond)
        

        try:
            plot = Analysis.plots[const]
        except KeyError:    
            plot = []
                
        print "Constraint: ",const,"  Conditions: ",cond
        for j in range(len(constRes)):
            massarray = []
            massarray.append([x for x in constRes[j][0][0]])
            massarray.append([x for x in constRes[j][0][1]])
            sigmalimit = SMSgetlimit.GetPlotLimit(massarray,plot,Analysis)
            print "Results for mass array ",constRes[j][0]," :"
            print "    Theorical value = ",constRes[j][1],"Conditions = ",condRes[j][1]
            print "    Experimental Limits: ",sigmalimit
        print '\n'
    print "-----------------------"
    print '\n \n'    
    sys.exit()  # Just print first analysis (avoid too much output)       



sys.exit()

