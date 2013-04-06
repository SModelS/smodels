#!/usr/bin/env python

import SMSglobals, SMSanalyses, sys, SMSmethods, SMSxsec, SMSgetlimit, copy
from prettytable import PrettyTable
#import SMSxsec

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Initialize global variables:
SMSglobals.initglob()
#Creat analyses list:
SMSanalyses.load()

#Generate events and compute cross-sections:
nevts = 10000
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

EvTop_table = PrettyTable(["Topology","#Vertices B[0]","#Vertices B[1]", "#Insertions B[0]","#Insertions B[1]", "#Elements", "Sum of weights"])
EvElement_table = PrettyTable(["Topology","Element","Particles B[0]","Particles B[1]", "Masses B[0]","Masses B[1]","Element Weight"])

   

eltot = 0
#Print Results:
for i in range(len(SMSTopList)):
    sumw = {}
    sumw.update(weight)
    for w in sumw.keys():
        sumw[w] = 0.
    for j in range(len(SMSTopList[i].WeightList)):
        for w in SMSTopList[i].WeightList[j].keys():
            sumw[w] += SMSTopList[i].WeightList[j][w]            
    EvTop_table.add_row([i,SMSTopList[i].B[0].vertnumb,SMSTopList[i].B[1].vertnumb,SMSTopList[i].B[0].vertparts,SMSTopList[i].B[1].vertparts,len(SMSTopList[i].B[0].ElList),sumw])
    eltot += len(SMSTopList[i].B[0].ElList)

 
            
#Print element list for Topology[i]:    
    if i >= 0:           
        for j in range(len(SMSTopList[i].B[0].ElList)):
            EvElement_table.add_row([i,j,SMSTopList[i].B[0].ElList[j].particles,SMSTopList[i].B[1].ElList[j].particles,SMSTopList[i].B[0].ElList[j].masses,SMSTopList[i].B[1].ElList[j].masses,SMSTopList[i].WeightList[j]])
        EvElement_table.add_row(["---","---","---","---","---","---","---"])    
            
            
print "Number of Global topologies = ",len(SMSTopList)          
print(EvTop_table)
print "Total Number of Elements = ",eltot
#print(EvElement_table)


print '\n \n \n'


#Add event topologies to analyses:
for Analysis in SMSglobals.ListOfAnalyses:
    SMSmethods.AddToAnalysis(SMSTopList,Analysis)

#Print analyses output:
AnElement_table = PrettyTable(["Analyses","Element","Masses","Element Weight"])  

for Ana in SMSglobals.ListOfAnalyses:
    for i in range(len(Ana.Top.B[0].ElList)):
        for j in range(len(Ana.Top.B[0].ElList[i].masses)):
            AnElement_table.add_row([Ana.label,[Ana.Top.B[0].ElList[i].particles,Ana.Top.B[1].ElList[i].particles],[Ana.Top.B[0].ElList[i].masses[j],Ana.Top.B[1].ElList[i].masses[j]],Ana.Top.WeightList[i][j]])
    AnElement_table.add_row(["---","---","---","---"])  
        

print(AnElement_table)

print '\n \n \n'

         
#Compute theoretical predictions to analyses results:
for Analysis in SMSglobals.ListOfAnalyses:
    print "---------------Analysis Label = ",Analysis.label
    Results_table = PrettyTable(["Result","Conditions","Mass","Theoretical Value","Experimental Limits"])
    Results_table.float_format = "4.2"

    for res in Analysis.results.keys():
        theoRes = SMSmethods.EvalRes(res,Analysis)  #Theoretical values for result and conditions
#        sys.exit()
        try:
            plot = Analysis.plots[res]
        except KeyError:    
            plot = []
                
#        print "Constraint: ",res,"  Conditions: ",Analysis.results[res]
        for imass in range(len(theoRes)):
            mass = theoRes[imass]['mass']
            tvalue = theoRes[imass]['result']
            conds = theoRes[imass]['conditions']
            massarray = copy.deepcopy(mass)
            sigmalimit = SMSgetlimit.GetPlotLimit(massarray,plot,Analysis)
            
            if sigmalimit and len(sigmalimit) > 0:
                Results_table.add_row([res,conds,mass,tvalue,sigmalimit[0]])
                for ilim in range(1,len(sigmalimit)):
                    Results_table.add_row(["","","","",sigmalimit[ilim]])
            else:
                Results_table.add_row([res,conds,mass,tvalue,sigmalimit])
            Results_table.add_row(["---------","---------","---------","---------","---------"])
#            print "Results for mass array ",mass," :"
#            print "    Theorical value = ",res,"Conditions = ",conds
#            print "    Experimental Limits: ",sigmalimit
#        print '\n'
#    print "-----------------------"
#    print '\n \n'    
    print(Results_table)
#    sys.exit()  # Just print first analysis (avoid too much output)       



sys.exit()

