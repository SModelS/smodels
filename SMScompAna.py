#!/usr/bin/env python

import sys
from Experiment import SMSResults
import SMSanalyses ## SuK/AL descriptions
import SMSAnalysisFactory ## UND/WW descriptions

#Creat analyses list:
ListOfAnalysesU = SMSAnalysisFactory.load()
ListOfAnalysesA = SMSanalyses.load()



#Compare analyses:
print "-----------------------------"
print "Comparing matching entries:"    
for Ana in ListOfAnalysesA:
    for key in Ana.results.keys():
        conA = Ana.results[key].replace(" ","")
        conA = conA.split("]]],[[[")
        TxA = Ana.plots[key][0]
        plots = Ana.plots[key][1]
        diff = False              
        for plot in plots:
            match = False
            strplot = plot+':'+TxA
            for AnaU in ListOfAnalysesU:
                if AnaU.label != strplot: continue
                for keyU in AnaU.results.keys():
                    if AnaU.plots[keyU][0] == TxA and plot in AnaU.plots[keyU][1]:
                        match = True
                        conB = AnaU.results[keyU].replace(" ","")
                        conB = conB.split("]]],[[[")                        
                        if keyU.replace(" ","") != key.replace(" ","") or conA != conB:
                            
                            if not diff:
                                diff = True
                                print "\nConditions or constraint differ for",TxA
                                print "from SMSAnalyses.py:"
                                print key.replace(" ","")
                                for ix in range(len(conA)):
                                    x = conA[ix]
                                    if x == "None":
                                        print x
                                        continue
                                    if ix != 0 and ix != len(conA)-1:
                                        print "[[["+x+"]]]"
                                    elif ix == 0:
                                        print x+"]]]"
                                    elif ix == len(conA)-1:
                                        print "[[["+x

                            
                            
                            print "from Factory entry",AnaU.label,":"
                            print keyU.replace(" ","")
                            for ix in range(len(conB)):
                                x = conB[ix]
                                if x == "None":
                                    print x
                                    continue

                                if ix != 0 and ix != len(conB)-1:
                                    print "[[["+x+"]]]"
                                elif ix == 0:
                                    print x+"]]]"
                                elif ix == len(conB)-1:
                                    print "[[["+x
#                if match:
#                    AnaU.label = "done"
            
            if not diff and match:
                print "Analyses",plot,"for topology",TxA,": OK"
            if not match:                
                print "Analyses",plot,"for topology",TxA,"not found in Factory"
        if diff: print "\n"

print "\n\n"  

#adone = 0
#for Ana in ListOfAnalysesU:
#    if Ana.label == "done":
#        adone += 1
#    else:
#        print Ana.label
#        
#print "done=",adone
#sys.exit()            


print "-----------------------------"
print "Checking for missing analysis:"
Ucount = 0
for Ana in ListOfAnalysesU:
    TxU = Ana.label.split(":")[1]
    plotU = Ana.label.split(":")[0]
    match = False
    for key in AnaU.results.keys():
        Ucount += 1
    for AnaA in ListOfAnalysesA:
        for key in AnaA.results.keys():
            TxA = AnaA.plots[key][0]
            if TxA != TxU: continue
            plots = AnaA.plots[key][1]
            for plot in plots:
                if plot == plotU: match = True
       
    if not match:
        if SMSResults.isPublic(plotU):
            print "Analysis",plotU," for ",TxU," missing in SMSanalyses.py"
        else:
            print "Analysis",plotU," for ",TxU," missing in SMSanalyses.py (NOT PUBLIC)"


print "-----------------------------"
Acount = 0
for Ana in ListOfAnalysesA:
    for key in Ana.results.keys():
        TxA = Ana.plots[key][0]
        plots = Ana.plots[key][1]
        for plot in plots:
            Acount += 1
            match = False
            strplot = plot+':'+TxA
            for AnaU in ListOfAnalysesU:
                if AnaU.label != strplot: continue
                for keyU in AnaU.results.keys():
                    if AnaU.plots[keyU][0] == TxA and plot in AnaU.plots[keyU][1]:
                        match = True
            if not match:
                print "Analysis",plot," for",TxA," missing in Factory"


print "\n"  

print "-----------------------------"
print "Checking for duplicate entries:"
for iAna in range(len(ListOfAnalysesU)-1):
    Ana = ListOfAnalysesU[iAna]
    TxU = Ana.label.split(":")[1]
    plotU = Ana.label.split(":")[0]
    for jAna in range(iAna+1,len(ListOfAnalysesU)):
        AnaB = ListOfAnalysesU[jAna]
        TxUB = AnaB.label.split(":")[1]
        plotUB = AnaB.label.split(":")[0]
        if TxUB == TxU and plotUB == plotU:
            print "Duplicate entry for analysis",plotU,"and topology",TxU," in Factory"
            Ucount -= 1
            
for iAna in range(len(ListOfAnalysesA)-1):
    Ana = ListOfAnalysesA[iAna]
    for key in Ana.results.keys():
        TxA = Ana.plots[key][0]
        plots = Ana.plots[key][1]
        for plot in plots:
            if plots.count(plot) > 1:
                print "Duplicate plot entry",plot,"in Analyses",Ana.label," in SMSanalyses.py"
                Acount -= 1
            for jAna in range(len(ListOfAnalysesA)):
                AnaB = ListOfAnalysesA[jAna]
                for key2 in AnaB.results.keys():
                    if iAna == jAna and key == key2: continue
                    TxAB = AnaB.plots[key2][0]
                    plotsB = AnaB.plots[key2][1]
                    for plotB in plotsB:
                        if TxAB == TxA and plotB == plot:
                            print "Duplicate entry for analysis",plot,"and topology",TxA," in SMSanalyses.py"
                            Acount -= 1


print "\n"            
print "-----------------------------"
print "Total of (non-duplicated) analyses (including possible empty entries) in:"
print "  SMSanalyses.py = ",Acount
print "  Factory = ",Ucount


sys.exit()

