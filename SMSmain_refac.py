#!/usr/bin/env python

import time
from prettytable import PrettyTable
from theory import slhaDecomposer
from tools.PhysicsUnits import addunit
from tools import SMSPrettyPrinter
from tools.SMSPrettyPrinter import wrap
from experiment import smsanalysisFactory
from theory.theoryPrediction import theoryPredictionFor
import logging

logger = logging.getLogger(__name__)

def main():
    # useXsec = CrossSection.XSectionInfo()
    # useXsec.sqrts = addunit(8,'TeV')
    # useXsec.order = 2
    # useXsec.label = 'tev8'
    # UseXSecs = [useXsec]
    
    listOfAnalyses = smsanalysisFactory.load()
    printer=SMSPrettyPrinter.SMSPrettyPrinter()
    slhafile = "slha/andrePT4.slha"
    # lhefile = "lhe/ued_1.lhe"
    # lhefile = "lhe/TChiChipmSlepL_1.lhe"
    # nevts = 10000
    DoCompress = True
    DoInvisible = True
    minmassgap = addunit(5.,'GeV')
    sigmacut = addunit(0.1,'fb')
    t1 = time.time()
    SMSTopList = slhaDecomposer.decompose(slhafile,sigmacut,DoCompress,DoInvisible,minmassgap)
    # SMSTopList = lheDecomposer.decompose(lhefile,None,None,DoCompress,DoInvisible,minmassgap)
    
    EvTop_table = PrettyTable(["Topology","#Vertices", "#Insertions", "#Elements", "Sum of weights"])
    EvElement_table = PrettyTable(["Topology","Element","Particles B[0]","Particles B[1]", "Masses B[0]","Masses B[1]","Element Weight"])
    
    eltot = 0
    # totweight = []
    # Print Results:
    # for i in range(len(SMSTopList)):
    for i,topo in enumerate(SMSTopList):
        sumw = topo.getTotalWeight().getDictionary()
        EvTop_table.add_row([i,topo.vertnumb,topo.vertparts,len(topo.ElList),wrap(printer.pformat(sumw),width=30)])
        eltot += len(topo.ElList)
    
     
          
    #Print element list for Topology[i]:  
        if i == 0:
            for j,el in enumerate(topo.ElList):
                if el.getParticles() != [[['b','b']],[['b','b']]]:
                    continue
                EvElement_table.add_row([i,j,el.getParticles()[0],el.getParticles()[1],wrap(printer.pformat(el.getMasses()[0]),width=25),wrap(printer.pformat(el.getMasses()[1]),width=25),wrap(printer.pformat(el.weight.getDictionary()),width=30)])
            EvElement_table.add_row(["---","---","---","---","---","---","---"])  
    
         
    print "Number of Global topologies = ",len(SMSTopList)      
    print(EvTop_table)
    print "Total Number of Elements = ",eltot
    print "Total weight = ",SMSTopList.getTotalWeight()
    # print(EvElement_table)
    
    print '\n \n \n'
    print 'slhaDecomposer done in',time.time()-t1,'s'
    
    for ana in listOfAnalyses:    
        preds = theoryPredictionFor(ana,SMSTopList)
        if not preds:
            continue
        print ana.label
        for pred in preds:
            print 'mass=',pred.mass
            print 'theory prediction=',pred.value
            print 'theory conditions:'
            if not pred.conditions:
                print pred.conditions
            else:
                for cond in pred.conditions:
                    print pred.conditions[cond]
            print '\n'
    

if __name__ == '__main__':
    main()
