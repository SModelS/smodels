#!/usr/bin/env python

"""
.. module:: pythia8Comp
   :synopsis: Compares the output of pythia8 cross-sections with pythia 6 cross-sections.
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

import sys,os,glob

sys.path.insert(0,"../")
from smodels.theory.crossSection import getXsecFromSLHAFile
from smodels.tools import xsecComputer
from smodels.tools.xsecComputer import LO
from smodels.tools.physicsUnits import TeV,fb,pb
from unum import Unum
import logging
import pyslha
from math import sqrt
import tempfile
import multiprocessing

Unum.VALUE_FORMAT = "%.8f"  

squarks = [1000000 + i for i in range(1,5)]
squarks += [2000000 + i for i in range(1,5)]

logger = logging.getLogger()
logger.setLevel(logging.WARNING)


def compareXSections(dictA,dictB,nevts,relError = 0.1):    

    missingXsecs = set(dictA.keys()).symmetric_difference(set(dictB.keys()))
    commonXsecs = set(dictA.keys()).intersection(set(dictB.keys()))
    totXsecA = sum([x.values()[0].asNumber(fb) for x in dictA.values()])
    totXsecB = sum([x.values()[0].asNumber(fb) for x in dictB.values()])
    totXsec = max(totXsecA,totXsecB)*fb
    mcError = totXsec/sqrt(float(nevts))
    
    diffXsecs = []

    for xsec in commonXsecs:
        diff = abs(dictA[xsec].values()[0] - dictB[xsec].values()[0])
        if diff > 2.*mcError and (diff/(abs(dictA[xsec].values()[0] + dictB[xsec].values()[0]))).asNumber() > relError:
            logger.debug('Cross-section for %s differ by (%s +- %s) fb' %(str(xsec),str(diff.asNumber(fb)),str(mcError.asNumber(fb))))
            logger.debug('   %s   %s' %(dictA[xsec].values()[0],dictB[xsec].values()[0]))
            diffXsecs.append(xsec)

    for xsec in missingXsecs:
        if xsec in dictA:
            if dictA[xsec].values()[0] > 2.*mcError:
                logger.debug('Cross-section for %s missing' %str(xsec))
                logger.debug('    %s' %(dictA[xsec].values()[0]))
                diffXsecs.append(xsec)
        else:
            if dictB[xsec].values()[0] > 2.*mcError:
                logger.debug('Cross-section for %s missing' %str(xsec))
                logger.debug('    %s' %(dictB[xsec].values()[0]))
                diffXsecs.append(xsec)
                
    return diffXsecs

def debugFile(slhafile,nevts=10000):
    #Individual file debugging:
    computer6 = xsecComputer.XSecComputer(LO, nevts, 6)
    computer8 = xsecComputer.XSecComputer(LO, nevts, 8)
    w6 = computer6.compute(8*TeV, slhafile, pythiacard = './my_pythia6.card').getDictionary()
    w8 = computer8.compute(8*TeV, slhafile, pythiacard = './my_pythia8.cfg').getDictionary()
    
#     print 'Pythia 6:'
#     for key,val in sorted(w6.items()):
#         print key,val.values()[0]
#             
#     print 'Pythia 8:'
#     for key,val in sorted(w8.items()):
#         print key,val.values()[0]
    
    comp = compareXSections(w6,w8,nevts,relError=0.1)
        
    return comp 


def checkFiles(slha6,slha8,Nevents = 50000):
    
    xsecs8 = getXsecFromSLHAFile(slha8)
    if not xsecs8:
        return (slha6,slha8,True)
    xsecs6 = getXsecFromSLHAFile(slha6)    
    w6 = xsecs6.getXsecsFor('8 TeV (LO)').getDictionary()
    w8 = xsecs8.getXsecsFor('8 TeV (LO)').getDictionary()

    comp = compareXSections(w6,w8,Nevents,relError=0.1)
    
    if not comp:
        return (slha6,slha8,True)
    
    
    #Remove degenerate squarks
    f = pyslha.readSLHAFile(slha6)
    masses = f.blocks['MASS']
    squarksMasses = [abs(mass) for pid,mass in masses.items() if pid in squarks]
    avgmass = sum(squarksMasses)/len(squarksMasses)
    if abs(max(squarksMasses)-avgmass) > 0.1 or  abs(min(squarksMasses)-avgmass) > 0.1:
        for pid in squarks:
            f.blocks['MASS'][pid] = avgmass
            
        slhaF,slhafile = tempfile.mkstemp(suffix='.slha', dir='./')
        os.write(slhaF,f.write())            
        os.close(slhaF)
        logger.warning("Testing degenerate squarks for %s with average mass %s" %(slha6,avgmass))
        comp = debugFile(slhafile, nevts=Nevents)
        os.remove(slhafile)
    
    #Remove the antisbottom-gluino xsec (seems to be missing in Pythia 8):
    if (-1000005, 1000021) in comp:
        comp.remove((-1000005, 1000021))
    
    if comp:
        logger.error(str(comp))
        return (slha6,slha8,False)
    else:        
        return (slha6,slha8,True)


def checkFolders(slha6Folder,slha8Folder):
    
    badFiles = []
    
    slhaFiles = []
    for slha6 in glob.glob(os.path.join(slha6Folder,'*.slha')):
#         if not '21725753' in slha6: continue
        slha8 = slha6.replace('.slha','_new.slha').replace(slha6Folder,slha8Folder)
        if not os.path.isfile(slha8):
            continue
        slhaFiles.append((slha6,slha8))          
    
    
    pool = multiprocessing.Pool(processes=5)
    jobs = [pool.apply_async(checkFiles,args=slha) for slha in slhaFiles]
    for job in jobs:
        slha6,slha8,good = job.get(timeout=500)
        if not good:
            logger.error('Files %s and %s differ' %(slha6,slha8))
            badFiles.append(slha6)
        else:
            os.remove(slha8)
            logger.warning('File %s OK' %os.path.basename(slha6))
            
    return badFiles            
            
if __name__ == "__main__":    

    if len(sys.argv) == 2:
        print debugFile(sys.argv[1])
    elif len(sys.argv) == 3:
        print checkFolders(sys.argv[1],sys.argv[2])
        
            