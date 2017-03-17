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
logger.setLevel(logging.DEBUG)


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

def debugFile(slhafile,nevts=50000,forceDegenerate=False):
    #Individual file debugging:
    
    if forceDegenerate:
        f = pyslha.readSLHAFile(slhafile)
        masses = f.blocks['MASS']
        squarksMasses = [abs(mass) for pid,mass in masses.items() if pid in squarks]
        avgmass = sum(squarksMasses)/len(squarksMasses)
        if abs(max(squarksMasses)-avgmass) > 0.1 or  abs(min(squarksMasses)-avgmass) > 0.1:
            for pid in squarks:
                f.blocks['MASS'][pid] = avgmass
                
            slhaF,slhafile_new = tempfile.mkstemp(suffix='.slha', dir='./')
            os.write(slhaF,f.write())            
            os.close(slhaF)
            logger.warning("Testing degenerate squarks for %s with average mass %s" %(slhafile,avgmass))
        comp = debugFile(slhafile_new, nevts=nevts,forceDegenerate=False)
        os.remove(slhafile_new)
        return comp
    
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
    
    
    #Remove the antisbottom-gluino xsec (seems to be missing in Pythia 8):
#    if (-1000005, 1000021) in w6:
#        w6.pop((-1000005, 1000021))
    #Remove the antisdown-gluino xsec (seems to be missing in Pythia 8):
#    if (-1000001, 1000021) in w6:
#        w6.pop((-1000001, 1000021))
#    if (-2000001, 1000021) in w6:
#        w6.pop((-2000001, 1000021))
#    if (-1000003, 1000021) in w6:
#        w6.pop((-1000003, 1000021))        
    #Remove the antisbottom-gluino xsec (seems to be missing in Pythia 8):
#    if (-1000024, 1000021) in w6:
#        totxsec = w6[(-1000024, 1000021)].values()[0]
#        if (1000021, 1000024) in w6:
#            totxsec += w6[(1000021, 1000024)].values()[0]
#        w6.pop((-1000024, 1000021))
#        w6[(1000021, 1000024)] = {'8 TeV (LO)' : totxsec}    
    
    comp = compareXSections(w6,w8,nevts,relError=0.1)

        
    return comp 


def checkFiles(slha6,slha8,Nevents = 50000):
    
    xsecs8 = getXsecFromSLHAFile(slha8)
    if not xsecs8:
        return (slha6,slha8,True)
    xsecs6 = getXsecFromSLHAFile(slha6)    
    w6 = xsecs6.getXsecsFor('8 TeV (LO)').getDictionary()
    w8 = xsecs8.getXsecsFor('8 TeV (LO)').getDictionary()


    #Remove the antisbottom-gluino xsec (seems to be missing in Pythia 8):
#    if (-1000005, 1000021) in w6:
#        w6.pop((-1000005, 1000021))
    #Remove the antisdown-gluino xsec (seems to be missing in Pythia 8):
#    if (-1000001, 1000021) in w6:
#        w6.pop((-1000001, 1000021))
#    if (-2000001, 1000021) in w6:
#        w6.pop((-2000001, 1000021))
#    if (-1000003, 1000021) in w6:
#        w6.pop((-1000003, 1000021))
    #Remove the antisbottom-gluino xsec (seems to be missing in Pythia 8):
#    if (-1000024, 1000021) in w6:
#        totxsec = w6[(-1000024, 1000021)].values()[0]
#        if (1000021, 1000024) in w6:
#            totxsec += w6[(1000021, 1000024)].values()[0]
#        w6.pop((-1000024, 1000021))
#        w6[(1000021, 1000024)] = {'8 TeV (LO)' : totxsec}

    comp = compareXSections(w6,w8,Nevents,relError=0.1)
    
    if not comp:
        return (slha6,slha8,True)
    else:
        #Check for the case of degenerate squarks
        comp = debugFile(slha6,nevts=Nevents,forceDegenerate=True)
        
    if comp:
        logger.error(str(comp))
        return (slha6,slha8,False)
    else:        
        return (slha6,slha8,True)


def checkFolders(slha6Folder,slha8Folder):
    
    badFiles = []
    
    slhaFiles = []
    for slha6 in glob.glob(os.path.join(slha6Folder,'*.slha')):

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
        logger.setLevel(logging.DEBUG)
        print debugFile(sys.argv[1])
#        print debugFile(sys.argv[1],forceDegenerate=True)
    elif len(sys.argv) >= 3:
        if len(sys.argv) == 4:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.WARNING)
        print checkFolders(sys.argv[1],sys.argv[2])
        
            
