#!/usr/bin/env python

"""
.. module:: nllFast
    :synopsis: This module provides methods to access the nllfast grid and
    compute k-factors (when available) to SUSY pair production cross-sections

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import commands, os, sys
import set_path
from theory import pyslha2 as pyslha
from tools.physicsUnits import rmvunit
import numpy
import logging
import operator
from tools import toolBox
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

squarks = [1000001,2000001,1000002,2000002,1000003,2000003,1000004,2000004]
antisquarks = map ( operator.neg, squarks )
third = [1000005,2000005,1000006,2000006]
gluinos = [1000021]

def getKfactorsFor(pIDs,sqrts,slhafile,pdf='cteq'):
    """Reads the NLLfast grid and returns a pair of k-factors (NLO and NLL) for
        the production of a pID pair.  If NLLfast does not contain the process,
        return k-factors = 1.  Uses the slhafile to obtain the SUSY
        spectrum.
    """


    if not os.path.isfile(slhafile):
        logger.error("SLHA file %s not found" %slhafile)
        return False
        
#Get process name (in NLLfast notation)
    process = getProcessName(pIDs)
    if not process: return (1.,1.)  #Return k-factors =1, if NLLfast does not have the process

#Obtain relevant masses:
    readfile=pyslha.readSLHAFile(slhafile)
    gluinomass = abs(readfile[0]['MASS'].entries[1000021])
    squarkmass = sum([abs(readfile[0]['MASS'].entries[pid]) for pid in squarks])/8.
    pid1,pid2 = sorted(pIDs)    
    if pid1 in antisquarks and pid2 in squarks: squarkmass = (abs(readfile[0]['MASS'].entries[abs(pid1)])+abs(readfile[0]['MASS'].entries[pid2]))/2        
    elif pid1 in squarks and pid2 in squarks: squarkmass = (abs(readfile[0]['MASS'].entries[pid1])+abs(readfile[0]['MASS'].entries[pid2]))/2
    elif abs(pid1) == pid2 and pid2 in third: squarkmass = abs(readfile[0]['MASS'].entries[abs(pid1)])

    # Set up NLLfast run, the old way
    sqrtS = float(rmvunit(sqrts,'TeV'))
    energy = str(int(sqrtS))+'TeV'
    #if sqrtS == 8.: nllpath = basedir + "/nllfast-2.1"
    #elif sqrtS == 7.: nllpath = basedir + "/nllfast-1.2"
    #else:
    #    logger.warning("No NLLfast data for sqrts = " + energy)
    #    return (1.,1.)
    toolname="nllfast%d" % int(sqrtS)
    box=toolBox.ToolBox()
    tool=box.get(toolname)
    if tool==None:
        logger.warning("No NLLfast data for sqrts = " + str(sqrts) )
        return (1.,1.)
    nllpath=tool.installDirectory()
    nllexec=tool.pathOfExecutable()
    ret=tool.checkInstallation()
    #if not os.path.isfile(nllexec):
    #    logger.error("Missing NLL executable: " + nllexec)
    #    return False
    if process == "st": 
        nll_run = "./nllfast_"+energy+" %s %s %s" % (process,pdf,squarkmass)
    else: 
        nll_run ="./nllfast_"+energy+" %s %s %s %s" % \
                    (process,pdf,squarkmass,gluinomass)

#Run NLLfast
    nll_output = runNLLfast(nll_run,nllpath)

#If run was successful, return k-factors:    
    if "K_NLO" in nll_output:
        kFacs = getKfactorsFrom(nll_output)  #NLLfast ran ok, try to get the k-factors 
        if not kFacs or min(kFacs) <= 0.:
            logger.warning("Error obtaining k-factors")
            return (1.,1.)
        else: return kFacs
#If run was not successful, check for decoupling error messages:
    elif not "too low/high" in nll_output.lower():
        logger.warning("Error running NLLfast")
        return (1.,1.)
     
#Deal with decoupling regimes:
    gluino_dcp,squark_dcp = False,False
    if "too low/high gluino" in nll_output.lower(): gluino_dcp = True
    elif "too low/high squark" in nll_output.lower(): squark_dcp = True

#If produced particles are too heavy, return k-factors = 1        
    if gluino_dcp and ('g' in process or gluinomass < 500.):
        logger.warning("Gluino mass out of NLLfast grid for " + process)
        return (1.,1.)
    elif squark_dcp and (process != 'gg' or squarkmass < 500.):
        logger.warning("Squark mass out of NLLfast grid for " + process)
        return (1.,1.)

#If virtual particles are too heavy, interpolate to the decoupling limit        
    kFacsVector = []
    xmass = max(squarkmass,gluinomass) #To be interpolated
    while len(kFacsVector) < 2 and xmass > 500.:
        xmass -= 100.  #Reduce decoupled mass, until NLLfast produces results
        if squark_dcp: nll_run = "./nllfast_"+energy+" %s %s %s %s" % (process,pdf,xmass,gluinomass)
        elif gluino_dcp: nll_run = "./nllfast_"+energy+" %s %s %s %s" % (process,pdf,squarkmass,xmass)        
        nll_output = runNLLfast(nll_run,nllpath)                
        if "K_NLO" in nll_output:
            kfacs = getKfactorsFrom(nll_output)
            if kfacs and min(kfacs) > 0.: kFacsVector.append([xmass,kfacs])

    if len(kFacsVector) < 2:
        logger.warning("Not enough points for interpolation in the decoupling limit")
        return (1.,1.)
    else:
        if kFacsVector[0][0]/min(squarkmass,gluinomass) > 10.:
            kFacs = kFacsVector[0][1]  #Decoupling limit is satisfied, do not interpolate
        else:
            kFacs = interpolateKfactors(kFacsVector,max(squarkmass,gluinomass))  #Interpolate k-factors
    return kFacs



def getProcessName(pIDs):
    """Returns the process name (in NLLfast notation) for the pair production of pIDs.
    If the particle ID pair is not contained in NLLfast, returns None."""
    
    pid1,pid2 = sorted(pIDs)
    process = None
    
#Obtain the type of process:
#gluino-gluino production = gg, squark-antisquark = sb, squark-squark = ss, squark-gluino = sg, 
#stop-stop or sbottom-sbottom = st   
    if pid1 in antisquarks and pid2 in squarks: process = 'sb'      
    elif pid1 in squarks and pid2 in squarks: process = 'ss'
    elif pid1 == pid2 and pid1 in gluinos: process = 'gg'
    elif pid1 in squarks and pid2 == 1000021 or pid2 in squarks and pid1 == 1000021: process = 'sg'
    elif abs(pid1) == pid2 and pid2 in third: process = 'st'
    
    return process
        
 
def runNLLfast(nll_run,nllpath):
    """Simple method to run NLLfast with command nll_run at nllpath. Returns NLLfast output as a string. """
    
    current_dir = os.getcwd()
    os.chdir(nllpath)
    nll_output = commands.getoutput(nll_run)
    os.chdir(current_dir)
    return nll_output

    
def getKfactorsFrom(output):
    """Simple method to read NLLfast output and return the k-factors."""
    
#First check if run was successful:
    if not output: return False
    else:
        lines = output.split('\n')
        il = 0
        line = lines[il]
        while not "K_NLO" in line and il < len(lines)-2:
            if "process" in line: 
                process = line[line.find("process:")+8:].replace(" ","")
            il += 1
            line = lines[il]        
        if not process: return False
        line = lines[il]  #Line with header        
        nmass = line.count('GeV')  #Count number of mass entries
        line = lines[il+2]  #Line with values
        data = [eval(x) for x in line.split()]
        if len(data) != nmass + 11: return False
        else: kFacs = data[9+nmass:]

    return kFacs


def interpolateKfactors(kFacsVector,xval):
    """Simple method to interpolate a list of k-factor 
        values from kFacsVector = [[x0,[k1,k2,..]], [x1,[k1,k2,..],...].
    Returns the list of interpolated k-factor values at x-value xval."""
    
    kFacs = []
    
    xpts = [x[0] for x in kFacsVector]
    ypts = [x[1] for x in kFacsVector]
    coeffs = numpy.matrix.transpose(numpy.polyfit(xpts, ypts, len(xpts)-1))    
    for ik in range(len(ypts[0])):
        kfac = 0.
        for ip,coeff in enumerate(coeffs[ik]): kfac += coeff*xval**(len(xpts)-1-ip)
        if kfac <= 0.: kfac = 1.
        kFacs.append(kfac)
    
    return kFacs

if __name__ == "__main__":
    """ called as script, we get the k factors for some pid pair """
    import set_path
    import SModelS
    from physicsUnits import addunit
    slhaF=SModelS.installDirectory()+"inputFiles/slha/T1.slha"
    kNLO,kNLL = getKfactorsFor((1000021,1000021),addunit(13.,"TeV"),slhaF)
    print "[nllFast.py] nlo,nll=",kNLO,kNLL
