#!/usr/bin/env python

"""
.. module:: XSecComputer
        :synopsis: The unit responsible for the computation of reference ("theory") \
            production cross sections
        
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""
from Tools.PhysicsUnits import rmvunit
import os, commands, shutil
import crossSection, nllFast
import cStringIO
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def addXSecToFile(sqrts,maxOrder,nevts,slhafile,lhefile=None,basedir=None):
    """ Runs pythia at sqrts and compute SUSY cross-sections for the input SLHA file.
    Write cross-sections to file
    :param sqrts: sqrt{s} to run Pythia
    :param maxOrder: maximum order to compute the cross-section
    if maxOrder = 0, compute only LO pythia xsecs
    if maxOrder = 1, apply NLO K-factors fron NLLfast (if available)
    if maxOrder = 2, apply NLO+NLL K-factors fron NLLfast (if available)
    :param nevts: number of events for pythia run  
    :param lhefile: LHE file. If None, do not write pythia output to file \
    if file does not exist, write pythia output to this file name \
    if file exists, read LO xsecs from this file (does not run pythia)
    :param basedir: is where pythia_lhe is to be found, as well as the data and etc folders \
    if not defined, it is set as current folder
    """

    if not basedir: basedir = os.getcwd()

#Check if SLHA file exists
    if not os.path.isfile(slhafile):
        logger.error("SLHA file not found.")
        return None
#Check if file already contain cross-section blocks:
    infile = open(slhafile,'r')
    slhadata = infile.read()
    infile.close()
    if 'XSECTION' in slhadata:
        logger.warning("SLHA file already contains a XSECTION block. Appending new cross-sections.")

#Check i lhefile exists:    
    if lhefile:
        if os.path.isfile(lhefile): logger.warning("Using LO cross-sections from "+lhefile)                        
        else: logger.warning("Writing pythia LHE output to "+lhefile) 
#Get file with lhe Events:
    if not lhefile or not os.path.isfile(lhefile):
        lheFile = runPythia(slhafile,nevts,rmvunit(sqrts,'TeV'),lhefile,basedir)      
    else:
        lheFile = open(lhefile,'r')
            
#Get LO cross-sections from LHE events     
    LOxsecs = crossSection.getXsecFromLHEFile(lheFile)
    xsecs = LOxsecs
#Set cross-section label and order:
    wlabel = str(int(rmvunit(sqrts,'TeV')))+' TeV'
    if maxOrder == 0:  wlabel += ' (LO)'
    elif maxOrder == 1: wlabel += ' (NLO)'
    elif maxOrder == 2: wlabel += ' (NLL)'
    for ixsec,xsec in enumerate(xsecs):
        xsecs[ixsec].info.label = wlabel
        xsecs[ixsec].info.order = maxOrder
#If maxOrder > 0, apply k-factors:
    if maxOrder > 0:
        pIDs = LOxsecs.getPIDpairs()   # Get particle ID pairs for all xsecs
        for pID in pIDs:
            k = 1.
            kNLO,kNLL = nllFast.getKfactorsFor(pID,sqrts,slhafile)
            if maxOrder == 1 and kNLO: k = kNLO
            elif maxOrder == 2 and kNLL and kNLO: k = kNLO*kNLL
            else:
                logger.warning("Unkown xsec order, using NLL+NLO k-factor (if available)")
                k = kNLO*kNLL
            k = float(k)
            for i,xsec in enumerate(xsecs):
                if set(xsec.pid) == set(pID): xsecs[i] = xsec*k   #Apply k-factor

#Write cross-sections to file    
    outfile = open(slhafile,'append')
    for xsec in xsecs: outfile.write(xsecToBlock(xsec,inPDGs=(2212,2212),comment="Nevts="+str(nevts))+"\n")
    outfile.close()    
    
    return True

def xsecToBlock(xsec,inPDGs=(2212,2212),comment=None):
    """Generates a string for a XSECTION block in the SLHA format from a XSection() object.
    inPDGs defines the PDGs of the incoming states (default = 2212,2212).
    comment is added at the end of the header as a comment"""
    
    if type(xsec) != type(crossSection.XSection()):
        logger.error("Wrong input")
        return False
    header = "XSECTION  "+str(rmvunit(xsec.info.sqrts,'GeV'))  #Sqrt(s) in GeV
    for pdg in inPDGs: header += " "+str(pdg)  #PDGs of incoming states
    header += " "+str(len(xsec.pid))   #Number of outgoing states
    for pid in xsec.pid: header += " "+str(pid)  #PDGs of outgoing states
    header += "   # "+str(comment)   #Comment
    entry = "0  "+str(xsec.info.order)+"  0  0  0  0  "+str(rmvunit(xsec.value,'fb'))+" SModelS "+version()
    
    return header + "\n" + entry
    

def version():
    """Returns SModelS version. To be replaced with a proper version method."""
    return "1.0"

def runPythia(slhafile,nevts,sqrts,lhefile=None,basedir=None):
    """ run pythia_lhe with n events, at sqrt(s)=sqrts. Returns a file object with the lhe events
        :param slhafile: input SLHA file
        :param nevts: number of events to be generated
        :param sqrts: center of mass sqrt{s} (in TeV)
        :param lhefile: option to write LHE output to file. If None, do not write output to disk
        :param basedir: is where pythia_lhe is to be found, as well as the data and etc folders\
        if not defined, it is set as current folder
    """
    
    if not basedir: basedir = os.getcwd()
    
    datadir = basedir+"/data/"
    etcdir = basedir+"/etc/"
#Check if slhafile, pythia_lhe, data and etc folders exist in basedir:
    if not os.path.isfile(slhafile):
        logger.error("File %s no found."%slhafile)
        return False        
    else:
        shutil.copyfile(slhafile, datadir+"fort.61")    
    if not os.path.isfile(basedir+"/pythia_lhe") or not os.access(basedir+"/pythia_lhe",os.X_OK):
        logger.error("pythia_lhe file no found in "+basedir)
        return False
    elif not os.path.isdir(datadir):
        logger.error("data folder no found in "+basedir)
        return False
    elif not os.path.isfile(etcdir+"external_lhe.template"):
        logger.error("external_lhe.template file no found in "+etcdir)
        return False    
    
#Read pythia options:    
    f=open(etcdir+"external_lhe.template")
    pythiaOpts = f.readlines()
    f.close()
#Create pythia par file with options    
    pythiaParFile = open(datadir+"external_lhe.dat","w")
    for option in pythiaOpts:
        if "MSTP(163)=" in option: option = "MSTP(163)=6\n"    #Switches output to screen, so no file is written to disk
        option = option.replace("NEVENTS",str(nevts)).replace("SQRTS",str(1000*sqrts))
        pythiaParFile.write(option)
    pythiaParFile.close()
    
    executable="%s/pythia_lhe" % basedir
    lhedata = commands.getoutput("cd %s; %s < external_lhe.dat" % (datadir, executable))    
    if not "<LesHouchesEvents" in lhedata:
        logger.error("LHE events not found in pythia output")
        return False

#Generates file object with lhe events:
    if lhefile:
        lheFile = open(lhefile,'w')
        lheFile.write(lhedata)
        lheFile.close()
        lheFile = open(lhefile,'r')       
    else:
        lheFile = cStringIO.StringIO(lhedata)   #Creates memory only file object

    return lheFile

