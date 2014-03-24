#!/usr/bin/env python

"""
.. module:: xsecComputer
        :synopsis: The unit responsible for the computation of reference ("theory") \
            production cross sections

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
from tools import setPath
from physicsUnits import rmvunit,addunit
import os, commands, shutil
from theory import crossSection
import nllFast
import cStringIO
import logging

logger = logging.getLogger(__name__)


def addXSecToFile(sqrts,maxOrder,nevts,slhafile,lhefile=None,externaldir=None):
    """ Runs pythia at sqrts and compute SUSY cross-sections for the input SLHA file.
    Writes cross-sections to slha file.
    :param sqrts: sqrt{s} to run Pythia
    :param maxOrder: maximum order to compute the cross-section
    if maxOrder = 0, compute only LO pythia xsecs
    if maxOrder = 1, apply NLO K-factors fron NLLfast (if available)
    if maxOrder = 2, apply NLO+NLL K-factors fron NLLfast (if available)
    :param nevts: number of events for pythia run
    :param lhefile: LHE file. If None, do not write pythia output to file \
    if file does not exist, write pythia output to this file name \
    if file exists, read LO xsecs from this file (does not run pythia)
    :param externaldir: location of pythia6 and nllfast folders
    if not defined, it is set as current folder
    """

    if not externaldir: externaldir = os.getcwd() + '/external/'

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
        lheFile = runPythia(slhafile,nevts,rmvunit(sqrts,'TeV'),lhefile)
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
    import SModelS
    entry = "0  "+str(xsec.info.order)+"  0  0  0  0  "+str(rmvunit(xsec.value,'fb'))+" SModelS "+SModelS.version()

    return header + "\n" + entry

def runPythia(slhafile,nevts,sqrts,lhefile=None):
    """ run pythia_lhe with n events, at sqrt(s)=sqrts. Returns a file object with the lhe events
        :param slhafile: input SLHA file
        :param nevts: number of events to be generated
        :param sqrts: center of mass sqrt{s} (in TeV)
        :param lhefile: option to write LHE output to file. If None, do not write output to disk
    """

    import toolBox
    box=toolBox.ToolBox()
    tool=box.get("pythia6")
    pythiadir=tool.installDirectory()+"/"
    tool.checkInstallation()

#Check if slhafile, pythia_lhe, data and etc folders exist in pythiadir:
    if not os.path.isfile(slhafile):
        logger.error("File %s no found."%slhafile)
        return False
    else:
        shutil.copyfile(slhafile, pythiadir+"fort.61")
    #if not os.path.isdir(pythiadir):
    #    logger.error("pythia folder " + pythiadir +" not found in ")
    #    return False
    #elif not os.path.isfile(pythiadir+"/pythia_lhe") or not os.access(pythiadir+"/pythia_lhe",os.X_OK):
    #    logger.error("pythia_lhe file no found in "+pythiadir)
    #    return False
    #elif not os.path.isfile(pythiadir+"pythia.card"):
    #    logger.error("pythia.card file no found in "+pythiadir)
    #    return False

#Read pythia options:
    f=open(pythiadir+"pythia.card")
    pythiaOpts = f.readlines()
    f.close()
#Create pythia par file with options
    pythiaParFile = open(pythiadir+"pythia_card.dat","w")
    for option in pythiaOpts:
        if "MSTP(163)=" in option: option = "MSTP(163)=6\n"    #Switches output to screen, so no file is written to disk
        option = option.replace("NEVENTS",str(nevts)).replace("SQRTS",str(1000*sqrts))
        pythiaParFile.write(option)
    pythiaParFile.close()

    executable="%s/pythia_lhe" % pythiadir
    lhedata = commands.getoutput("cd %s; %s < pythia_card.dat" % (pythiadir, executable))
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

if __name__ == "__main__":
    """ called as script, we compute the cross section of a given slha file """
    import argparse, types, sys
    argparser = argparse.ArgumentParser(description='computes the cross section of a file')
    argparser.add_argument('file', type=types.StringType, nargs=1,
                           help='the slha or lhe file to compute cross section for')
    argparser.add_argument ( '-s', '--sqrts', nargs='+', action='append', help='sqrt(s) [TeV]. Can supply more than one value.', type=int, default=[] )
    argparser.add_argument ( '-e', '--nevents', help='number of events to be simulated.', type=int, default=100 )
    argparser.add_argument('-f','--tofile',help='write cross sections also to file', action='store_true')
    argparser.add_argument('-S','--slha',help='input file is slha file', action='store_true')
    argparser.add_argument('-n','--NLO',help='compute at the NLO level (default is LO)', action='store_true')
    argparser.add_argument('-N','--NLL',help='compute at the NLL level (takes precedence over NLL, default is LO)', action='store_true')
    argparser.add_argument('-L','--lhe',help='input file is lhe file', action='store_true')
    args=argparser.parse_args()

    import SModelS


    sqrtses=[item for sublist in args.sqrts for item in sublist]
    if len(sqrtses)==0: sqrtses=[8] ## default is: we compute for 8 tev!
    sqrtses.sort()
    sqrtses=set(sqrtses) ## unique values!
    order=0
    if args.NLO: order=1
    if args.NLL: order=2
    if order > 0:
        for sqrts in sqrtses:
            if not sqrts in [ 7,8,13,14,30,100]:
                logger.error("Cannot compute NLO or NLL xsecs for sqrts = %d TeV!" % sqrts)
                sqrtses.remove(sqrts)
    File=args.file[0]
    if not os.path.exists(File):
        logger.error("File ``%s'' does not exist." % File)
        sys.exit(1)
    if File[-5:].lower()==".slha" or args.slha:
        if args.tofile:
            logger.info("Computing slha cross section from %s, and adding to slha file." % File)
            for s in sqrtses:
                ss=addunit(s, 'TeV')
                external_dir = SModelS.installDirectory() + "/tools/external"
                addXSecToFile ( ss,order,args.nevents,File,externaldir=external_dir)
            logger.info("done.")
            sys.exit(0)
        else:
            logger.error("compute slha cross section, print out, but dont add to file. FIXME not yet implemented.")
            sys.exit(0)
    if File[-4:].lower()==".lhe" or args.lhe:
        if args.tofile:
            logger.error("Compute lhe section, and add to file. FIXME I guess we dont need this case?")
            sys.exit(0)
        else:
            logger.error("Compute lhe section, print out, but dont add to file. FIXME not yet implemented.")
            sys.exit(0)

