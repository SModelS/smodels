#!/usr/bin/env python3

"""
.. module:: particlesLoader
   :synopsis: Loads the file Defining the list of Z2-even and Z2-odd particles to be used.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Matthias Wolf <matthias.wolf@wot.at>

"""

from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
import os,sys

def getParticlesFromSLHA(slhafile):
    """
    Defines BSM particles from the QNUMBERS blocks in the slhafile.
    OBS: If the QNUMBERS block does not include a 5-th entry, the particle
    will be assumed to be Z2-odd

    :param slhafile: Path to the SLHA file

    :return: List with Particle objects
    """

    from smodels.theory.particle import Particle
    from smodels.installation import installDirectory

    checkDirs =  [os.path.join(installDirectory(),"smodels","share","models"),installDirectory(),
                os.path.join(installDirectory(),"smodels")]

    filename = slhafile
    #If file does not exist, check if it is in any of the default folders:
    if not os.path.isfile(slhafile):
        for dirPath in checkDirs:
            if os.path.isfile(os.path.join(dirPath,slhafile)):
                filename = os.path.join(dirPath,slhafile)
                break

    if not os.path.isfile(filename):
        logger.error("Model file %s not found." %slhafile)
        raise SModelSError()

    logger.debug("Trying to define BSM particles from SLHA input file %s" %filename)

    #Read file and extract blocks 
    # (can not use pyslha, since it only allows for a single block per block name):
    with open(filename,'r') as f:
        data = f.read()
    data = data.lower()
    qnumberBlocks = []
    qBlock = False
    for l in data.splitlines():
        l = l.strip()
        # Ignore empty lines and comments
        if not l or l.startswith('#'):
            continue
        # Beginning of QNUMBERS block
        if l.startswith('block qnumbers'):
            qBlock = True
            qnumberBlocks.append([l])
            continue
        # If any other block is starting set qBlock to False
        elif l.startswith('block'):
            qBlock = False
            continue
        # If a cross-section or decay block is starting set qBlock to False
        elif l.startswith('xsection') or l.startswith('decay'):
            qBlock = False
            continue                   
        
        # If current block is not a qnumbers block, skip
        if not qBlock:
            continue
        # Add line to qnumbers block
        qnumberBlocks[-1].append(l)

    if not qnumberBlocks:
        logger.error("No QNUMBERS blocks were found in %s" %slhafile)
        raise SModelSError()

    #Build list of BSM particles:
    BSMList = []
    for b in qnumberBlocks:
        headerInfo = [x for x in b[0].replace('block','').replace('qnumbers','').split() 
                      if x != '#']
        if headerInfo[0].replace('-','').replace('+','').isdigit():
            pdg = eval(headerInfo[0])
        else:
            logger.error("Error obtaining PDG number from QNUMBERS block:\n %s \n" %b)

        if any(p.pdg == pdg for p in BSMList):
            logger.warning("Particle with pdg %i appears multiple times in QNUMBERS blocks")
            continue

        if len(headerInfo) > 1:
            label = headerInfo[1].strip()
        else:
            label = str(pdg)
            logger.debug("Could not find label for particle %i, will use its PDG number" %pdg)
        try:
            numbers = [l[:l.find('#')].lstrip().split() for l in b[1:]]
            numbers = dict([x for x in numbers if x])
            numbers = dict([[eval(x),eval(y)] for x,y in numbers.items()])
        except:
            logger.error("Error reading quantum numbers from block: \n %s \n" %b)
            continue
        if any(not x in numbers for x in [1,2,3]):
            logger.error("Missing quantum numbers in block:\n %s\n" %b)
            continue
        newParticle = Particle(Z2parity=-1, label=label, pdg=pdg,
                                eCharge=numbers[1]/3.,
                                colordim=numbers[3],
                                spin=(numbers[2]-1.)/2.)

        #Allows an additional quantum number defining the Z2-parity:
        if 11 in numbers:
            newParticle.Z2parity = int((-1)**int(numbers[11]))
        BSMList.append(newParticle)
        if numbers[4]: #Particle is not its own anti-particle
            newParticleC = newParticle.chargeConjugate()
            if any(p.pdg == newParticleC.pdg for p in BSMList):
                continue
            BSMList.append(newParticleC)

    return BSMList


def getParticlesFromModule(modelFile):
    """
    Reads the python model file and retrieves the list of BSM particles (BSMList)

    :param modelFile: Name/path to the python module containing the BSM particle definitions

    :return: a list of Particle objects
    """

    from smodels.installation import installDirectory
    from importlib import import_module
    fulldir = os.path.join(installDirectory(), "smodels", "share", "models")
    sys.path.insert(0, installDirectory())
    sys.path.insert(0, os.path.join(installDirectory(), "smodels"))
    sys.path.insert(0, fulldir)
    sys.path.insert(0, ".")

    logger.debug("Trying to load model file: %s" % modelFile)

    fname = modelFile[:]
    if "/" in fname:
        import shutil
        filename = os.path.basename(fname)
        shutil.copy(fname, filename)
    else:
        filename = fname

    
    if filename.endswith(".py"):
        importName = filename[:-3]
    else:
        importName = filename

    pM=import_module(importName, package='smodels')
    logger.debug("Found model file at %s" % pM.__file__)
    if filename != fname:
        os.remove(filename)
    BSMList = pM.BSMList

    return BSMList

def load():

    from smodels.tools.runtime import modelFile

    try:
        BSMList = getParticlesFromModule(modelFile)
    #If failed, assume the input is an SLHA file:
    except (ImportError,AttributeError,SModelSError):
        try:
            BSMList = getParticlesFromSLHA(modelFile)
        except SModelSError:
            logger.error("Could not load input model from %s. The file should be either a python module with particle definitions\
            or a SLHA file with QNUMBERS blocks." %modelFile)
            raise SModelSError()
    return BSMList

BSMList = load()
