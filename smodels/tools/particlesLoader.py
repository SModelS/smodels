#!/usr/bin/env python3

"""
.. module:: particlesLoader
   :synopsis: Loads the file defining the list of BSM particles to be used.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import os
import sys
from smodels.base.exceptions import SModelSBaseError as SModelSError
from smodels.base.smodelsLogging import logger
from smodels.base.particle import Particle
from smodels.installation import installDirectory
from smodels.share.models.SMparticles import SMList
from importlib import import_module



def getParticlesFromSLHA(slhaData: str) -> list[Particle]:
    """
    Defines BSM particles from the QNUMBERS blocks in the slhafile.

    :param slhaData: Path to the SLHA file or the a string with the SLHA data

    :return: List with Particle objects
    """


    # Create a list of SM PDGs, so if a QNUMBERS block for a SM particle
    # is present, it will be ignored.
    SMpdgs = set()
    for ptc in SMList:
        if isinstance(ptc.pdg,(int,float)):
            SMpdgs.add(int(abs(ptc.pdg)))
        else:
            for pdg in ptc.pdg:
                SMpdgs.add(int(abs(pdg)))
    SMpdgs = list(SMpdgs)

    if 'block qnumbers' in slhaData.lower():
        data = slhaData[:]
    else:
        filename = slhaData
        #If file does not exist, check if it is in any of the default folders:
        checkDirs = [os.path.join(installDirectory(), "smodels", "share", "models"),
                     installDirectory(),
                     os.path.join(installDirectory(), "smodels")]
        if not os.path.isfile(slhaData):
            for dirPath in checkDirs:
                if os.path.isfile(os.path.join(dirPath, slhaData)):
                    filename = os.path.join(dirPath, slhaData)
                    break

        if not os.path.isfile(filename):
            logger.error(f"Model file {slhaData} not found.")
            raise SModelSError()

        logger.debug(f"Trying to define BSM particles from SLHA input file {filename}")

        #Read file and extract blocks:
        with open(filename, 'r') as f:
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
        logger.error("No QNUMBERS blocks were found")
        raise SModelSError()

    #Build list of BSM particles:
    BSMList = []
    for b in qnumberBlocks:
        headerInfo = [x for x in b[0].replace('block','').replace('qnumbers','').split()
                      if x != '#']
        if headerInfo[0].replace('-','').replace('+','').isdigit():
            pdg = eval(headerInfo[0])
        else:
            logger.error(f"Error obtaining PDG number from QNUMBERS block:\n {b} \n")
            raise SModelSError()

        if any(p.pdg == pdg for p in BSMList):
            logger.warning("Particle with pdg %i appears multiple times in QNUMBERS blocks" %pdg)
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
            logger.error(f"Error reading quantum numbers from block: \n {b} \n")
            continue
        if any(x not in numbers for x in [1,2,3]):
            logger.error(f"Missing quantum numbers in block:\n {b}\n")
            continue

        # Ignore SM blocks:
        if abs(pdg) in SMpdgs:
            continue

        # If it is not a SM particle, assume it is a BSM particle
        newParticle = Particle(isSM=False, label=label, pdg=pdg,
                               eCharge=numbers[1]/3.,
                               colordim=numbers[3],
                               spin=(numbers[2]-1.)/2.)

        BSMList.append(newParticle)
        if numbers[4]:  # Particle is not its own anti-particle
            newParticleC = newParticle.chargeConjugate()
            if any(p.pdg == newParticleC.pdg for p in BSMList):
                continue
            BSMList.append(newParticleC)

    return BSMList

def getParticlesFromLHE(lhefile: str) -> list[Particle]:
    """
    Defines BSM particles from the QNUMBERS blocks in the <slha> block contained in the lhefile.

    :param lhefile: Path to the LHE file containing a <slha> block

    :return: List with Particle objects
    """

    checkDirs = [os.path.join(installDirectory(), "smodels", "share", "models"), installDirectory(),
                os.path.join(installDirectory(), "smodels")]


    filename = lhefile
    #If file does not exist, check if it is in any of the default folders:
    if not os.path.isfile(lhefile):
        for dirPath in checkDirs:
            if os.path.isfile(os.path.join(dirPath, lhefile)):
                filename = os.path.join(dirPath, lhefile)
                break

    if not os.path.isfile(filename):
        logger.error(f"Model file {lhefile} not found.")
        raise SModelSError()

    logger.debug(f"Trying to define BSM particles from LHE input file {filename}")


    import re
    import xml.etree.ElementTree as ET
    with open(filename,'r') as f:
        lhe_text = f.read()

    # Locate the first <slha>...</slha> block
    m = re.search( r'<slha\b[^>]*>(?P<inner>.*?)</slha>', lhe_text,
                    flags=re.DOTALL|re.IGNORECASE )
    if not m:
        raise ValueError("No <slha>...</slha> block found.")
    inner = m.group('inner')
    # Wrap the inner content to parse the mixed text+elements safely
    wrapper_xml = f'<__wrap__>{inner}</__wrap__>'
    # Default parser is fine; we don't need to preserve comments as "tags"
    wrapper = ET.fromstring(wrapper_xml)
    slhaStr = wrapper.text
    if not slhaStr or 'block qnumbers' not in slhaStr.lower():
        raise ValueError(f"No qnumbers block found in {lhefile}.")

    BSMList = getParticlesFromSLHA(slhaStr)

    return BSMList

def getParticlesFromModule(modelFile: str) -> list[Particle]:
    """
    Reads the python model file and retrieves the list of BSM particles (BSMList)

    :param modelFile: Name/path to the python module containing the BSM particle definitions

    :return: a list of Particle objects
    """

    fulldir = os.path.join(installDirectory(), "smodels", "share", "models")
    sys.path.insert(0, installDirectory())
    sys.path.insert(0, os.path.join(installDirectory(), "smodels"))
    sys.path.insert(0, fulldir)
    sys.path.insert(0, ".")

    logger.debug(f"Trying to load model file: {modelFile}")

    fname = modelFile[:]
    if "/" in fname:
        import shutil
        filename = os.path.basename(fname)
        if not os.path.exists ( filename ) or not os.path.samefile ( fname, filename ):
            shutil.copy(fname, filename)
    else:
        filename = fname

    if filename.endswith(".py"):
        importName = filename[:-3]
    else:
        importName = filename

    pM=import_module(importName, package='smodels')
    logger.debug(f"Found model file at {pM.__file__}")
    if filename != fname:
        os.remove(filename)
    BSMList = pM.BSMList

    return BSMList



def load() -> list[Particle]:
    """
    Load the BSM particle definitions from the default model file.
    Tries to load as a Python module first, then as SLHA or LHE.

    :return: list of BSM Particle objects
    """
    from smodels.base.runtime import modelFile

    BSMList = None
    try:
        BSMList = getParticlesFromModule(modelFile)
    #If failed, assume the input is an SLHA or LHE file:
    except (ImportError, AttributeError, ValueError, SModelSError):
        checkDirs = [os.path.join(installDirectory(), "smodels", "share", "models"),
                      installDirectory(),
                      os.path.join(installDirectory(), "smodels")]

        filename = modelFile
        #If file does not exist, check if it is in any of the default folders:
        if not os.path.isfile(modelFile):
            for dirPath in checkDirs:
                if os.path.isfile(os.path.join(dirPath, modelFile)):
                    filename = os.path.join(dirPath, modelFile)
                    break

        if not os.path.isfile(filename):
            logger.error(f"Model file {modelFile} not found.")
            raise SModelSError()

        with open(filename,'r') as f:
            data = f.read().lower()
        if ('<slha>' in data) and ('</slha>' in data):
            try:
                BSMList = getParticlesFromLHE(modelFile)
            except (ImportError, AttributeError, ValueError, SModelSError):
                pass
        elif 'block qnumbers' in data:
            try:
                BSMList = getParticlesFromSLHA(modelFile)
            except (ImportError, AttributeError, ValueError, SModelSError):
                pass

    if BSMList is None:
        logger.error(f"Could not load input model from {modelFile}. The file should be either a python module with particle definitions, a SLHA file with QNUMBERS blocks or a LHE file with a <slha> block.")
        raise SModelSError()

    return BSMList

