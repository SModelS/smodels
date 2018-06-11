#!/usr/bin/env python3

"""
.. module:: particlesLoader
   :synopsis: Loads the file Defining the list of R-even and R-odd particles to be used.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Matthias Wolf <matthias.wolf@wot.at>

   (Parameter descriptions taken from Andre Lessa's module particles)
   :parameter rOdd: dictionary with PDG codes for the rOdd (Z2-odd) particles
   and their respective labels
   :parameter rEven: dictionary with PDG codes for the rEven (Z2-eveb)
   particles and their respective labels
   :parameter ptcDic: dictionary with inclusive labels to help defining group
   of particles in the analysis database
   
"""

def load ():
    import os, sys
    from smodels.tools.runtime import modelFile
    from smodels.tools.smodelsLogging import logger
    from smodels.installation import installDirectory
    fulldir = os.path.join(installDirectory(),"inputFiles","models")
    sys.path.insert(0,fulldir)
    sys.path.insert(0,".")

    logger.debug ( "Trying to load model file: %s" % modelFile )

    if "/" in modelFile:
        import shutil
        filename=os.path.basename(modelFile)
        shutil.copy ( modelFile, filename )
        modelFile=filename

    if modelFile.endswith(".py"):
        modelFile=modelFile[:-3]

    from importlib import import_module
    try:
        pM=import_module (modelFile, package='smodels')
        logger.debug ( "Found model file at %s" % pM.__file__ )
        return pM
    except ModuleNotFoundError as e:
        logger.error ( "Model file %s not found." % modelFile )
        sys.exit()

pM = load()
        
rOdd = pM.rOdd
rEven = pM.rEven
qNumbers = pM.qNumbers
