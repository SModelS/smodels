#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: uesed to create info.txt, txname.txt, twiki.txt and sms.py.

.. moduleauthor:: Michael Traub <michael.traub@gmx.at>

"""   
import sys
import os
import argparse
import types

argparser = argparse.ArgumentParser(description = \
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath', \
help = 'path to the package smodels_utils',\
type = types.StringType)
args = argparser.parse_args()
    
if args.utilsPath:
    utilsPath = args.utilsPath
else:
    databaseRoot = '../../../'
    sys.path.append(os.path.abspath(databaseRoot))
    from utilsPath import utilsPath 
    utilsPath = databaseRoot + utilsPath

sys.path.append(os.path.abspath(utilsPath))
from smodels_utils.dataPreparation.inputObjects import TxNameInput, MetaInfoInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.origPlotObjects import x, y


info = MetaInfoInput('ATLAS-SUSY-2013-05')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/'
info.sqrts = '8.0*TeV'
info.lumi = 20.1
info.prettyname = 'ATLAS 2b'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1308.2631'
info.publication = 'http://link.springer.com/article/10.1007/JHEP10%282013%29189'
info.comment = 'upper limits for T6bbWWC150 are not public'
info.supersedes = 'ATLAS_CONF_2013_001;ATLAS_CONF_2013_053'
info.implemented_by = 'MT'
info.comment = 'upper limits for T6bbWWC150 are not public'

#+++++++++++ add new txName +++++++++++++++++

T2bb = TxNameInput('T2bb')
T2bb.on.checked = 'VM'
T2bb.off.checked = 'VM'
#-----constraint,condition,....------------------------------
T2bb.on.constraint = "[[['b']],[['b']]]"
T2bb.on.condition = None
T2bb.on.fuzzycondition = None

#------ add new massplane-------

T2bb_1 = T2bb.addMassPlane(motherMass = x, lspMass = y)
T2bb_1.figure = 'Fig.(aux) 4'
T2bb_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/figaux_04.png'
#-----limits------------------------------
T2bb_1.obsUpperLimit.setSource('orig/T2bb_2014-09-22.dat','txt')
T2bb_1.obsUpperLimit.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1247462/d32'
#-----exclusions------------------------------
T2bb_1.obsExclusion.setSource('orig/T2bb_exc.dat','txt')
T2bb_1.obsExclusionM1.setSource('orig/T2bb_excMinusSigma.dat','txt')
T2bb_1.obsExclusionP1.setSource('orig/T2bb_excPlusSigma.dat','txt')
T2bb_1.expExclusion.setSource('orig/T2bb_excExpected.dat','txt')
T2bb_1.expExclusionM1.setSource('orig/T2bb_excExpectedMinusSigma.dat','txt')
T2bb_1.expExclusionP1.setSource('orig/T2bb_excExpectedPlusSigma.dat','txt')
#-----exclusions.dataUrl------------------------------
T2bb_1.obsExclusion.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d1'
T2bb_1.obsExclusionM1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d2'
T2bb_1.obsExclusionP1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d3'
T2bb_1.expExclusion.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d4'
T2bb_1.expExclusionM1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d5'
T2bb_1.expExclusionP1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d5'


#+++++++++++ add new txName +++++++++++++++++

T6bbWW = TxNameInput('T6bbWW')
T6bbWW.on.checked = 'VM'
T6bbWW.off.checked = 'VM'
#-----constraint,condition,....------------------------------
T6bbWW.on.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.off.constraint = "[[['b'],['L','nu']],[['b'],['L','nu']]] + [[['b'],['L','nu']],[['b'],['jet','jet']]] + [[['b'],['jet','jet']],[['b'],['jet','jet']]]"
T6bbWW.on.condition = None
T6bbWW.off.condition = None
T6bbWW.on.fuzzycondition = None
T6bbWW.off.fuzzycondition = None

#------ add new massplane-------

T6bbWWD005 = T6bbWW.addMassPlane(motherMass = x, interMass0 = y + 5., lspMass = y)
T6bbWWD005.figure = 'Fig.(aux) 8a'
T6bbWWD005.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/figaux_08a.png'
#-----limits------------------------------
T6bbWWD005.obsUpperLimit.setSource('orig/T6bbWWoffD005_2014-09-22.dat','txt')
T6bbWWD005.obsUpperLimit.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1247462/d37'
#-----exclusions------------------------------
T6bbWWD005.obsExclusion.setSource('orig/T6bbWWoffD005_exc.dat','txt')
T6bbWWD005.obsExclusionM1.setSource('orig/T6bbWWoffD005_excMinusSigma.dat','txt')
T6bbWWD005.obsExclusionP1.setSource('orig/T6bbWWoffD005_excPlusSigma.dat','txt')
T6bbWWD005.expExclusion.setSource('orig/T6bbWWoffD005_excExpected.dat','txt')
T6bbWWD005.expExclusionM1.setSource('orig/T6bbWWoffD005_excExpectedMinusSigma.dat','txt')
T6bbWWD005.expExclusionP1.setSource('orig/T6bbWWoffD005_excExpectedPlusSigma.dat','txt')
#-----exclusions.dataUrl------------------------------
T6bbWWD005.obsExclusion.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d19'
T6bbWWD005.obsExclusionM1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d20'
T6bbWWD005.obsExclusionP1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d21'
T6bbWWD005.expExclusion.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d22'
T6bbWWD005.expExclusionM1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d23'
T6bbWWD005.expExclusionP1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d24'

#------ add new massplane-------

T6bbWWD020 = T6bbWW.addMassPlane(motherMass = x, interMass0 = y + 20., lspMass = y)
T6bbWWD020.figure = 'Fig.(aux) 8b'
T6bbWWD020.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/figaux_08b.png'
#-----limits------------------------------
T6bbWWD020.obsUpperLimit.setSource('orig/T6bbWWoffD020_2014-09-22.dat','txt')
T6bbWWD020.obsUpperLimit.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1247462/d38'
#-----exclusions------------------------------
T6bbWWD020.obsExclusion.setSource('orig/T6bbWWoffD020_exc.dat','txt')
T6bbWWD020.obsExclusionM1.setSource('orig/T6bbWWoffD020_excMinusSigma.dat','txt')
T6bbWWD020.obsExclusionP1.setSource('orig/T6bbWWoffD020_excPlusSigma.dat','txt')
T6bbWWD020.expExclusion.setSource('orig/T6bbWWoffD020_excExpected.dat','txt')
T6bbWWD020.expExclusionM1.setSource('orig/T6bbWWoffD020_excExpectedMinusSigma.dat','txt')
T6bbWWD020.expExclusionP1.setSource('orig/T6bbWWoffD020_excExpectedPlusSigma.dat','txt')
#-----exclusions.dataUrl------------------------------
T6bbWWD020.obsExclusion.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d25'
T6bbWWD020.obsExclusionM1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d26'
T6bbWWD020.obsExclusionP1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d27'
T6bbWWD020.expExclusion.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d28'
T6bbWWD020.expExclusionM1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d29'
T6bbWWD020.expExclusionP1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d30'

#------ add new massplane-------

T6bbWWM1300 = T6bbWW.addMassPlane(motherMass = 300., interMass0 = x, lspMass = y)
T6bbWWM1300.figure = 'Fig.(aux) 9'
T6bbWWM1300.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/figaux_09.png'
#-----limits------------------------------
T6bbWWM1300.obsUpperLimit.setSource('orig/T6bbWWM1300_2014-09-22.dat','txt')
T6bbWWM1300.obsUpperLimit.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1247462/d39'
#-----exclusions------------------------------
T6bbWWM1300.obsExclusion.setSource('orig/T6bbWWM1300_exc.dat','txt')
T6bbWWM1300.obsExclusionM1.setSource('orig/T6bbWWM1300_excMinusSigma.dat','txt')
T6bbWWM1300.obsExclusionP1.setSource('orig/T6bbWWM1300_excPlusSigma.dat','txt')
T6bbWWM1300.expExclusion.setSource('orig/T6bbWWM1300_excExpected.dat','txt')
T6bbWWM1300.expExclusionM1.setSource('orig/T6bbWWM1300_excExpectedMinusSigma.dat','txt')
T6bbWWM1300.expExclusionP1.setSource('orig/T6bbWWM1300_excExpectedPlusSigma.dat','txt')
#-----exclusions.dataUrl------------------------------
T6bbWWM1300.obsExclusion.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d7'
T6bbWWM1300.obsExclusionM1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d8'
T6bbWWM1300.obsExclusionP1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d9'
T6bbWWM1300.expExclusion.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d10'
T6bbWWM1300.expExclusionM1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d11'
T6bbWWM1300.expExclusionP1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d12'

#------ add new massplane-------

T6bbWWC150 = T6bbWW.addMassPlane(motherMass = x, interMass0 = 150., lspMass = y)
T6bbWWC150.figure = 'Fig.(aux) 6b'
T6bbWWC150.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/figaux_06b.png'
#-----limits------------------------------
#T6bbWWC150.obsUpperLimit.setSource('???','txt') NO DATA !!
#T6bbWWC150.obsUpperLimit.dataUrl = '???' NO DATA !!
#-----exclusions------------------------------
T6bbWWC150.obsExclusion.setSource('orig/T6bbWWM2150_exc.dat','txt')
T6bbWWC150.obsExclusionM1.setSource('orig/T6bbWWM2150_excMinusSigma.dat','txt')
T6bbWWC150.obsExclusionP1.setSource('orig/T6bbWWM2150_excPlusSigma.dat','txt')
T6bbWWC150.expExclusion.setSource('orig/T6bbWWM2150_excExpected.dat','txt')
T6bbWWC150.expExclusionM1.setSource('orig/T6bbWWM2150_excExpectedMinusSigma.dat','txt')
T6bbWWC150.expExclusionP1.setSource('orig/T6bbWWM2150_excExpectedPlusSigma.dat','txt')
#-----exclusions.dataUrl------------------------------
T6bbWWC150.obsExclusion.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d13'
T6bbWWC150.obsExclusionM1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d14'
T6bbWWC150.obsExclusionP1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d15'
T6bbWWC150.expExclusion.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d16'
T6bbWWC150.expExclusionM1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d17'
T6bbWWC150.expExclusionP1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1247462/d18'


databaseCreator.create()
