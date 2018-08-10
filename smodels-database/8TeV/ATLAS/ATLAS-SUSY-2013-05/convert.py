#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse
import types

argparser = argparse.ArgumentParser(description =  
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath', 
help = 'path to the package smodels_utils',\
type = str)
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = str)
args = argparser.parse_args()

if args.utilsPath:
    utilsPath = args.utilsPath
else:
    databaseRoot = '../../../'
    sys.path.append(os.path.abspath(databaseRoot))
    from utilsPath import utilsPath
    utilsPath = databaseRoot + utilsPath
if args.smodelsPath:
    sys.path.append(os.path.abspath(args.smodelsPath))

sys.path.append(os.path.abspath(utilsPath))
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z



#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2013-05')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/'
info.sqrts = '8.0*TeV'
info.lumi = 20.1
info.prettyName = '0 leptons + 2 b-jets + Etmiss'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1308.2631'
info.publication = 'http://link.springer.com/article/10.1007/JHEP10%282013%29189'
info.comment = 'upper limits for T6bbWWC150 are not public'
info.supersedes = 'ATLAS_CONF_2013_001;ATLAS_CONF_2013_053'
info.implementedBy = 'MT'
info.comment = 'upper limits for T6bbWWC150 are not public'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2bb = dataset.addTxName('T2bb')
T2bb.checked = 'VM'
T2bb.constraint = "[[['b']],[['b']]]"
T2bb.conditionDescription = None
T2bb.condition = None
T2bb.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2bb_1 = T2bb.addMassPlane(2*[[x, y]])
T2bb_1.figure = 'Fig.(aux) 4'
T2bb_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/figaux_04.png'
T2bb_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1247462/d32'
T2bb_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/T2bb_excExpected.dat', 'orig/T2bb_excExpectedMinusSigma.dat', 'orig/T2bb_excExpectedPlusSigma.dat', 'orig/T2bb_exc.dat', 'orig/T2bb_excMinusSigma.dat', 'orig/T2bb_excPlusSigma.dat', 'orig/T2bb_2014-09-22.dat'],
                 dataFormats= ['txt', 'txt', 'txt', 'txt', 'txt', 'txt', 'txt'])

#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.checked = 'VM'
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription = None
T6bbWW.condition = None
T6bbWW.source = "ATLAS"
T6bbWW.massConstraint = None
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.checked = 'VM'
T6bbWWoff.constraint = "[[['b'],['L','nu']],[['b'],['L','nu']]] + [[['b'],['L','nu']],[['b'],['jet','jet']]] + [[['b'],['jet','jet']],[['b'],['jet','jet']]]"
T6bbWWoff.conditionDescription = None
T6bbWWoff.condition = None
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6bbWWD005 = T6bbWWoff.addMassPlane(2*[[x, y+5., y]]) #This is a purely off-shell plane, do not add it to the on-shell txname
T6bbWWD005.figure = 'Fig.(aux) 8a'
T6bbWWD005.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/figaux_08a.png'
T6bbWWD005.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1247462/d37'
T6bbWWD005.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/T6bbWWoffD005_excExpected.dat', 'orig/T6bbWWoffD005_excExpectedMinusSigma.dat', 'orig/T6bbWWoffD005_excExpectedPlusSigma.dat', 'orig/T6bbWWoffD005_exc.dat', 'orig/T6bbWWoffD005_excMinusSigma.dat', 'orig/T6bbWWoffD005_excPlusSigma.dat', 'orig/T6bbWWoffD005_2014-09-22.dat'],
                 dataFormats= ['txt', 'txt', 'txt', 'txt', 'txt', 'txt', 'txt'])
#+++++++ next mass plane block ++++++++++++++
T6bbWWD020 = T6bbWWoff.addMassPlane(2*[[x, y+20., y]]) #This is a purely off-shell plane, do not add it to the on-shell txname
T6bbWWD020.figure = 'Fig.(aux) 8b'
T6bbWWD020.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/figaux_08b.png'
T6bbWWD020.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1247462/d38'
T6bbWWD020.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/T6bbWWoffD020_excExpected.dat', 'orig/T6bbWWoffD020_excExpectedMinusSigma.dat', 'orig/T6bbWWoffD020_excExpectedPlusSigma.dat', 'orig/T6bbWWoffD020_exc.dat', 'orig/T6bbWWoffD020_excMinusSigma.dat', 'orig/T6bbWWoffD020_excPlusSigma.dat', 'orig/T6bbWWoffD020_2014-09-22.dat'],
                 dataFormats= ['txt', 'txt', 'txt', 'txt', 'txt', 'txt', 'txt'])
#+++++++ next mass plane block ++++++++++++++
T6bbWWM1300 = T6bbWW.addMassPlane(2*[[300., x, y]])
T6bbWWM1300.figure = 'Fig.(aux) 9'
T6bbWWM1300.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-05/figaux_09.png'
T6bbWWM1300.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1247462/d39'
T6bbWWM1300.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/T6bbWWM1300_excExpected.dat', 'orig/T6bbWWM1300_excExpectedMinusSigma.dat', 'orig/T6bbWWM1300_excExpectedPlusSigma.dat', 'orig/T6bbWWM1300_exc.dat', 'orig/T6bbWWM1300_excMinusSigma.dat', 'orig/T6bbWWM1300_excPlusSigma.dat', 'orig/T6bbWWM1300_2014-09-22.dat'],
                 dataFormats= ['txt', 'txt', 'txt', 'txt', 'txt', 'txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWM1300)



databaseCreator.create()
