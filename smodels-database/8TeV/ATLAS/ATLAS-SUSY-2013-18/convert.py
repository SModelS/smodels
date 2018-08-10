#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse

argparser = argparse.ArgumentParser(description =  
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath', 
help = 'path to the package smodels_utils',\
type = str )
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = str )
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
info = MetaInfoInput('ATLAS-SUSY-2013-18')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/'
info.sqrts = '8.0*TeV'
info.lumi = 20.1
info.prettyName = '0-1 leptons + >= 3 b-jets + Etmiss'
info.private =False
info.arxiv = 'http://arxiv.org/abs/1407.0600'
info.contact =''
info.publication ='http://link.springer.com/article/10.1007/JHEP10(2014)024'
info.comment ='Used 1lep+3b UL result, since they give the two analyses UL separately and no combination.T1btbt result lost wrt ATLAS-CONF-2013-061'
info.supersedes ='ATLAS-CONF-2013-061'
info.supersededBy =''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked = ''
T1bbbb.constraint = "[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane(2*[[x, y]])
T1bbbb_1.figure = "Figure 13a"
T1bbbb_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_03.png'
T1bbbb_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1304457/all'
T1bbbb_1.histoDataUrl = 'http://hepdata.cedar.ac.uk/view/ins1304457/all'
T1bbbb_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1304457/d28/input'
T1bbbb_1.exclusionDataUrl = 'http://hepdata.cedar.ac.uk/view/ins1304457/all'
T1bbbb_1.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T1bbbb_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'orig/T1bbbb_ATLAS-SUSY-2013-18_UL.dat'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked = ''
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure ='Figure 12b aux'
T1tttt_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_12b.png'
T1tttt_1.dataUrl = ''
T1tttt_1.histoDataUrl = 'http://hepdata.cedar.ac.uk/view/ins1304457/all'
T1tttt_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1304457/d41/input'
T1tttt_1.exclusionDataUrl = 'http://hepdata.cedar.ac.uk/view/ins1304457/all'
T1tttt_1.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T1tttt_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'orig/T1tttt_ATLAS-SUSY-2013-18_UL.dat'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])



databaseCreator.create()
