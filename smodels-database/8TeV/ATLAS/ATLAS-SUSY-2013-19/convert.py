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
info = MetaInfoInput('ATLAS-SUSY-2013-19')
info.comment = 'T2tt UL are from DF channel only, no combined UL map available'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication = 'http://link.springer.com/article/10.1007/JHEP06(2014)124'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-19/'
info.arxiv = 'http://arxiv.org/abs/1403.4853'
info.prettyName = '2 OS leptons + (b)jets + Etmiss (leptonic/hadronic m_T2)'
info.supersedes = 'ATLAS-CONF-2013-048'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint ="[[['t+']],[['t-']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2tt = T2tt.addMassPlane(2*[[x, y]])
T2tt.dataUrl = "http://hepdata.cedar.ac.uk/view/ins1286444/d72"
T2tt.histoDataUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-19/figaux_10a.png"
T2tt.figure = "fig 10a"
T2tt.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-19/figaux_10a.png"
T2tt.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusionline_T2tt_DF.txt', 'orig/T2tt.txt'],
                 dataFormats= ['txt', 'txt'])

#+++++++ next txName block ++++++++++++++
T2bbWW = dataset.addTxName('T2bbWW')
T2bbWW.constraint ="[[['b','W+']],[['b','W-']]]"
T2bbWW.conditionDescription ="None"
T2bbWW.condition ="None"
T2bbWW.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2bbWW = T2bbWW.addMassPlane(2*[[x, y]])
T2bbWW.figure = 'Fig.(aux) 3e'
T2bbWW.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-19/figaux_03e.png'
T2bbWW.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1286444/d42'
T2bbWW.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusionline_T2bbWW.txt', 'orig/T2bbWW.txt'],
                 dataFormats= ['txt', 'txt'])

#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.checked ="VM"
T6bbWW.constraint ="[[['b'],['W+']],[['b'],['W-']]]"
T6bbWW.conditionDescription ="None"
T6bbWW.condition ="None"
T6bbWW.source = "ATLAS"
T6bbWW.massConstraint = None
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.constraint ="22*([[['b'],['l+','nu']],[['b'],['l-','nu']]])"
T6bbWWoff.conditionDescription="[[['b'],['l+','nu']],[['b'],['l-','nu']]] > 2*[[['b'],['e+','nu']],[['b'],['e-','nu']]]"
T6bbWWoff.condition="Cgtr([[['b'],['l+','nu']],[['b'],['l-','nu']]],2*[[['b'],['e+','nu']],[['b'],['e-','nu']]])"
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6bbWWLSP001 = T6bbWW.addMassPlane(2*[[x, y, 1.0]])
T6bbWWLSP001.figure = 'Fig.(aux) 3a'
T6bbWWLSP001.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-19/figaux_3a.png'
T6bbWWLSP001.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1286444/d30'
T6bbWWLSP001.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusionline_T6bbWWLSP001.txt', 'orig/T6bbWWLSP001.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWLSP001)
#+++++++ next mass plane block ++++++++++++++
T6bbWWD010 = T6bbWW.addMassPlane(2*[[x, x-10.0, y]])
T6bbWWD010.figure = "fig(aux) 3b"
T6bbWWD010.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-19/figaux_3b.png'
T6bbWWD010.dataUrl = 'Not defined'
T6bbWWD010.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusionline_T6bbWWD10.txt', 'orig/T6bbWWD010.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWD010)
#+++++++ next mass plane block ++++++++++++++
T6bbWWM1300 = T6bbWW.addMassPlane(2*[[300.0, x, y]])
T6bbWWM1300.figure = 'Fig.(aux) 3c'
T6bbWWM1300.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-19/fig_16.png'
T6bbWWM1300.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1286444/d36'
T6bbWWM1300.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusionline_T6bbWWM1300.txt', 'orig/T6bbWWM1300.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWM1300)
#+++++++ next mass plane block ++++++++++++++
T6bbWWC106 = T6bbWW.addMassPlane(2*[[x, 106.0, y]])
T6bbWWC106.figure = 'Fig.(aux) 3f'
T6bbWWC106.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-19/fig_20.png'
T6bbWWC106.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1286444/d68'
T6bbWWC106.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusionline_T6bbWWC106.txt', 'orig/T6bbWWC106.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWC106)
#+++++++ next mass plane block ++++++++++++++
T6bbWWx200 = T6bbWW.addMassPlane(2*[[x, y*2.0, y]])
T6bbWWx200.figure = 'Fig.(aux) 3d'
T6bbWWx200.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-19/fig_17.png'
T6bbWWx200.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1286444/d39'
T6bbWWx200.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusionline_T6bbWWx200.txt', 'orig/T6bbWWx200.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWx200)



databaseCreator.create()
