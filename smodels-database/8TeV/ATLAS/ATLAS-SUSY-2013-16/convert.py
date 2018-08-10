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
info = MetaInfoInput('ATLAS-SUSY-2013-16')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-16/'
info.sqrts = '8'
info.lumi = '20.1'
info.prettyName = '0 lepton + 6 (2 b-)jets + Etmiss'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1406.1122'
info.contact = ''
info.publication = 'http://link.springer.com/article/10.1007%2FJHEP09%282014%29015'
info.supersedes = 'ATLAS-CONF-2013-024'
info.comment =''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked = ""
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription = "None"
T2tt.condition = "None"
T2tt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure ="fig_08"
T2tt_1.figureUrl ="https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-16/fig_08.png"
T2tt_1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1299143'
T2tt_1.histoDataUrl ='http://hepdata.cedar.ac.uk/view/ins1299143'
T2tt_1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1299143'
T2tt_1.exclusionDataUrl ='http://hepdata.cedar.ac.uk/view/ins1299143'
T2tt_1.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2tt_Exclusion.txt', 'orig/T2tt_UL.dat'],
                 dataFormats= ['txt', 'txt'])



databaseCreator.create()
