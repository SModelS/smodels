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
info = MetaInfoInput('ATLAS-SUSY-2015-01')
info.url = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2015-01/"
info.sqrts = 13
info.lumi = 3.2
info.prettyName = "2b"
info.private = False
info.arxiv =  'https://arxiv.org/abs/1606.08772v2'
info.contact = 'ATLAS collaboration'
info.publication ='http://link.springer.com/article/10.1140/epjc/s10052-016-4382-4'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2bb = dataset.addTxName('T2bb')
T2bb.checked = ''
T2bb.constraint ="[[['b']],[['b']]]"
T2bb.conditionDescription = None
T2bb.condition = None
T2bb.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2bb_1 = T2bb.addMassPlane(2*[[x, y]])
T2bb_1.figure = 'Fig.4'
T2bb_1.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2015-01/fig_04.png'
T2bb_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1472822/d4/input'
T2bb_1.histoDataUrl = 'http://hepdata.cedar.ac.uk/view/ins1472822/d4/input'
T2bb_1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1472822/d4/input'
T2bb_1.exclusionDataUrl = 'http://hepdata.cedar.ac.uk/view/ins1472822/d4/input'
T2bb_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2bb_Exp_Excl.dat', 'orig/T2bb_Obs_Excl.dat', 'orig/T2bb_Obs_UL_fixed.dat'],
                 dataFormats= ['txt', 'txt', 'txt'])



databaseCreator.create()
