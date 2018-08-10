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
info = MetaInfoInput('ATLAS-SUSY-2013-23')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-23/'
info.sqrts = 8
info.lumi = 20.3
info.prettyName = '1 lepton + 2 b-jets (or 2photons) + Etmiss (mbb = mH)'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1501.07110'
info.publication = 'http://link.springer.com/article/10.1140/epjc/s10052-015-3408-7'
info.comment = "Used combined limit"
info.supersedes = 'ATLAS-CONF-2013-093'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiWH = dataset.addTxName('TChiWH')
TChiWH.checked =''
TChiWH.constraint = "[[['higgs']],[['W']]]"
TChiWH.conditionDescription = None
TChiWH.condition = None
TChiWH.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiWH_1 = TChiWH.addMassPlane(2*[[x, y]])
TChiWH_1.figure = 'Fig. Aux1'
TChiWH_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-23/fig_08d.png'
TChiWH_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1341609'
TChiWH_1.histoDataUrl ='http://hepdata.cedar.ac.uk/view/ins1341609'
TChiWH_1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1341609/d34/input'
TChiWH_1.exclusionDataUrl = 'http://hepdata.cedar.ac.uk/view/ins1341609'
TChiWH_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/TChiWH-ATLAS-SUSY-2013-23-Exp_Excl.dat', 'orig/TChiWH-ATLAS-SUSY-2013-23-Obs_Excl.dat', 'orig/TChiWH-ATLAS-SUSY-2013-23-UL.dat'],
                 dataFormats= ['txt', 'txt', 'txt'])



databaseCreator.create()
