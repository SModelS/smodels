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
info = MetaInfoInput('ATLAS-SUSY-2015-02')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2015-02/'
info.sqrts = 13
info.lumi = 3.2
info.prettyName = 'Top 1l'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1606.03903'
info.publication ='http://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.052009'
info.comment = 'Very weird exclusion due to best CLs selection, and fluctuation in the Obs data. Digitised by FA'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription = None
T2tt.condition = None
T2tt.source = "ATLAS"
T2tt.massConstraint = None
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint =  "[[['b','W']],[['b','W']]]"
T2ttoff.conditionDescription = None
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0'], ['dm <= 169.0']]
T2ttoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure = 'Aux.Fig.5b'
T2tt_1.figureUrl ='http://hepdata.cedar.ac.uk/resource/1469069/figAuxiliaryFigure5b.png'
T2tt_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1469069/all'
T2tt_1.histoDataUrl = 'http://hepdata.cedar.ac.uk/view/ins1469069/d24/input'
T2tt_1.dataUrl ='http://hepdata.cedar.ac.uk/view/ins1469069/d24/input'
T2tt_1.exclusionDataUrl =  'http://hepdata.cedar.ac.uk/resource/1469069/figFigure8b.png'
T2tt_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2tt_ATLAS_SUSY_2015_02_Exp.txt', 'orig/T2tt_ATLAS_SUSY_2015_02_Obs.txt', 'orig/T2tt_ATLAS_SUSY_2015_02_Obs_UL.dat'],
                 dataFormats= ['txt', 'txt', 'txt'])
T2ttoff.addMassPlane(T2tt_1)



databaseCreator.create()
