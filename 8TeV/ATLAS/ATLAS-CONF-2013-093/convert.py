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
info = MetaInfoInput('ATLAS-CONF-2013-093')
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-093/'
info.supersededBy = 'ATLAS-SUSY-2013-23'
info.prettyName = '1 lepton + 2 b-jets + Etmiss (mbb = mH)'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiWH = dataset.addTxName('TChiWH')
TChiWH.checked ="VM"
TChiWH.constraint ="[[['higgs']],[['W']]]"
TChiWH.conditionDescription ="None"
TChiWH.condition ="None"
TChiWH.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiWH = TChiWH.addMassPlane(2*[[x, y]])
TChiWH.figure = "fig 10"
TChiWH.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-093/fig_10.png" 
TChiWH.dataUrl = 'Not defined'
TChiWH.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/TChiWH_exc.dat', 'orig/TChiWH.dat'],
                 dataFormats= ['txt', 'txt'])



databaseCreator.create()
