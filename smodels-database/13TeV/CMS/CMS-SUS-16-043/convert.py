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
type = types.StringType)
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
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
if args.smodelsPath:
    sys.path.append(os.path.abspath(args.smodelsPath))

sys.path.append(os.path.abspath(utilsPath))
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z



#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-16-043')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-043/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'EWK WH(bb)'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1706.09933'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'JHEP 11 (2017) 029 , http://dx.doi.org/10.1007/JHEP11(2017)029'
info.comment = 'Moriond 2017'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++

TChiWH=dataset.addTxName('TChiWH')
TChiWH.checked=''
TChiWH.constraint="[[['W']],[['higgs']]]"
TChiWH.condition=None
TChiWH.conditionDescription = None
TChiWH.source="CMS"
TChiWH.massConstraint=None

#++++++++++mass plane++++++++++++++++++

TChiWH_1 = TChiWH.addMassPlane(2*[[x,y]])
TChiWH_1.figure='Fig. 6'
TChiWH_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-043/CMS-SUS-16-043_Figure_006-b.png'
TChiWH_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-043/CMS-SUS-16-043_Figure_006-b.root'
TChiWH_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-043/CMS-SUS-16-043_Figure_006-b.root'
TChiWH_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-043/CMS-SUS-16-043_Figure_006-b.root'
TChiWH_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-043_Figure_006-b.root','orig/CMS-SUS-16-043_Figure_006-b.root','orig/CMS-SUS-16-043_Figure_006-b.root','orig/CMS-SUS-16-043_Figure_006-b.root','orig/CMS-SUS-16-043_Figure_006-b.root','orig/CMS-SUS-16-043_Figure_006-b.root','orig/CMS-SUS-16-043_Figure_006-b.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['gExp;1','gExpDn;1','gExpDn;1','gObs;1','gObsDn;1','gObsUp;1','hObsXsec;1'],
                    units=[None,None,None,None,None,None,'pb'])



databaseCreator.create()
