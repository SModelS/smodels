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
info = MetaInfoInput('CMS-SUS-16-034')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = '2 OSSF leptons'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1709.08908'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Submitted to J. High Energy Phys. '
info.comment = 'Moriond 2017. Fig. 10 with long cascade decay via chargino and slepton not implemented because of assumtions on masses and BRs.'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++
TChiWZ=dataset.addTxName('TChiWZ')
TChiWZ.checked=''
TChiWZ.constraint="[[['W']],[['Z']]]"
TChiWZ.condition=None
TChiWZ.conditionDescription = None
TChiWZ.source="CMS"
TChiWZ.massConstraint=None

#++++++ mass plane block+++++++++

TChiWZ_1 = TChiWZ.addMassPlane(2*[[x,y]])
TChiWZ_1.figure='Fig. 8'
TChiWZ_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_008.png'
TChiWZ_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_008.root'
TChiWZ_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_008.root'
TChiWZ_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_008.root'
TChiWZ_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-034_Figure_008.root','orig/CMS-SUS-16-034_Figure_008.root','orig/CMS-SUS-16-034_Figure_008.root','orig/CMS-SUS-16-034_Figure_008.root','orig/CMS-SUS-16-034_Figure_008.root','orig/CMS-SUS-16-034_Figure_008.root','orig/CMS-SUS-16-034_Figure_008.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['gExp;1','gExpDn;1','gExpUp;1','gObs;1','gObsDn;1','gObsUp;1','hObsXsec;1'],
                    units=[None,None,None,None,None,None,'pb'])

#+++++txName block +++++++++++++++++
T5ZZ=dataset.addTxName('T5ZZ')
T5ZZ.checked=''
T5ZZ.constraint="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.condition=None
T5ZZ.conditionDescription = None
T5ZZ.source="CMS"
T5ZZ.massConstraint=None

#++++++ mass plane block+++++++++

T5ZZ_1 = T5ZZ.addMassPlane(2*[[x,y,1.0]])
T5ZZ_1.figure='Fig. 7'
T5ZZ_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_007.png'
T5ZZ_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_007.root'
T5ZZ_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_007.root'
T5ZZ_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_007.root'
T5ZZ_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['gExp;1','gExpDn2;1','gExpUp2;1','gObs;1','gObsDn;1','gObsUp;1','hObsXsec;1'],
                    units=[None,None,None,None,None,None,'pb'])

T5ZZ_2 = T5ZZ.addMassPlane(2*[[x,y,0.0]])
T5ZZ_2.figure='Fig. 7'
T5ZZ_2.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_007.png'
T5ZZ_2.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_007.root'
T5ZZ_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_007.root'
T5ZZ_2.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-034/CMS-SUS-16-034_Figure_007.root'
T5ZZ_2.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root','orig/CMS-SUS-16-034_Figure_007.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['gExp;1','gExpDn2;1','gExpUp2;1','gObs;1','gObsDn;1','gObsUp;1','hObsXsec;1'],
                    units=[None,None,None,None,None,None,'pb'])
databaseCreator.create()
