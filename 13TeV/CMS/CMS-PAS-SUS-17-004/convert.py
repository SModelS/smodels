#!/usr/bin/env python3

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
info = MetaInfoInput('CMS-PAS-SUS-17-004')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-17-004/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'Multilepton EWK searches'
info.private = False
info.arxiv = ''
info.contact = ''
info.publication = ''
info.comment = 'TChiHH, TChiZZ mass planes not provided. TChiWH/Z not implemented as BR(Chi20->H Chi10) = BR(Chi20->Chi10 Z)=0.5 provided for a single mass plane.'


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

#offshell txName block

TChiWZoff=dataset.addTxName('TChiWZoff')
TChiWZoff.checked=''
TChiWZoff.constraint="71.*([[['mu+','mu-']],[['l','nu']]] + [[['e+','e-']],[['l','nu']]])"
TChiWZoff.condition = "cGtr([[['mu+','mu-']],[['l','nu']]],[[['e+','e-']],[['l','nu']]])"
TChiWZoff.massConstraint = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription=None
TChiWZoff.source="CMS"

#++++++next mass plane block+++++++++

TChiWZ_1 = TChiWZ.addMassPlane(2*[[x,y]])
TChiWZ_1.figure='Fig. 8-a'
TChiWZ_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/PAS-SUS-17-004/CMS-PAS-SUS-17-004_Figure_008-a.png'
TChiWZ_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/PAS-SUS-17-004/CMS-PAS-SUS-17-004_Figure_008-a.root'
TChiWZ_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/PAS-SUS-17-004/CMS-PAS-SUS-17-004_Figure_008-a.root'
TChiWZ_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/PAS-SUS-17-004/CMS-PAS-SUS-17-004_Figure_008-a.root'
TChiWZ_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-PAS-SUS-17-004_Figure_008-a.root','orig/CMS-PAS-SUS-17-004_Figure_008-a.root','orig/CMS-PAS-SUS-17-004_Figure_008-a.root','orig/CMS-PAS-SUS-17-004_Figure_008-a.root','orig/CMS-PAS-SUS-17-004_Figure_008-a.root','orig/CMS-PAS-SUS-17-004_Figure_008-a.root','orig/CMS-PAS-SUS-17-004_Figure_008-a.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['TChiWZ;1','TChiWZ;1','TChiWZ;1','TChiWZ;1','TChiWZ;1','TChiWZ;1','TChiWZ;1'],
                    indices= [4, 6, 5, 7, 9, 8, 2],units=[None,None,None,None,None,None,'pb'])

TChiWZoff.addMassPlane(TChiWZ_1)

#++++++next txName block+++++++++++++++

TChiWH=dataset.addTxName('TChiWH')
TChiWH.checked=''
TChiWH.constraint="[[['W']],[['higgs']]]"
TChiWH.condition=None
TChiWH.conditionDescription=None
TChiWH.source="CMS"


#++++++next mass plane block++++++++
TChiWH_1 = TChiWH.addMassPlane(2*[[x,y]])
TChiWH_1.figure='Fig. 8-b'
TChiWH_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/PAS-SUS-17-004/CMS-PAS-SUS-17-004_Figure_008-b.png'
TChiWH_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/PAS-SUS-17-004/CMS-PAS-SUS-17-004_Figure_008-b.root'
TChiWH_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/PAS-SUS-17-004/CMS-PAS-SUS-17-004_Figure_008-b.root'
TChiWH_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/PAS-SUS-17-004/CMS-PAS-SUS-17-004_Figure_008-c.root'
TChiWH_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-PAS-SUS-17-004_Figure_008-b.root','orig/CMS-PAS-SUS-17-004_Figure_008-b.root','orig/CMS-PAS-SUS-17-004_Figure_008-b.root','orig/CMS-PAS-SUS-17-004_Figure_008-b.root','orig/CMS-PAS-SUS-17-004_Figure_008-b.root','orig/CMS-PAS-SUS-17-004_Figure_008-b.root','orig/CMS-PAS-SUS-17-004_Figure_008-b.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1'],
                    indices= [4, 6, 5, 7, 9, 8, 2],units=[None,None,None,None,None,None,'pb'])

databaseCreator.create()

