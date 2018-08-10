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
info = MetaInfoInput('CMS-SUS-17-001')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-001/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'Stop search in dilepton + jets + Etmiss final state'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1711.00752'
info.contact = ''
info.publication = 'to appear in Phys. Rev. D '
info.comment = 'Moriond 2017. A long cascade decay not implemented: stop -> b chargino, chargino-> slepton v, slepton -> l lsp with chargino fixed halfway between the stop and lsp and 3 mass planes provided for the intermediate slepton at x = 0.05,0.5 and 0.95 the mass of the chargino. '


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)


#+++++txName block +++++++++++++++++

T2tt=dataset.addTxName('T2tt')
T2tt.checked=''
T2tt.constraint="[[['t']],[['t']]]"
T2tt.condition=None
T2tt.conditionDescription = None
T2tt.source="CMS"
T2tt.massConstraint=None

T2ttoff=dataset.addTxName('T2ttoff')
T2ttoff.checked=''
T2ttoff.constraint="[[['b','W']],[['b','W']]]"
T2ttoff.condition=None
T2ttoff.conditionDescription = None
T2ttoff.source="CMS"
T2ttoff.massConstraint=[['80<dm<169'],['80<dm<169']]

#++++++next mass plane block+++++++++

T2tt_1 = T2tt.addMassPlane(2*[[x,y]])
T2tt_1.figure='Fig. 11-a'
T2tt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-001/CMS-SUS-17-001_Figure_011-a.png'
T2tt_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-001/CMS-SUS-17-001_Figure_011-a.root'
T2tt_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-001/CMS-SUS-17-001_Figure_011-a.root'
T2tt_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-001/CMS-SUS-17-001_Figure_011-a.root'
T2tt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-17-001_Figure_011-a.root','orig/CMS-SUS-17-001_Figure_011-a.root','orig/CMS-SUS-17-001_Figure_011-a.root','orig/CMS-SUS-17-001_Figure_011-a.root','orig/CMS-SUS-17-001_Figure_011-a.root','orig/CMS-SUS-17-001_Figure_011-a.root','orig/CMS-SUS-17-001_Figure_011-a.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1','cCONT_;1'],indices=[3,5,4,6,8,7,2],units=[None,None,None,None,None,None,'pb'])

T2ttoff.addMassPlane(T2tt_1);
#+++++txName block +++++++++++++++++

T6bbWW=dataset.addTxName('T6bbWW')
T6bbWW.checked=''
T6bbWW.constraint="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.condition=None
T6bbWW.conditionDescription = None
T6bbWW.source="CMS"
T6bbWW.massConstraint=None
#++++++next mass plane block+++++++++

T6bbWW_1 = T6bbWW.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
T6bbWW_1.figure='Fig. 11-b'
T6bbWW_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-001/CMS-SUS-17-001_Figure_011-b.root'
T6bbWW_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-001/CMS-SUS-17-001_Figure_011-b.root'
T6bbWW_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-001/CMS-SUS-17-001_Figure_011-b.root'
T6bbWW_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-001/CMS-SUS-17-001_Figure_011-b.root'
T6bbWW_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-17-001_Figure_011-b.root','orig/CMS-SUS-17-001_Figure_011-b.root','orig/CMS-SUS-17-001_Figure_011-b.root','orig/CMS-SUS-17-001_Figure_011-b.root','orig/CMS-SUS-17-001_Figure_011-b.root','orig/CMS-SUS-17-001_Figure_011-b.root','orig/CMS-SUS-17-001_Figure_011-b.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['cCONT_asdf;1','cCONT_asdf;1','cCONT_asdf;1','cCONT_asdf;1','cCONT_asdf;1','cCONT_asdf;1','cCONT_asdf;1'],indices=[3,5,4,6,8,7,2],units=[None,None,None,None,None,None,'pb'])


databaseCreator.create()
