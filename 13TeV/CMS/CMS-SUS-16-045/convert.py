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
info = MetaInfoInput('CMS-SUS-16-045')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-045/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'Sbottom searches in bHbH and H -> photon photon final states'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1709.00384'
info.contact = ''
info.publication = 'Submitted to Phys. Lett. B '
info.comment = 'TChiWH topology has a wrong root file. The exclusion plot 5b does not match the exclusion lines in the corresponding root file. Renamed root files to CMS-SUS-16-045_Figure_005-1.root and CMS-SUS-16-045_Figure_005-2.root '


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++

T6bbHH=dataset.addTxName('T6bbHH')
T6bbHH.checked=''
T6bbHH.constraint="[[['b'],['higgs']],[['b'],['higgs']]]"
T6bbHH.condition=None
T6bbHH.conditionDescription = None
T6bbHH.source="CMS"
T6bbHH.massConstraint = None

#++++++++++mass plane++++++++++++++++++

T6bbHH_1 = T6bbHH.addMassPlane(2*[[x,y+130.0,y]])
T6bbHH_1.figure='Fig. 5-a'
T6bbHH_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-045/CMS-SUS-16-045_Figure_005-a.png'
T6bbHH_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-045/CMS-SUS-16-045_Figure_005.root'
T6bbHH_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-045/CMS-SUS-16-045_Figure_005.root'
T6bbHH_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-045/CMS-SUS-16-045_Figure_005.root'
T6bbHH_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-045_Figure_005-1.root','orig/CMS-SUS-16-045_Figure_005-1.root','orig/CMS-SUS-16-045_Figure_005-1.root','orig/CMS-SUS-16-045_Figure_005-1.root','orig/CMS-SUS-16-045_Figure_005-1.root','orig/CMS-SUS-16-045_Figure_005-1.root','orig/CMS-SUS-16-045_Figure_005-1.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['expL;1','expLD;1','expLU;1','obsL;1','obsLD;1','obsLU;1','xsecLimit;1'],
                    units=[None,None,None,None,None,None,'pb'])


#+++++txName block +++++++++++++++++

TChiWH=dataset.addTxName('TChiWH')
TChiWH.checked=''
TChiWH.constraint="[[['W']],[['higgs']]]"
TChiWH.condition=None
TChiWH.conditionDescription = None
TChiWH.source="CMS"
TChiWH.massConstraint = None

#++++++++++mass plane++++++++++++++++++

TChiWH_1 = TChiWH.addMassPlane(2*[[x,y]])
TChiWH_1.figure='Fig. 5-b'
TChiWH_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-045/CMS-SUS-16-045_Figure_005-b.png'
TChiWH_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-045/CMS-SUS-16-045_Figure_005-b.root'
TChiWH_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-045/CMS-SUS-16-045_Figure_005-b.root'
TChiWH_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-045/CMS-SUS-16-045_Figure_005-b.root'
TChiWH_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-045_Figure_005-b.root','orig/CMS-SUS-16-045_Figure_005-b.root','orig/CMS-SUS-16-045_Figure_005-b.root','orig/CMS-SUS-16-045_Figure_005-b.root','orig/CMS-SUS-16-045_Figure_005-b.root','orig/CMS-SUS-16-045_Figure_005-b.root','orig/CMS-SUS-16-045_Figure_005-b.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['expL;1','expLD;1','expLU;1','obsL;1','obsLD;1','obsLU;1','xsecLimit;1'],
                    units=[None,None,None,None,None,None,'pb'])


databaseCreator.create()
