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
info = MetaInfoInput('CMS-SUS-16-046')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'Photon + Etmiss'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1711.08008'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Submitted to Phys. Lett. B '
info.comment = 'Moriond 2017. Several UL maps in the context of gauge-mediated SSB. Implemented T5gg and T6gg; others with mixed BRs are not useable. NB results only for high gluino mass.'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++

T5gg=dataset.addTxName('T5gg')
T5gg.checked=''
T5gg.constraint="[[['jet','jet'],['photon']],[['jet','jet'],['photon']]]"
T5gg.condition=None
T5gg.conditionDescription = None
T5gg.source="CMS"


#++++++next mass plane block+++++++++

T5gg_1 = T5gg.addMassPlane(2*[[x,y,1.0]])
T5gg_1.figure='Fig. 9-a'
T5gg_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_009-a.png'
T5gg_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_009-a.root'
T5gg_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_009-a.root'
T5gg_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_009-a.root'
T5gg_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['exp;1','exp-1;1','exp+1;1','obs;1','obs-1;1','obs+1;1','h_obs_xs;1'],units=[None,None,None,None,None,None,'pb'])

T5gg_2 = T5gg.addMassPlane(2*[[x,y,0.0]])
T5gg_2.figure='Fig. 9-a'
T5gg_2.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_009-a.png'
T5gg_2.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_009-a.root'
T5gg_2.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_009-a.root'
T5gg_2.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_009-a.root'
T5gg_2.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root','orig/CMS-SUS-16-046_Figure_009-a.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['exp;1','exp-1;1','exp+1;1','obs;1','obs-1;1','obs+1;1','h_obs_xs;1'],units=[None,None,None,None,None,None,'pb'])

#+++++txName block +++++++++++++++++

T6gg=dataset.addTxName('T6gg')
T6gg.checked=''
T6gg.constraint="[[['jet'],['photon']],[['jet'],['photon']]]"
T6gg.condition=None
T6gg.conditionDescription = None
T6gg.source="CMS"


#++++++next mass plane block+++++++++

T6gg_1 = T6gg.addMassPlane(2*[[x,y,1.0]])
T6gg_1.figure='Fig. 10-a'
T6gg_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_010-a.png'
T6gg_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_010-a.root'
T6gg_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_010-a.root'
T6gg_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_010-a.root'
T6gg_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['exp;1','exp-1;1','exp+1;1','obs;1','obs-1;1','obs+1;1','h_obs_xs;1'],units=[None,None,None,None,None,None,'pb'])

T6gg_2 = T6gg.addMassPlane(2*[[x,y,0.0]])
T6gg_2.figure='Fig. 10-a'
T6gg_2.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_010-a.png'
T6gg_2.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_010-a.root'
T6gg_2.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_010-a.root'
T6gg_2.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-046/CMS-SUS-16-046_Figure_010-a.root'
T6gg_2.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                  dataFiles=['orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root','orig/CMS-SUS-16-046_Figure_010-a.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['exp;1','exp-1;1','exp+1;1','obs;1','obs-1;1','obs+1;1','h_obs_xs;1'],units=[None,None,None,None,None,None,'pb'])


databaseCreator.create()
