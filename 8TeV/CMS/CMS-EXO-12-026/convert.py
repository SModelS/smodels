#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.dat and the <txname>.dat files.

"""
import sys
import os
import argparse
import types

argparser = argparse.ArgumentParser(description =  
'create info.dat, txname.dat, twiki.dat and sms.py')
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
info = MetaInfoInput('CMS-EXO-12-026')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-12-026/index.html'
info.sqrts = 8
info.lumi = 18.8
info.prettyName ='hscp search'
info.private = False
info.contact ='cms-pag-conveners-exotica@cern.ch'
info.arxiv = 'https://arxiv.org/abs/1305.0491'
info.comment ='Upper limits digitized from Track+TOF png plots. Used conservative gluino bounds (50% gluino ball probability).  '

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ txnames ++++++++++++++++++++
#+++++++ next txName block ++++++++++++++
HSCPM1 = dataset.addTxName('THSCPM1b')
HSCPM1.checked =''
HSCPM1.constraint = "[[],[]]"
HSCPM1.condition =None
HSCPM1.finalState = ['HSCP','HSCP']
HSCPM1.massConstraints = None
HSCPM1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-12-026/CMS-EXO-12-026_Figure_008-d.png'
HSCPM1.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
plane = HSCPM1.addMassPlane([[x],[x]])
plane.setSources(dataLabels=['upperLimits'],
                    dataFiles=['orig/CMS-EXO-12-026_Figure_008-d_stau.dat'],
                    dataFormats=['txt'],units=['pb'])


#+++++++ next txName block ++++++++++++++
RHadGM1 = dataset.addTxName('TRHadGM1')
RHadGM1.checked =''
RHadGM1.constraint = "[[],[]]"
RHadGM1.condition =None
RHadGM1.finalState = ['RHadronG','RHadronG']
RHadGM1.massConstraints = None
RHadGM1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-12-026/CMS-EXO-12-026_Figure_008-d.png'
RHadGM1.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
plane = RHadGM1.addMassPlane([[x],[x]])
plane.setSources(dataLabels=['upperLimits'],
                    dataFiles=['orig/CMS-EXO-12-026_Figure_008-d_gluino50.dat'],
                    dataFormats=['txt'],units=['pb'])

#+++++++ next txName block ++++++++++++++
RHadQM1 = dataset.addTxName('TRHadQM1')
RHadQM1.checked =''
RHadQM1.constraint = "[[],[]]"
RHadQM1.condition =None
RHadQM1.finalState = ['RHadronQ','RHadronQ']
RHadQM1.massConstraints = None
RHadQM1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-12-026/CMS-EXO-12-026_Figure_008-d.png'
RHadQM1.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
plane = RHadQM1.addMassPlane([[x],[x]])
plane.setSources(dataLabels=['upperLimits'],
                    dataFiles=['orig/CMS-EXO-12-026_Figure_008-d_stop.dat'],
                    dataFormats=['txt'],units=['pb'])


databaseCreator.create()

