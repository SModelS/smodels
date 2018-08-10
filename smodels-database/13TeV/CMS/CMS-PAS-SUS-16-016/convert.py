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
info = MetaInfoInput('CMS-PAS-SUS-16-016')
info.url = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
info.sqrts = 13
info.lumi = 12.9
info.prettyName = '>= 1 jet + Etmiss, alpha_T'
info.private = False
info.comment = "Only CDS entry https://cds.cern.ch/record/2205163.  Superseded by CMS-SUS-16-033 and CMS-SUS-16-036."
info.supersededBy = 'CMS-SUS-16-033'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked =''
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription =None
T1bbbb.condition =None
T1bbbb.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane(2*[[x, y]])
T1bbbb_1.figure = "Figure 5-a"
T1bbbb_1.figureUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/CMS-PAS-SUS-16-016_Figure_005-a.png"
T1bbbb_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T1bbbb_1.histoDataUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T1bbbb_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/CMS-PAS-SUS-16-016_Figure_005-a.root"
T1bbbb_1.exclusionDataUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T1bbbb_1.setSources(dataLabels= ['expExclusion', 'expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T1bbbbExp.dat', 'orig/CMS-PAS-SUS-16-016_Figure_005-a.root', 'orig/T1bbbbObs.dat', 'orig/CMS-PAS-SUS-16-016_Figure_005-a.root'],
                 dataFormats= ['txt', 'root', 'txt', 'root'],objectNames= ['None', 'expected_limit_xs_noInterpolation', 'None', 'observed_limit_xs_noInterpolation'])

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked =''
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription =None
T1tttt.condition =None
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked =''
T1ttttoff.constraint = "[[['b','W','b','W']],[['b','W','b','W']]]" 
T1ttttoff.conditionDescription =None
T1ttttoff.condition =None
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure = "Figure 5-b"
T1tttt_1.figureUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html#Figure_005-b"
T1tttt_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T1tttt_1.histoDataUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T1tttt_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/CMS-PAS-SUS-16-016_Figure_005-b.root"
T1tttt_1.exclusionDataUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T1ttttExp.dat', 'orig/CMS-PAS-SUS-16-016_Figure_005-b.root', 'orig/T1ttttObs.dat', 'orig/CMS-PAS-SUS-16-016_Figure_005-b.root'],
                 dataFormats= ['txt', 'root', 'txt', 'root'],objectNames= ['None', 'expected_limit_xs_noInterpolation', 'None', 'observed_limit_xs_noInterpolation'])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2bb')
T2.checked =''
T2.constraint ="[[['b']],[['b']]]"
T2.conditionDescription =None
T2.condition =None
T2.source = "CMS"
T2.massConstraint = None
"""
T2off = dataset.addTxName('T2off')
T2off.checked =''
T2off.constraint =None
T2off.conditionDescription =None
T2off.condition =None
T2off.massConstraint = [['dm >= 0.0'], ['dm >= 0.0']]
T2off.source = "CMS"
"""
#+++++++ next mass plane block ++++++++++++++
T2_1 = T2.addMassPlane(2*[[x, y]])
T2_1.figure = "Figure 5-c"
T2_1.figureUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/CMS-PAS-SUS-16-016_Figure_005-c.png"
T2_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T2_1.histoDataUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T2_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/CMS-PAS-SUS-16-016_Figure_005-c.root"
T2_1.exclusionDataUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T2_1.setSources(dataLabels= ['expExclusion', 'expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2bbExp.dat', 'orig/CMS-PAS-SUS-16-016_Figure_005-c.root', 'orig/T2bbObs.dat', 'orig/CMS-PAS-SUS-16-016_Figure_005-c.root'],
                 dataFormats= ['txt', 'root', 'txt', 'root'],objectNames= ['None', 'expected_limit_xs_noInterpolation', 'None', 'observed_limit_xs_noInterpolation'])
# T2off.addMassPlane(T2_1)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked =''
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.conditionDescription =None
T2tt.condition =None
T2tt.source = "CMS"
T2tt.massConstraint = None
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked =''
T2ttoff.constraint = "[[['b','W']],[['b','W']]]"
T2ttoff.conditionDescription =None
T2ttoff.condition =None
T2ttoff.massConstraint = [['dm <= 169.0'], ['dm <= 169.0']]
T2ttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2tt_2 = T2tt.addMassPlane(2*[[x, y]])
T2tt_2.figure = "Figure 5-d"
T2tt_2.figureUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/CMS-PAS-SUS-16-016_Figure_005-d.png"
T2tt_2.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T2tt_2.histoDataUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T2tt_2.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/CMS-PAS-SUS-16-016_Figure_005-d.root"
T2tt_2.exclusionDataUrl ="http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-16-016/index.html"
T2tt_2.setSources(dataLabels= ['expExclusion', 'expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2ttExp.dat', 'orig/CMS-PAS-SUS-16-016_Figure_005-d.root', 'orig/T2ttObs.dat', 'orig/CMS-PAS-SUS-16-016_Figure_005-d.root'],
                 dataFormats= ['txt', 'root', 'txt', 'root'],objectNames= ['None', 'expected_limit_xs_noInterpolation', 'None', 'observed_limit_xs_noInterpolation'])
T2ttoff.addMassPlane(T2tt_2)



databaseCreator.create()
