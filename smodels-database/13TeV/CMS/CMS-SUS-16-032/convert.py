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
info = MetaInfoInput('CMS-SUS-16-032')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-032/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'Sbottom and compressed stop (jets + Etmiss)'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1707.07274'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Submitted to Phys. Lett. B.'
info.comment = 'Moriond 2017. The name for upper limits in the root files are: hXsec_exp_corr;1 and hXsec_exp_corr;2 with both containing same data (checked by Federico). Here used hXsec_exp_corr;2 for implementation.'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)


#+++++txName block +++++++++++++++++

T2bb=dataset.addTxName('T2bb')
T2bb.checked=''
T2bb.constraint="[[['b']],[['b']]]"
T2bb.condition=None
T2bb.conditionDescription = None
T2bb.source="CMS"
T2bb.massConstraint=None

#++++++next mass plane block+++++++++

T2bb_1 = T2bb.addMassPlane(2*[[x,y]])
T2bb_1.figure='Fig. 5'
T2bb_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-032/CMS-SUS-16-032_Figure_005.png'
T2bb_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-032/CMS-SUS-16-032_Figure_005.root'
T2bb_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-032/CMS-SUS-16-032_Figure_005.root'
T2bb_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-032/CMS-SUS-16-032_Figure_005.root'
T2bb_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-032_Figure_005.root','orig/CMS-SUS-16-032_Figure_005.root','orig/CMS-SUS-16-032_Figure_005.root','orig/CMS-SUS-16-032_Figure_005.root','orig/CMS-SUS-16-032_Figure_005.root','orig/CMS-SUS-16-032_Figure_005.root','orig/CMS-SUS-16-032_Figure_005.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['smoothed_Exp;1','smoothed_ExpM;1','smoothed_ExpP;1','smoothed_obs;1','smoothed_obsM;1','smoothed_obsP;1','hXsec_exp_corr;2'],units=[None,None,None,None,None,None,'pb'])


#+++++txName block +++++++++++++++++

T2cc=dataset.addTxName('T2cc')
T2cc.checked=''
T2cc.constraint="[[['c']],[['c']]]"
T2cc.condition=None
T2cc.conditionDescription = None
T2cc.source="CMS"

#++++++next mass plane block+++++++++

T2cc_1 = T2cc.addMassPlane(2*[[x,y]])
T2cc_1.figure='Fig. 6'
T2cc_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-032/CMS-SUS-16-032_Figure_006.png'
T2cc_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-032/CMS-SUS-16-032_Figure_006.root'
T2cc_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-032/CMS-SUS-16-032_Figure_006.root'
T2cc_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-032/CMS-SUS-16-032_Figure_006.root'
T2cc_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-032_Figure_006.root','orig/CMS-SUS-16-032_Figure_006.root','orig/CMS-SUS-16-032_Figure_006.root','orig/CMS-SUS-16-032_Figure_006.root','orig/CMS-SUS-16-032_Figure_006.root','orig/CMS-SUS-16-032_Figure_006.root','orig/CMS-SUS-16-032_Figure_006.root'],
                    dataFormats=['root','root','root','root','root','root','root'],objectNames=['smoothed_Exp;1','smoothed_ExpM;1','smoothed_ExpP;1','smoothed_obs;1','smoothed_obsM;1','smoothed_obsP;1','hXsec_exp_corr;2'],units=[None,None,None,None,None,None,'pb'])

databaseCreator.create()
