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
info = MetaInfoInput('CMS-PAS-SUS-13-015')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13015'
info.sqrts = 8
info.lumi = 19.4
info.prettyName = '>= 5(1b-)jets + Etmiss'
info.private = False
info.arxiv = ''
info.contact =''
info.publication ='http://cds.cern.ch/record/1635353/files/SUS-13-015-pas.pdf'
info.comment ='PAS so no arxiv publication - only cds.cern. SuperseedeBy CMS-SUS-14-001 which has no public data.'
info.supersedes =''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked =''
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription =None
T2tt.condition =None
T2tt.source = "CMS"
T2tt.massConstraint = None
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked =''
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription =None
T2ttoff.condition =None
T2ttoff.massConstraint = [['dm <= 169.0'], ['dm <= 169.0']]
T2ttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure ='Fig. 8'
T2tt_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015/T2tt_alt_combined_all_fullToysXSEC.png'
T2tt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13015'
T2tt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015/T2tt_alt_combined_all_fullToysXSEC.C'
T2tt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015/T2tt_alt_combined_all_fullToysXSEC.C'
T2tt_1.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015/T2tt_alt_combined_all_fullToysXSEC.C'
T2tt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/CMS-SUS-13-015_from_cMacro.root', 'orig/CMS-SUS-13-015_from_cMacro.root', 'orig/CMS-SUS-13-015_from_cMacro.root', 'orig/CMS-SUS-13-015_from_cMacro.root', 'orig/CMS-SUS-13-015_from_cMacro.root', 'orig/CMS-SUS-13-015_from_cMacro.root', 'orig/CMS-SUS-13-015_from_cMacro.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['combined_expExclOneTimesProspino', 'combined_obsExclMinusSysErrProspino', 'combined_obsExclPlusSysErrProspino', 'combined_obsExclOneTimesProspino', 'combined_obsExclMinusSysErrProspino', 'combined_obsExclPlusSysErrProspino', 'combined_expLimit'])
T2ttoff.addMassPlane(T2tt_1)



databaseCreator.create()
