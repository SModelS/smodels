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
argparser.add_argument ('-t', '--ntoys',
            help = 'number of toys to throw',\
            type = int, default=200000  )
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

## databaseCreator.ncpus = 2
DataSetInput.ntoys = args.ntoys

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-13-012')
info.comment = 'T1,T2,T1tttt official efficiency maps from the CMS collaboration; T1btbt,T1bbbb,T2bb,T2tt,T5,T5bbbb,T5tttt,T6bbWW,TChiWZ,TChiZZ,TChiWW,T5WW and T5ZZ created by the SModelS collaboration using MadAnalysis5'
info.sqrts = '8.0'
info.private = False
info.lumi = '19.5'
info.prettyName = 'jet multiplicity + HTmiss'
info.publication = 'JHEP06(2014)055'
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13012'
info.arxiv = 'http://arxiv.org/abs/1402.4770'
info.contact = 'cms-pag-conveners-sus@NOSPAMcernSPAMNOT.ch, SModelS collaboration'
info.implementedBy = 'Federico A.'

def add ( dataset ):
    #+++++++ next txName block ++++++++++++++
    T5WW = dataset.addTxName('T5WW')
    T5WW.checked = ' '
    T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
    T5WW.conditionDescription ="None"
    T5WW.condition ="None"
    T5WW.massConstraint = None
    T5WW.source = 'SModelS'
    T5WW.dataUrl = None
    #+++++++ next txName block ++++++++++++++
    T5WWoff = dataset.addTxName('T5WWoff')
    T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
    T5WWoff.conditionDescription ="None"
    T5WWoff.condition =None
    T5WWoff.dataUrl = None
    T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
    T5WWoff.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++

    T5WW_x05 = T5WW.addMassPlane([[x,0.5*(x+y),y]]*2)
    #T5WW_x05.dataUrl = None
    T5WW_x05.addSource('obsExclusion', "orig/CMS_T5VV_x05.dat", "txt")
    T5WW_x05.addSource('efficiencyMap',"orig/T5WW_x05/MA5_EM_T5WW_x05_%s.dat" % dataset._name, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T5WW_x005 = T5WW.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
    #T5WW_x005.dataUrl = None
    T5WW_x005.addSource('efficiencyMap',"orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_%s.dat" % dataset._name, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T5WW_x095 = T5WW.addMassPlane([[x,0.95*x + 0.05*y,y]]*2)
    #T5WW_x095.dataUrl = None
    T5WW_x095.addSource('efficiencyMap',"orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_%s.dat" % dataset._name, "txt")
    T5WW_d075 = T5WW.addMassPlane([[x,y+75,y]]*2)
    #T5WW_d075.dataUrl = None
    T5WW_d075.addSource('efficiencyMap',"orig/cms_sus_13_012_T5WW_75_EM_MAPS/MA5_EM_T5WW_75_%s.dat" % dataset._name, "txt")
    T5WW_d010 = T5WW.addMassPlane([[x,y+10,y]]*2)
    #T5WW_d010.dataUrl = None
    T5WW_d010.addSource('efficiencyMap',"orig/cms_sus_13_012_T5WW_10_EM_MAPS/MA5_EM_T5WW_10_%s.dat" % dataset._name, "txt")
    T5WWoff.addMassPlane(T5WW_d075)
    T5WWoff.addMassPlane(T5WW_d010)
    T5WWoff.addMassPlane(T5WW_x05)
    T5WWoff.addMassPlane(T5WW_x005)
    T5WWoff.addMassPlane(T5WW_x095)
    TChiZZ = dataset.addTxName('TChiZZ')
    TChiZZ.checked = ''
    TChiZZ.dataUrl = None
    TChiZZ.constraint = "[[['Z']],[['Z']]]"
    TChiZZ.conditionDescription ="None"
    TChiZZ.condition ="None"
    TChiZZ.massConstraint = None
    TChiZZ.source = 'SModelS'
    TChiZZ.massConstraint = None
    TChiZZ_1 = TChiZZ.addMassPlane( [[x,y]]*2 )
    TChiZZ_1.addSource ( "efficiencyMap", "orig/cms_sus_13_012_TChiZZ_1_EM_MAPS/MA5_EM_TChiZZ_1_%s.dat" % dataset, "txt", objectName ="None", index = None )
    TChiZZ_1.dataUrl = None
    TChiWZ = dataset.addTxName('TChiWZ')
    TChiWZ.checked = ''
    TChiWZ.constraint = "[[['W']],[['Z']]]"
    TChiWZ.conditionDescription ="None"
    TChiWZ.condition ="None"
    TChiWZ.source = 'SModelS'
    TChiWZ.massConstraint = None
    TChiWZ_1 = TChiWZ.addMassPlane( [[x,y]]*2 )
    TChiWZ_1.addSource ( "efficiencyMap", "orig/cms_sus_13_012_TChiWZ_1_EM_MAPS/MA5_EM_TChiWZ_1_%s.dat" % dataset, "txt", objectName ="None", index = None )
    TChiWZ_1.dataUrl = None
    TChiWW = dataset.addTxName('TChiWW')
    TChiWW.checked = ''
    TChiWW.constraint = "[[['W']],[['W']]]"
    TChiWW.massConstraint = None
    TChiWW.conditionDescription ="None"
    TChiWW.condition ="None"
    TChiWW.source = 'SModelS'
    TChiWW_1 = TChiWW.addMassPlane( [[x,y]]*2 )
    TChiWW_1.addSource ( "efficiencyMap", "orig/cms_sus_13_012_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_%s.dat" % dataset, "txt", objectName ="None", index = None )
    TChiWW_1.figureUrl = "FIXME"
    TChiWW_1.dataUrl = None
    T5 = dataset.addTxName('T5')
    T5.checked = ''
    T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
    T5.conditionDescription ="None"
    T5.condition ="None"
    T5.massConstraint = None
    T5.source = 'SModelS'
    T5.dataUrl = None
    T5_x005 = T5.addMassPlane( [[x,0.05*x + 0.95*y,y]]*2 )
    T5_x005.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5_x005_EM_MAPS/MA5_EM_T5_x005_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T5_x05 = T5.addMassPlane( [[x,0.5*x + 0.5*y,y]]*2 )
    T5_x05.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5_x05_EM_MAPS/MA5_EM_T5_x05_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T5_x095 = T5.addMassPlane( [[x,0.95*x + 0.05*y,y]]*2 )
    T5_x095.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5_x095_EM_MAPS/MA5_EM_T5_x095_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T5bbbb = dataset.addTxName('T5bbbb')
    T5bbbb.checked = ''
    T5bbbb.constraint ="[[['b'],['b']],[['b'],['b']]]"
    T5bbbb.conditionDescription ="None"
    T5bbbb.condition ="None"
    T5bbbb.massConstraint = None
    T5bbbb.source = 'SModelS'
    T5bbbb.dataUrl = None
    T5bbbb_x005 = T5bbbb.addMassPlane( [[x,0.05*x + 0.95*y,y]]*2 )
    T5bbbb_x005.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5bbbb_x005_EM_MAPS/MA5_EM_T5bbbb_x005_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T5bbbb_x05 = T5bbbb.addMassPlane( [[x,0.5*x + 0.5*y,y]]*2 )
    T5bbbb_x05.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5bbbb_x05_EM_MAPS/MA5_EM_T5bbbb_x05_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T5bbbb_x095 = T5bbbb.addMassPlane( [[x,0.95*x + 0.05*y,y]]*2 )
    T5bbbb_x095.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5bbbb_x095_EM_MAPS/MA5_EM_T5bbbb_x095_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T5tttt = dataset.addTxName('T5tttt')
    T5tttt.checked = ''
    T5tttt.constraint ="[[['t'],['t']],[['t'],['t']]]"
    T5tttt.conditionDescription ="None"
    T5tttt.condition ="None"
    T5tttt.massConstraint = None
    T5tttt.source = 'SModelS'
    T5tttt.dataUrl = None
    T5tttt_x05 = T5tttt.addMassPlane( [[x,0.5*x + 0.5*y,y]]*2 )
    T5tttt_x05.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5tttt_x05_EM_MAPS/MA5_EM_T5tttt_x05_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T5tttt_p177 = T5tttt.addMassPlane( [[x, x-177.,y]]*2 )
    T5tttt_p177.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5tttt_DiffGluStop177_EM_MAPS/MA5_EM_T5tttt_DiffGluStop177_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T5tttt_m177 = T5tttt.addMassPlane( [[x, y+177.,y]]*2 )
    T5tttt_m177.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5tttt_DiffStopNeu177_EM_MAPS/MA5_EM_T5tttt_DiffStopNeu177_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T6bbWW = dataset.addTxName('T6bbWW')
    T6bbWW.checked = ''
    T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
    T6bbWW.conditionDescription ="None"
    T6bbWW.condition ="None"
    T6bbWW.massConstraint = None
    T6bbWW.source = 'SModelS'
    T6bbWW.dataUrl = None
    T6bbWW_x05 = T6bbWW.addMassPlane( [[x,0.5*x + 0.5*y,y]]*2 )
    T6bbWW_x05.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T6bbWW_x01 = T6bbWW.addMassPlane( [[x,0.1 *x + 0.9 *y,y]]*2  )
    T6bbWW_x01.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T6bbWW_x09 = T6bbWW.addMassPlane( [[x,0.9 *x + 0.1 *y,y]]*2  )
    T6bbWW_x09.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T6bbWW_d10 = T6bbWW.addMassPlane( [[x, y+10. , y ]]*2  )
    T6bbWW_d10.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_DiffChargNeu10_EM_MAPS/MA5_EM_T6bbWW_DiffChargNeu10_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T6bbWW_d75 = T6bbWW.addMassPlane( [[x, y+75. , y ]]*2  )
    T6bbWW_d75.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_75_EM_MAPS/MA5_EM_T6bbWW_75_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T6bbWWoff = dataset.addTxName('T6bbWWoff')
    T6bbWWoff.checked = ''
    T6bbWWoff.constraint = "2.23*[[['b'],['jet','jet']],[['b'],['jet','jet']]]"
    T6bbWWoff.conditionDescription ="None"
    T6bbWWoff.condition ="None"
    ## T6bbWWoff.massConstraint = None
    T6bbWWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
    T6bbWWoff.source = 'SModelS'
    T6bbWWoff.dataUrl = None
    T6bbWWoff_x05 = T6bbWWoff.addMassPlane( [[x,0.5*x + 0.5*y,y]]*2 )
    T6bbWWoff_x05.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T6bbWWoff_x01 = T6bbWWoff.addMassPlane( [[x,0.1 *x + 0.9 *y,y]]*2  )
    T6bbWWoff_x01.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T6bbWWoff_x09 = T6bbWWoff.addMassPlane( [[x,0.9 *x + 0.1 *y,y]]*2  )
    T6bbWWoff_x09.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T6bbWWoff_d10 = T6bbWWoff.addMassPlane( [[x, y+10. , y ]]*2  )
    T6bbWWoff_d10.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_DiffChargNeu10_EM_MAPS/MA5_EM_T6bbWW_DiffChargNeu10_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T6bbWWoff_d75 = T6bbWWoff.addMassPlane( [[x, y+75. , y ]]*2  )
    T6bbWWoff_d75.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_75_EM_MAPS/MA5_EM_T6bbWW_75_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T2 = dataset.addTxName('T2')
    T2.checked =''
    T2.constraint ="[[['jet']],[['jet']]]"
    T2.conditionDescription ="None"
    T2.condition ="None"
    T2.source = 'CMS'
    #+++++++ next mass plane block ++++++++++++++
    T2qq = T2.addMassPlane([[x,y]]*2)
    T2qq.figure = "Fig_7a"
    T2qq.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7a.pdf"
    T2qq.dataUrl = None
    T2qq.addSource('obsExclusion', "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_obsExclOneTimesProspino")
    T2qq.addSource('obsExclusionM1', "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_obsExclMinusSysErrProspino")
    T2qq.addSource('obsExclusionP1',"orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_obsExclPlusSysErrProspino")
    T2qq.addSource('expExclusion', "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_expExclOneTimesProspino")
    T2qq.addSource('expExclusionM1',"orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_expExclMinusOneSigmaProspino")
    T2qq.addSource('expExclusionP1',"orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_expExclPlusOneSigmaProspino")
    T2qq.addSource('efficiencyMap',"orig/SUS13012_XsecLimits_T2qq.root", "root", objectName ="h_EffAcc_%s" % dataset )
    #+++++++ next txName block ++++++++++++++
    T1tttt = dataset.addTxName('T1tttt')
    T1tttt.constraint ="[[['t','t']],[['t','t']]]"
    T1tttt.conditionDescription ="None"
    T1tttt.condition ="None"
    T1tttt.massConstraint = None
    T1tttt.source = 'CMS'
    #+++++++ next txName block ++++++++++++++
    T1ttttoff = dataset.addTxName('T1ttttoff')
    T1ttttoff.constraint = "[[['b','W','b','W']],[['b','W','b','W']]]"
    T1ttttoff.conditionDescription ="None"
    T1ttttoff.condition = "None"
    T1ttttoff.massConstraint = [['dm <= 338.']]*2
    T1ttttoff.source = 'CMS'
    #+++++++ next mass plane block ++++++++++++++
    T1bbbb = dataset.addTxName('T1bbbb')
    T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
    T1bbbb.conditionDescription ="None"
    T1bbbb.condition ="None"
    T1bbbb.massConstraint = None
    T1bbbb.source = 'SModelS'
    T1btbt = dataset.addTxName('T1btbt')
    T1btbt.constraint ="[[['b','t']],[['b','t']]]"
    T1btbt.conditionDescription ="None"
    T1btbt.condition ="None"
    T1btbt.massConstraint = None
    T1btbt.source = 'SModelS'
    T2bb = dataset.addTxName('T2bb')
    T2bb.constraint ="[[['b']],[['b']]]"
    T2bb.conditionDescription ="None"
    T2bb.condition ="None"
    T2bb.massConstraint = None
    T2bb.source = 'SModelS'
    T2bb.dataUrl = None
    T2tt = dataset.addTxName('T2tt')
    T2tt.constraint ="[[['t']],[['t']]]"
    T2tt.conditionDescription ="None"
    T2tt.condition ="None"
    T2tt.massConstraint = None
    T2tt.source = 'SModelS'
    T2ttoff = dataset.addTxName('T2ttoff')
    T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
    T2ttoff.conditionDescription ="None"
    T2ttoff.condition ="None"
    T2ttoff.massConstraint = [[ 'dm <= 169.' ]]*2
    T2ttoff.source = 'SModelS'
    T2ttoff.dataUrl = None
    T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
    T1tttt_1.figure = "Fig_7c"
    T1tttt_1.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7c.pdf"
    T1tttt_1.dataUrl = None
    T1tttt_1.addSource('obsExclusion', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_obsExclOneTimesProspino")
    T1tttt_1.addSource('obsExclusionM1', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_obsExclMinusSysErrProspino")
    T1tttt_1.addSource('obsExclusionP1', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_obsExclPlusSysErrProspino")
    T1tttt_1.addSource('expExclusion', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_expExclOneTimesProspino")
    T1tttt_1.addSource('expExclusionM1', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_expExclMinusOneSigmaProspino")
    T1tttt_1.addSource('expExclusionP1', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_expExclPlusOneSigmaProspino")
    T1tttt_1.addSource('efficiencyMap',"orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName ="h_EffAcc_%s"% dataset )
    T1ttttoff_1 = T1ttttoff.addMassPlane([[x,y]]*2)
    T1ttttoff_1.figure = "Fig_7c"
    T1ttttoff_1.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7c.pdf"
    T1ttttoff_1.dataUrl = None
    T1ttttoff_1.addSource('obsExclusion', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_obsExclOneTimesProspino")
    T1ttttoff_1.addSource('obsExclusionM1', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_obsExclMinusSysErrProspino")
    T1ttttoff_1.addSource('obsExclusionP1', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_obsExclPlusSysErrProspino")
    T1ttttoff_1.addSource('expExclusion', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_expExclOneTimesProspino")
    T1ttttoff_1.addSource('expExclusionM1', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_expExclMinusOneSigmaProspino")
    T1ttttoff_1.addSource('expExclusionP1', "orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName = "combined_expExclPlusOneSigmaProspino")
    T1ttttoff_1.addSource('efficiencyMap',"orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName ="h_EffAcc_%s" % dataset )
    T1bbbb_1 = T1bbbb.addMassPlane([[x,y]]*2)
    T1bbbb_1.figure = "Fig_7c"
    T1bbbb_1.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7c.pdf"
    T1bbbb_1.dataUrl = None
    T1bbbb_1.addSource('efficiencyMap',"orig/cms_sus_13_012_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T1btbt_1 = T1btbt.addMassPlane([[x,y]]*2)
    T1btbt_1.figure = "Fig_7c"
    T1btbt_1.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7c.pdf"
    T1btbt_1.dataUrl = None
    T1btbt_1.addSource('efficiencyMap',"orig/cms_sus_13_012_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T2bb_1 = T2bb.addMassPlane([[x,y]]*2)
    T2bb_1.addSource('efficiencyMap',"orig/cms_sus_13_012_T2bb_1_EM_MAPS/MA5_EM_T2bb_1_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
    T2tt_1.figure = None
    T2tt_1.figureUrl = None
    T2tt_1.dataUrl = None
    T2tt_1.addSource('efficiencyMap',"orig/cms_sus_13_012_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_%s.dat" % dataset, "txt", objectName ="None", index = None )
    T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
    T2ttoff_1.figure = "FIXME"
    T2ttoff_1.figureUrl = "FIXME"
    T2ttoff_1.dataUrl = None
    T2ttoff_1.addSource('efficiencyMap',"orig/cms_sus_13_012_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_%s.dat" % dataset, "txt", objectName ="None", index = None )
    #+++++++ next txName block ++++++++++++++
    T1 = dataset.addTxName('T1')
    T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
    T1.conditionDescription ="None"
    T1.condition ="None"
    T1.source = 'CMS'
    #+++++++ next mass plane block ++++++++++++++
    T1qqqq = T1.addMassPlane([[x,y]]*2)
    T1qqqq.figure = "Fig_7b"
    T1qqqq.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7b.pdf"
    T1qqqq.dataUrl = None
    T1qqqq.addSource('obsExclusion', "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_obsExclOneTimesProspino")
    T1qqqq.addSource('obsExclusionM1', "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_obsExclMinusSysErrProspino")
    T1qqqq.addSource('obsExclusionP1', "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_obsExclPlusSysErrProspino")
    T1qqqq.addSource('expExclusion', "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_expExclOneTimesProspino")
    T1qqqq.addSource('expExclusionM1', "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_expExclMinusOneSigmaProspino")
    T1qqqq.addSource('expExclusionP1', "orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName = "combined_expExclPlusOneSigmaProspino")
    T1qqqq.addSource('efficiencyMap',"orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName ="h_EffAcc_%s" % dataset )
    #+++++++ next txName block ++++++++++++++
    T5ZZ = dataset.addTxName('T5ZZ')
    T5ZZ.checked = ' '
    T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
    T5ZZ.conditionDescription ="None"
    T5ZZ.condition ="None"
    T5ZZ.massConstraint = None
    T5ZZ.dataUrl = None
    T5ZZ.source = 'SModelS'
    #+++++++ next txName block ++++++++++++++
    #+++++++ next mass plane block ++++++++++++++
    T5ZZ_x05 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
    T5ZZ_x05.addSource('obsExclusion', "orig/CMS_T5VV_x05.dat", "txt")
    T5ZZ_x05.addSource('efficiencyMap',"orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_%s.dat" % dataset, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T5ZZ_x005 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
    T5ZZ_x005.addSource('efficiencyMap',"orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_%s.dat" % dataset, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T5ZZ_x095 = T5ZZ.addMassPlane([[x,0.95*x + 0.05*y,y]]*2)
    T5ZZ_x095.addSource('efficiencyMap',"orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_%s.dat" % dataset, "txt")

datasets = { "3NJet6_500HT800_450MHT600": ( 454, 418, 66 ), 
             "3NJet6_1250HT1500_450MHTinf": ( 23, 17.6, 4.1 ),
             "3NJet6_1250HT1500_300MHT450": ( 38, 42.8, 9.5 ),
             "6NJet8_1000HT1250_200MHT300": ( 67, 70, 16 ),
             "6NJet8_800HT1000_300MHT450": ( 35, 28.6, 6.9 ),
            "6NJet8_500HT800_300MHT450": ( 62, 52, 12 ),
            "3NJet6_800HT1000_200MHT300": ( 808, 777, 107 ),
            "3NJet6_500HT800_200MHT300": ( 6159, 6088, 665 ),
            "6NJet8_1250HT1500_450MHTinf": ( 2, 0.5, 2.6 ),
            "3NJet6_1000HT1250_300MHT450": ( 129, 137, 20 ),
            "8NJetinf_1250HT1500_200MHTinf": ( 5, 7.1, 3.8 ),
            "8NJetinf_1000HT1250_200MHTinf": ( 8, 5.6, 2.3 ),
            "3NJet6_800HT1000_300MHT450": ( 305, 330, 40 ),
            "8NJetinf_800HT1000_200MHTinf": ( 9, 8.3, 3.4 ),
            "6NJet8_1500HTinf_300MHTinf": ( 3, 7.9, 3.6 ),
            "3NJet6_800HT1000_600MHTinf": ( 52, 54.8, 9.7 ),
            "6NJet8_500HT800_450MHTinf": ( 9, 0.8, 3.3 ),
            "6NJet8_1000HT1250_300MHT450": ( 20, 21.6, 5.8 ),
            "6NJet8_800HT1000_450MHTinf": ( 4, 6.0, 2.8 ),
            "3NJet6_1500HTinf_200MHT300": ( 94, 86, 17 ),
            "3NJet6_500HT800_600MHTinf": ( 62, 57.4, 11.2 ),
            "3NJet6_1250HT1500_200MHT300": ( 98, 109, 18 ),
            "3NJet6_1000HT1250_600MHTinf": ( 32, 22.8, 5.2 ),
            "3NJet6_1000HT1250_200MHT300": ( 335, 305, 41 ),
            "8NJetinf_1500HTinf_200MHTinf": ( 2, 3.3, 4.7 ),
            "6NJet8_1500HTinf_200MHT300": ( 18, 21.1, 8.1 ),
            "6NJet8_1250HT1500_200MHT300": ( 24, 28., 8.2 ),
            "3NJet6_800HT1000_450MHT600": ( 124, 108, 15 ),
            "6NJet8_800HT1000_200MHT300": ( 111, 124, 29 ),
            "6NJet8_500HT800_200MHT300": ( 266, 290, 65 ),
            "3NJet6_500HT800_300MHT450": ( 2305, 2278, 266 ),
            "8NJetinf_500HT800_200MHTinf": ( 8, 4.8, 2.3 ),
            "3NJet6_1000HT1250_450MHT600": ( 34, 32.3, 6.1 ),
            "6NJet8_1000HT1250_450MHTinf": ( 4, 2.2, 3.8 ),
            "6NJet8_1250HT1500_300MHT450": ( 5, 9.4, 3.6 ),
            "3NJet6_1500HTinf_300MHTinf": ( 39, 29.7, 5.8 )
}

for name, numbers in datasets.items():
    #+++++++ dataset block ++++++++++++++
    dataset = DataSetInput( name )
    dataset.setInfo(dataType = 'efficiencyMap', dataId = name, 
            observedN = numbers[0], expectedBG = numbers[1], bgError = numbers[2] )
    add ( dataset )

databaseCreator.create()
