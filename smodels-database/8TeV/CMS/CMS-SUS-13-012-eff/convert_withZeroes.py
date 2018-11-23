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
            type = int, default= 10  )
args = argparser.parse_args()

'''
if args.utilsPath:
    utilsPath = args.utilsPath
else:
    databaseRoot = '../../../'
    sys.path.append(os.path.abspath(databaseRoot))
    from utilsPath import utilsPath
    utilsPath = databaseRoot + utilsPath
if args.smodelsPath:
    sys.path.append(os.path.abspath(args.smodelsPath))
'''
    
utilsPath = '../../../../smodels-utils'
    
sys.path.append(os.path.abspath(utilsPath))
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z

## databaseCreator.ncpus = 2
DataSetInput.ntoys = args.ntoys

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-13-012')
info.comment = 'T1,T2,T1tttt official efficiency maps from the CMS collaboration; T5WW and T5ZZ created by the SModelS collaboration using MadAnalysis5'
info.sqrts = '8.0'
info.private = False
info.lumi = '19.5'
info.publication = 'JHEP06(2014)055'
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13012'
info.arxiv = 'http://arxiv.org/abs/1402.4770'
info.contact = 'cms-pag-conveners-sus@NOSPAMcernSPAMNOT.ch, SModelS collaboration'
info.implementedBy = 'Federico A.'

def add ( dataset ):
    dataset_n = dataset._name.replace('SR_','')

    '''
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
    T5WW_x05.addSource('efficiencyMap',"orig/T5WW_x05/MA5_EM_T5WW_x05_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T5WW_x005 = T5WW.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
    #T5WW_x005.dataUrl = None
    T5WW_x005.addSource('efficiencyMap',"orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T5WW_x095 = T5WW.addMassPlane([[x,0.95*x + 0.05*y,y]]*2)
    #T5WW_x095.dataUrl = None
    T5WW_x095.addSource('efficiencyMap',"orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_%s.dat" % dataset_n, "txt")
    T5WW_d077 = T5WW.addMassPlane([[x,y+77,y]]*2)
    #T5WW_d077.dataUrl = None
    T5WW_d077.addSource('efficiencyMap',"orig/cms_sus_13_012_T5WW_77_EM_MAPS/MA5_EM_T5WW_77_%s.dat" % dataset_n, "txt")
    T5WW_d010 = T5WW.addMassPlane([[x,y+10,y]]*2)
    #T5WW_d010.dataUrl = None
    T5WW_d010.addSource('efficiencyMap',"orig/cms_sus_13_012_T5WW_10_EM_MAPS/MA5_EM_T5WW_10_%s.dat" % dataset_n, "txt")
    T5WWoff.addMassPlane(T5WW_d077)
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
    TChiZZ_1.addSource ( "efficiencyMap", "orig/cms_sus_13_012_TChiZZ_1_EM_MAPS/MA5_EM_TChiZZ_1_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    TChiZZ_1.dataUrl = None
    TChiWZ = dataset.addTxName('TChiWZ')
    TChiWZ.checked = ''
    TChiWZ.constraint = "[[['W']],[['Z']]]"
    TChiWZ.conditionDescription ="None"
    TChiWZ.condition ="None"
    TChiWZ.source = 'SModelS'
    TChiWZ.massConstraint = None
    TChiWZ_1 = TChiWZ.addMassPlane( [[x,y]]*2 )
    TChiWZ_1.addSource ( "efficiencyMap", "orig/cms_sus_13_012_TChiWZ_1_EM_MAPS/MA5_EM_TChiWZ_1_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    TChiWZ_1.dataUrl = None
    TChiWW = dataset.addTxName('TChiWW')
    TChiWW.checked = ''
    TChiWW.constraint = "[[['W']],[['W']]]"
    TChiWW.massConstraint = None
    TChiWW.conditionDescription ="None"
    TChiWW.condition ="None"
    TChiWW.source = 'SModelS'
    TChiWW_1 = TChiWW.addMassPlane( [[x,y]]*2 )
    TChiWW_1.addSource ( "efficiencyMap", "orig/cms_sus_13_012_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    TChiWW_1.figureUrl = "FIXME"
    TChiWW_1.dataUrl = None
    T5bbbb = dataset.addTxName('T5bbbb')
    T5bbbb.checked = ''
    T5bbbb.constraint ="[[['b'],['b']],[['b'],['b']]]"
    T5bbbb.conditionDescription ="None"
    T5bbbb.condition ="None"
    T5bbbb.massConstraint = None
    T5bbbb.source = 'SModelS'
    T5bbbb.dataUrl = None
    T5bbbb_x005 = T5bbbb.addMassPlane( [[x,0.05*x + 0.95*y,y]]*2 )
    T5bbbb_x005.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5bbbb_x005_EM_MAPS/MA5_EM_T5bbbb_x005_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T5bbbb_x05 = T5bbbb.addMassPlane( [[x,0.5*x + 0.5*y,y]]*2 )
    T5bbbb_x05.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5bbbb_x05_EM_MAPS/MA5_EM_T5bbbb_x05_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T5bbbb_x095 = T5bbbb.addMassPlane( [[x,0.95*x + 0.05*y,y]]*2 )
    T5bbbb_x095.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5bbbb_x095_EM_MAPS/MA5_EM_T5bbbb_x095_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    
    
    T5tttt = dataset.addTxName('T5tttt')
    T5tttt.checked = ''
    T5tttt.constraint ="[[['t'],['t']],[['t'],['t']]]"
    T5tttt.conditionDescription ="None"
    T5tttt.condition ="None"
    T5tttt.massConstraint = None
    T5tttt.source = 'SModelS'
    T5tttt.dataUrl = None
    T5tttt_x05 = T5tttt.addMassPlane( [[x,0.5*x + 0.5*y,y]]*2 )
    T5tttt_x05.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5tttt_x05_EM_MAPS/MA5_EM_T5tttt_x05_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T5tttt_p177 = T5tttt.addMassPlane( [[x, x-177,y]]*2 )
    T5tttt_p177.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5tttt_DiffGluStop177_EM_MAPS/MA5_EM_T5tttt_DiffGluStop177_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T5tttt_m177 = T5tttt.addMassPlane( [[x, y+177,y]]*2 )
    T5tttt_m177.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5tttt_DiffStopNeu177_EM_MAPS/MA5_EM_T5tttt_DiffStopNeu177_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    
    T6bbWW = dataset.addTxName('T6bbWW')
    T6bbWW.checked = ''
    T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
    T6bbWW.conditionDescription ="None"
    T6bbWW.condition ="None"
    T6bbWW.massConstraint = None
    T6bbWW.source = 'SModelS'
    T6bbWW.dataUrl = None
    T6bbWW_x05 = T6bbWW.addMassPlane( [[x,0.5*x + 0.5*y,y]]*2 )
    T6bbWW_x05.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T6bbWW_x01 = T6bbWW.addMassPlane( [[x,0.1 *x + 0.9 *y,y]]*2  )
    T6bbWW_x01.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T6bbWW_x09 = T6bbWW.addMassPlane( [[x,0.9 *x + 0.1 *y,y]]*2  )
    T6bbWW_x09.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T6bbWW_d10 = T6bbWW.addMassPlane( [[x, y+10 , y ]]*2  )
    T6bbWW_d10.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_DiffChargNeu10_EM_MAPS/MA5_EM_T6bbWW_DiffChargNeu10_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T6bbWW_d77 = T6bbWW.addMassPlane( [[x, y+77 , y ]]*2  )
    T6bbWW_d77.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_DiffChargNeu77_EM_MAPS/MA5_EM_T6bbWW_DiffChargNeu77_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
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
    T6bbWWoff_x05.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T6bbWWoff_x01 = T6bbWWoff.addMassPlane( [[x,0.1 *x + 0.9 *y,y]]*2  )
    T6bbWWoff_x01.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T6bbWWoff_x09 = T6bbWWoff.addMassPlane( [[x,0.9 *x + 0.1 *y,y]]*2  )
    T6bbWWoff_x09.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T6bbWWoff_d10 = T6bbWWoff.addMassPlane( [[x, y+10 , y ]]*2  )
    T6bbWWoff_d10.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_DiffChargNeu10_EM_MAPS/MA5_EM_T6bbWW_DiffChargNeu10_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T6bbWWoff_d77 = T6bbWWoff.addMassPlane( [[x, y+77 , y ]]*2  )
    T6bbWWoff_d77.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T6bbWW_DiffChargNeu77_EM_MAPS/MA5_EM_T6bbWW_DiffChargNeu77_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    
       
    T2 = dataset.addTxName('T2')
    T2.checked =''
    T2.constraint ="[[['jet']],[['jet']]]"
    T2.conditionDescription ="None"
    T2.condition ="None"
    T2.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++
    T2qq = T2.addMassPlane([[x,y]]*2)
    T2qq.figure = "None"
    T2qq.figureUrl = "None"
    T2qq.dataUrl = None
    T2qq.addSource('obsExclusion', "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_obsExclOneTimesProspino")
    T2qq.addSource('obsExclusionM1', "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_obsExclMinusSysErrProspino")
    T2qq.addSource('obsExclusionP1',"orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_obsExclPlusSysErrProspino")
    T2qq.addSource('expExclusion', "orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_expExclOneTimesProspino")
    T2qq.addSource('expExclusionM1',"orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_expExclMinusOneSigmaProspino")
    T2qq.addSource('expExclusionP1',"orig/SUS13012_XsecLimits_T2qq.root", "root", objectName = "combined_expExclPlusOneSigmaProspino")

    T2qq.addSource('efficiencyMap',"orig/MA5_T2/MA5_EM_T2_Results_SR%s.dat"% dataset_n, "txt" , objectName ="None", index = None)
    
    #T2qq.addSource('efficiencyMap',"orig/SUS13012_XsecLimits_T2qq.root", "root", objectName ="h_EffAcc_%s" % dataset )  ### old official CMS data
    TGQ = dataset.addTxName('T3GQon')
    TGQ.checked =''
    TGQ.constraint ="[[['jet']],[['jet'],['jet']]]"
    TGQ.conditionDescription ="None"
    TGQ.condition ="None"
    TGQ.source = 'SModelS'
    TGQ.dataUrl = None
    for num in  [name.split('_')[1] for name in os.listdir('orig/T3GQon') if 'TGQ' in name]: # listing the gluino masses from the name in the directories
    #for num in  ['200','300']:
        a = TGQ.addMassPlane( [ [x,y], [eval(num),x,y] ] )
        a.addSource('efficiencyMap',"orig/T3GQon/TGQ_"+num+"_cms_sus_13_012/MA5_EM_TGQon_MAPS_SR%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    

    T5 = dataset.addTxName('T5')                                                                                                                                              
    T5.checked = ''                                                                                                                                                            
    T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"                                                                                                                     
    T5.conditionDescription ="None"                                                                                                                                           
    T5.condition ="None"                                                                                                                                                      
    T5.massConstraint = None                                                                                                                                                   
    T5.source = 'SModelS'                                                                                                                                                      
    T5.dataUrl = None                                                                                                                                                          
    T5_x005 = T5.addMassPlane( [[x,0.05*x + 0.95*y,y]]*2 )                                                                                                                     
    T5_x005.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5_x005_EM_MAPS/MA5_EM_T5_x005_%s.dat" % dataset_n, "txt", objectName ="None", index = None )                      
    T5_x05 = T5.addMassPlane( [[x,0.5*x + 0.5*y,y]]*2 )                                                                                                                            
    T5_x05.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5_x05_EM_MAPS/MA5_EM_T5_x05_%s.dat" % dataset_n, "txt", objectName ="None", index = None )                         
    T5_x095 = T5.addMassPlane( [[x,0.95*x + 0.05*y,y]]*2 )                                                                                                                       
    T5_x095.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5_x095_EM_MAPS/MA5_EM_T5_x095_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    
    T5_DeltaGluSq_5 = T5.addMassPlane( [[x, x-5 ,y]]*2 )
    T5_DeltaGluSq_5.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5_DeltaGluSq_5/MA5_EM_T5_MAPS_SR%s.dat" % dataset_n, "txt", objectName ="None", index = None )

    T5_DeltaSqNeu_5 = T5.addMassPlane( [[x, y+5 ,y]]*2 )
    T5_DeltaSqNeu_5.addSource ( "efficiencyMap", "orig/cms_sus_13_012_T5_DeltaSqNeu_5/MA5_EM_T5_MAPS_SR%s.dat" % dataset_n, "txt", objectName ="None", index = None )

    

    
    
    
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
    T1tttt_1.addSource('efficiencyMap',"orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName ="h_EffAcc_%s"% dataset_n )
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
    T1ttttoff_1.addSource('efficiencyMap',"orig/SUS13012_XsecLimits_T1tttt.root", "root", objectName ="h_EffAcc_%s" % dataset_n )
    T1bbbb_1 = T1bbbb.addMassPlane([[x,y]]*2)
    T1bbbb_1.figure = "Fig_7c"
    T1bbbb_1.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7c.pdf"
    T1bbbb_1.dataUrl = None
    T1bbbb_1.addSource('efficiencyMap',"orig/cms_sus_13_012_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T1btbt_1 = T1btbt.addMassPlane([[x,y]]*2)
    T1btbt_1.figure = "Fig_7c"
    T1btbt_1.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7c.pdf"
    T1btbt_1.dataUrl = None
    T1btbt_1.addSource('efficiencyMap',"orig/cms_sus_13_012_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T2bb_1 = T2bb.addMassPlane([[x,y]]*2)
    T2bb_1.addSource('efficiencyMap',"orig/cms_sus_13_012_T2bb_1_EM_MAPS/MA5_EM_T2bb_1_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
    T2tt_1.figure = None
    T2tt_1.figureUrl = None
    T2tt_1.dataUrl = None
    T2tt_1.addSource('efficiencyMap',"orig/cms_sus_13_012_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
    T2ttoff_1.figure = "FIXME"
    T2ttoff_1.figureUrl = "FIXME"
    T2ttoff_1.dataUrl = None
    T2ttoff_1.addSource('efficiencyMap',"orig/cms_sus_13_012_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
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
    T1qqqq.addSource('efficiencyMap',"orig/SUS13012_XsecLimits_T1qqqq.root", "root", objectName ="h_EffAcc_%s" % dataset_n )
    
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
    T5ZZ_x05.addSource('efficiencyMap',"orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T5ZZ_x005 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
    T5ZZ_x005.addSource('efficiencyMap',"orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T5ZZ_x095 = T5ZZ.addMassPlane([[x,0.95*x + 0.05*y,y]]*2)
    T5ZZ_x095.addSource('efficiencyMap',"orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_%s.dat" % dataset_n, "txt")
    '''

    ### NEW MODELS
    '''    
    TGN = dataset.addTxName('TGN')                                                                                                                                              
    TGN.checked = ' '                                                                                                                                                             
    TGN.constraint ="[[['lsp'],['jet','jet']],[['lsp'],['jet','jet']]]"                                                                                                           
    TGN.conditionDescription ="None"                                                                                                                                             
    TGN.condition ="None"                                                                                                                                                        
    TGN.massConstraint = None                                                                                                                                                    
    TGN.dataUrl = None                                                                                                                                                            
    TGN.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++                                                                                                                             
    TGN_1 = TGN.addMassPlane([[y,x,y]*2])        
    TGN_1.addSource('efficiencyMap',"orig/cms_sus_13_012_TGN_1_EM_MAPS/MA5_EM_TGN_1_%s.dat" % dataset_n, "txt")
    '''
    """
    T2bt = dataset.addTxName('T2bt')
    T2bt.checked = ' '
    T2bt.constraint ="[[['b']],[['t']]]"
    T2bt.conditionDescription ="None"
    T2bt.condition ="None"
    T2bt.massConstraint = None
    T2bt.dataUrl = None
    T2bt.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                   
    T2bt_1 = T2bt.addMassPlane([[x,y]]*2)
    T2bt_1.addSource('efficiencyMap',"orig/T2bt_EM/cms_sus_13_012_T2bt_1_EM_MAPS/MA5_EM_T2bt_1_%s.dat" % dataset_n, "txt") 

    T6WWbb = dataset.addTxName('T6WWbb')
    T6WWbb.checked = ''
    T6WWbb.constraint ="[[['W'],['b']],[['W'],['b']]]"
    T6WWbb.conditionDescription ="None"
    T6WWbb.condition ="None"
    T6WWbb.massConstraint = None
    T6WWbb.dataUrl = None
    T6WWbb.source = 'SModelS'
    
    T6WWoffbb = dataset.addTxName('T6WWoffbb')                                                                                                                                 
    T6WWoffbb.checked = ''                                                                                                                                                          
    T6WWoffbb.constraint = "2.23*[[['jet','jet'],['b']],[['jet','jet'],['b']]]"                                                                                                     
    T6WWoffbb.conditionDescription ="None"
    T6WWoffbb.condition ="None"
    T6WWoffbb.massConstraint = [['dm <= 76.','dm >= 5.0',]]*2
    T6WWoffbb.dataUrl = None
    T6WWoffbb.source = 'SModelS'
    
    #+++++++ next mass plane block ++++++++++++++                                                                                                                               
    T6WWbb_x005 = T6WWbb.addMassPlane([[x, 0.05*x + (1-0.05)*y ,  y]]*2)
    T6WWbb_x005.addSource('efficiencyMap',"orig/T6WWbb_EM/cms_sus_13_012_T6WWbb_x005_EM_MAPS/MA5_EM_T6WWbb_x005_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                               
    T6WWbb_x075 = T6WWbb.addMassPlane([[x, 0.75*x + (1-0.75)*y ,  y]]*2)
    T6WWbb_x075.addSource('efficiencyMap',"orig/T6WWbb_EM/cms_sus_13_012_T6WWbb_x075_EM_MAPS/MA5_EM_T6WWbb_x075_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6WWbb_X050 = T6WWbb.addMassPlane([[x, 0.50*x + (1-0.50)*y ,  y]]*2)
    T6WWbb_X050.addSource('efficiencyMap',"orig/T6WWbb_EM/cms_sus_13_012_T6WWbb_X050_EM_MAPS/MA5_EM_T6WWbb_X050_%s.dat" % dataset_n, "txt") 
    #+++++++ next mass plane block ++++++++++++++                                                                                                                               
    T6WWbb_DeltaStopSbot10 = T6WWoffbb.addMassPlane([[x, x-10 , y]]*2)
    T6WWbb_DeltaStopSbot10.addSource('efficiencyMap',"orig/T6WWbb_EM/cms_sus_13_012_T6WWbb_DeltaStopSbot10_EM_MAPS/MA5_EM_T6WWbb_DeltaStopSbot10_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                
    T6WWbb_DeltaStopSbot40 = T6WWoffbb.addMassPlane([[x, x-40 , y]]*2)
    T6WWbb_DeltaStopSbot40.addSource('efficiencyMap',"orig/T6WWbb_EM/cms_sus_13_012_T6WWbb_DeltaStopSbot40_EM_MAPS/MA5_EM_T6WWbb_DeltaStopSbot40_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                  
    T6WWbb_DeltaStopSbot75 = T6WWoffbb.addMassPlane([[x, x-75 , y]]*2)
    T6WWbb_DeltaStopSbot75.addSource('efficiencyMap',"orig/T6WWbb_EM/cms_sus_13_012_T6WWbb_DeltaStopSbot75_EM_MAPS/MA5_EM_T6WWbb_DeltaStopSbot75_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                            
    T6WWbb_DeltaSbotNeu5 = T6WWbb.addMassPlane([[x, y+5 , y]]*2)
    T6WWbb_DeltaSbotNeu5.addSource('efficiencyMap',"orig/T6WWbb_EM/cms_sus_13_012_T6WWbb_DeltaSbotNeu5_EM_MAPS/MA5_EM_T6WWbb_DeltaSbotNeu5_%s.dat" % dataset_n, "txt")

# *************************************************************

    T6ZZbb = dataset.addTxName('T6ZZbb')
    T6ZZbb.checked = ''
    T6ZZbb.constraint ="[[['Z'],['b']],[['Z'],['b']]]"
    T6ZZbb.conditionDescription ="None"
    T6ZZbb.condition ="None"
    T6ZZbb.massConstraint = None
    T6ZZbb.dataUrl = None
    T6ZZbb.source = 'SModelS'

    T6ZZoffbb = dataset.addTxName('T6ZZoffbb')
    T6ZZoffbb.checked = ''                                                                                                                                                   
    T6ZZoffbb.constraint = "2.1*[[['jet','jet'],['b']],[['jet','jet'],['b']]]"    # using BR(z->hadrons)= 69.2                                                               
    T6ZZoffbb.conditionDescription ="None"
    T6ZZoffbb.condition ="None"
    T6ZZoffbb.massConstraint = [['dm <= 85.','dm >= 5.0',]]*2
    T6ZZoffbb.dataUrl = None
    T6ZZoffbb.source = 'SModelS'

    #+++++++ next mass plane block ++++++++++++++                                                                                                                             
    T6ZZbb_x010 = T6WWbb.addMassPlane([[x, 0.1*x + (1-0.10)*y ,  y]]*2)
    T6ZZbb_x010.addSource('efficiencyMap',"orig/T6ZZbb_EM/cms_sus_13_012_T6ZZbb_x010_EM_MAPS/MA5_EM_T6ZZbb_x010_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                               
    T6ZZbb_X050 = T6ZZbb.addMassPlane([[x, 0.50*x + (1-0.50)*y ,  y]]*2)
    T6ZZbb_X050.addSource('efficiencyMap',"orig/T6ZZbb_EM/cms_sus_13_012_T6ZZbb_x050_EM_MAPS/MA5_EM_T6ZZbb_x050_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                             
    T6ZZoffbb_DeltaSbot2Sbot40 = T6ZZoffbb.addMassPlane([[x, x-40 , y]]*2)
    T6ZZoffbb_DeltaSbot2Sbot40.addSource('efficiencyMap',"orig/T6ZZbb_EM/cms_sus_13_012_T6ZZbb_DeltaSbot2Sbot40_EM_MAPS/MA5_EM_T6ZZbb_DeltaSbot2Sbot40_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6ZZoffbb_DeltaSbot2Sbot5 = T6ZZoffbb.addMassPlane([[x, x-5 , y]]*2)
    T6ZZoffbb_DeltaSbot2Sbot5.addSource('efficiencyMap',"orig/T6ZZbb_EM/cms_sus_13_012_T6ZZbb_DeltaSbot2Sbot5_EM_MAPS/MA5_EM_T6ZZbb_DeltaSbot2Sbot5_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6ZZoffbb_DeltaSbot2Sbot85 = T6ZZoffbb.addMassPlane([[x, x-85 , y]]*2)
    T6ZZoffbb_DeltaSbot2Sbot85.addSource('efficiencyMap',"orig/T6ZZbb_EM/cms_sus_13_012_T6ZZbb_DeltaSbot2Sbot85_EM_MAPS/MA5_EM_T6ZZbb_DeltaSbot2Sbot85_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6ZZbb_DeltaSbot2Sbot95 = T6ZZbb.addMassPlane([[x, x-95 , y]]*2)
    T6ZZbb_DeltaSbot2Sbot95.addSource('efficiencyMap',"orig/T6ZZbb_EM/cms_sus_13_012_T6ZZbb_DeltaSbot2Sbot95_EM_MAPS/MA5_EM_T6ZZbb_DeltaSbot2Sbot95_%s.dat" % dataset_n, "txt")

# *************************************************************

    T6bbZZ = dataset.addTxName('T6bbZZ')
    T6bbZZ.checked = ''
    T6bbZZ.constraint ="[[['b'],['Z']],[['b'],['Z']]]"
    T6bbZZ.conditionDescription ="None"
    T6bbZZ.condition ="None"
    T6bbZZ.massConstraint = None
    T6bbZZ.dataUrl = None
    T6bbZZ.source = 'SModelS'
    
    T6bbZZoff = dataset.addTxName('T6bbZZoff')
    T6bbZZoff.checked = ''
    T6bbZZoff.constraint = "2.1*[[['b'],['jet','jet']],[['b'],['jet','jet']]]"    # using BR(z->hadrons)= 69.2                                                                    
    T6bbZZoff.conditionDescription ="None"
    T6bbZZoff.condition ="None"
    T6bbZZoff.massConstraint = [[ 'dm >= 5.0','dm <= 85.']]*2
    T6bbZZoff.dataUrl = None
    T6bbZZoff.source = 'SmodelS'
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                
    T6bbZZ_x025= T6bbZZ.addMassPlane([[x, 0.25*x + (1-0.25)*y ,  y]]*2)
    T6bbZZ_x025.addSource('efficiencyMap',"orig/T6bbZZ_EM/cms_sus_13_012_T6bbZZ_x025_EM_MAPS/MA5_EM_T6bbZZ_x025_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                              
    T6bbZZ_X095 = T6bbZZ.addMassPlane([[x, 0.95*x + (1-0.95)*y ,  y]]*2)
    T6bbZZ_X095.addSource('efficiencyMap',"orig/T6bbZZ_EM/cms_sus_13_012_T6bbZZ_x095_EM_MAPS/MA5_EM_T6bbZZ_x095_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6bbZZ_X050 = T6bbZZ.addMassPlane([[x, 0.50*x + (1-0.50)*y ,  y]]*2)
    T6bbZZ_X050.addSource('efficiencyMap',"orig/T6bbZZ_EM/cms_sus_13_012_T6bbZZ_x050_EM_MAPS/MA5_EM_T6bbZZ_x050_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++             
    T6bbZZoff_DeltaN2LSP50 = T6bbZZoff.addMassPlane([[x, y + 50.0 ,  y]]*2) 
    T6bbZZoff_DeltaN2LSP50.addSource('efficiencyMap',"orig/T6bbZZ_EM/cms_sus_13_012_T6bbZZ_DeltaN2LSP50_EM_MAPS/MA5_EM_T6bbZZ_DeltaN2LSP50_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6bbZZoff_DeltaN2LSP5 = T6bbZZoff.addMassPlane([[x, y + 5.0 ,  y]]*2) 
    T6bbZZoff_DeltaN2LSP5.addSource('efficiencyMap',"orig/T6bbZZ_EM/cms_sus_13_012_T6bbZZ_DeltaN2LSP5_EM_MAPS/MA5_EM_T6bbZZ_DeltaN2LSP5_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6bbZZoff_DeltaN2LSP75 = T6bbZZoff.addMassPlane([[x, y + 75.0 ,  y]]*2) 
    T6bbZZoff_DeltaN2LSP75.addSource('efficiencyMap',"orig/T6bbZZ_EM/cms_sus_13_012_T6bbZZ_DeltaN2LSP75_EM_MAPS/MA5_EM_T6bbZZ_DeltaN2LSP75_%s.dat" % dataset_n, "txt")

# *************************************************************

    T6WWtt = dataset.addTxName('T6WWtt')
    T6WWtt.checked = ''
    T6WWtt.constraint ="[[['W'],['t']],[['W'],['t']]]"
    T6WWtt.conditionDescription ="None"
    T6WWtt.condition ="None"
    T6WWtt.massConstraint = [['dm >= 80.0','dm >= 175.']]*2
    T6WWtt.dataUrl = None
    T6WWtt.source = 'SModelS'

    T6WWttoff = dataset.addTxName('T6WWttoff')
    T6WWttoff.checked = ''
    T6WWttoff.constraint = "[[['W'],['b','W']],[['W'],['b','W']]]"                                                                      
    T6WWttoff.conditionDescription ="None"
    T6WWttoff.condition ="None"
    T6WWttoff.massConstraint = [['dm >= 80.0','dm >= 85.']]*2
    T6WWttoff.dataUrl = None
    T6WWttoff.source = 'SmodelS'

    T6WWofftt = dataset.addTxName('T6WWofftt')
    T6WWofftt.checked = ''
    T6WWofftt.constraint = "2.23*[[['jet','jet'],['t']],[['jet','jet'],['t']]]"                                                                      
    T6WWofftt.conditionDescription ="None"
    T6WWofftt.condition ="None"
    T6WWofftt.massConstraint = [['dm <= 80.0','dm >= 175.']]*2
    T6WWofftt.dataUrl = None
    T6WWofftt.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                 
    T6WWtt_x025= T6WWtt.addMassPlane([[x, 0.25*x + (1-0.25)*y , y]]*2)
    T6WWtt_x025.addSource('efficiencyMap',"orig/T6WWtt_EM/cms_sus_13_012_T6WWtt_x025_EM_MAPS/MA5_EM_T6WWtt_x025_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                           
    T6WWtt_x050= T6WWtt.addMassPlane([[x, 0.50*x + (1-0.50)*y , y]]*2)
    T6WWtt_x050.addSource('efficiencyMap',"orig/T6WWtt_EM/cms_sus_13_012_T6WWtt_x050_EM_MAPS/MA5_EM_T6WWtt_x050_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                              
    T6WWtt_x075= T6WWtt.addMassPlane([[x, 0.75*x + (1-0.75)*y , y]]*2)
    T6WWtt_x075.addSource('efficiencyMap',"orig/T6WWtt_EM/cms_sus_13_012_T6WWtt_x075_EM_MAPS/MA5_EM_T6WWtt_x075_%s.dat" % dataset_n, "txt")
    
    #+++++++ next mass plane block ++++++++++++++
    T6WWttoff_DeltaStopNeu90 = T6WWttoff.addMassPlane([[x, y+90.0 , y]]*2)
    T6WWttoff_DeltaStopNeu90.addSource('efficiencyMap',"orig/T6WWtt_EM/cms_sus_13_012_T6WWtt_DeltaStopNeutralino90_EM_MAPS/MA5_EM_T6WWtt_DeltaStopNeutralino90_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6WWttoff_DeltaStopNeu130 = T6WWttoff.addMassPlane([[x, y+130.0 , y]]*2)
    T6WWttoff_DeltaStopNeu130.addSource('efficiencyMap',"orig/T6WWtt_EM/cms_sus_13_012_T6WWtt_DeltaStopNeutralino130_EM_MAPS/MA5_EM_T6WWtt_DeltaStopNeutralino130_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6WWttoff_DeltaStopNeu165 = T6WWttoff.addMassPlane([[x, y+165.0 , y]]*2)
    T6WWttoff_DeltaStopNeu165.addSource('efficiencyMap',"orig/T6WWtt_EM/cms_sus_13_012_T6WWtt_DeltaStopNeutralino165_EM_MAPS/MA5_EM_T6WWtt_DeltaStopNeutralino165_%s.dat" % dataset_n, "txt")

    #+++++++ next mass plane block ++++++++++++++                                                                             
    T6WWofftt_DeltaSbottomStop5 = T6WWofftt.addMassPlane([[x, x-5.0 , y]]*2)
    T6WWofftt_DeltaSbottomStop5.addSource('efficiencyMap',"orig/T6WWtt_EM/cms_sus_13_012_T6WWtt_DeltaSbottomStop5_EM_MAPS/MA5_EM_T6WWtt_DeltaSbottomStop5_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6WWofftt_DeltaSbottomStop40 = T6WWofftt.addMassPlane([[x, x-40.0 , y]]*2)
    T6WWofftt_DeltaSbottomStop40.addSource('efficiencyMap',"orig/T6WWtt_EM/cms_sus_13_012_T6WWtt_DeltaSbottomStop40_EM_MAPS/MA5_EM_T6WWtt_DeltaSbottomStop40_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                
    T6WWofftt_DeltaSbottomStop75 = T6WWofftt.addMassPlane([[x, x-75.0 , y]]*2)
    T6WWofftt_DeltaSbottomStop75.addSource('efficiencyMap',"orig/T6WWtt_EM/cms_sus_13_012_T6WWtt_DeltaSbottomStop75_EM_MAPS/MA5_EM_T6WWtt_DeltaSbottomStop75_%s.dat" % dataset_n, "txt")


# ********************************************************        
    
    T6ttWW = dataset.addTxName('T6ttWW')
    T6ttWW.checked = ''
    T6ttWW.constraint ="[[['t'],['W']],[['t'],['W']]]"
    T6ttWW.conditionDescription ="None"
    T6ttWW.condition ="None"
    T6ttWW.massConstraint = [['dm >= 175.0','dm >= 80.']]*2
    T6ttWW.dataUrl = None
    T6ttWW.source = 'SModelS'

    T6ttWWoff = dataset.addTxName('T6ttWWoff')
    T6ttWWoff.checked = ''
    T6ttWWoff.constraint = "[[['t'],['jet','jet']],['t'],['jet','jet']]]"                                                                      
    T6ttWWoff.conditionDescription ="None"
    T6ttWWoff.condition ="None"
    T6ttWWoff.massConstraint = [[ 'dm >= 175.0','dm <= 80.']]*2
    T6ttWWoff.dataUrl = None
    T6ttWWoff.source = 'SModelS'

    T6ttoffWW = dataset.addTxName('T6ttoffWW')
    T6ttoffWW.checked = ''
    T6ttoffWW.constraint = "[[['b','W'],['W']],[['b','W'],['W']]]"                                                                      
    T6ttoffWW.conditionDescription ="None"
    T6ttoffWW.condition ="None"
    T6ttoffWW.massConstraint = [[ 'dm >= 85.0','dm >= 80.']]*2
    T6ttoffWW.dataUrl = None
    T6ttoffWW.source = 'SModelS'

    #+++++++ next mass plane block ++++++++++++++                                                                                                                              
    T6ttWW_x025 = T6ttWW.addMassPlane([[x, 0.25*x + (1-0.25)*y , y]]*2)
    T6ttWW_x025.addSource('efficiencyMap',"orig/T6ttWW_EM/cms_sus_13_012_T6ttWW_x025_EM_MAPS/MA5_EM_T6ttWW_x025_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                            
    T6ttWW_x050 = T6ttWW.addMassPlane([[x, 0.50*x + (1-0.50)*y , y]]*2)
    T6ttWW_x050.addSource('efficiencyMap',"orig/T6ttWW_EM/cms_sus_13_012_T6ttWW_x050_EM_MAPS/MA5_EM_T6ttWW_x050_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                              
    T6ttWW_x075 = T6ttWW.addMassPlane([[x, 0.75*x + (1-0.75)*y , y]]*2)
    T6ttWW_x075.addSource('efficiencyMap',"orig/T6ttWW_EM/cms_sus_13_012_T6ttWW_x075_EM_MAPS/MA5_EM_T6ttWW_x075_%s.dat" % dataset_n, "txt")
    
    #+++++++ next mass plane block ++++++++++++++
    T6ttoffWW_DeltaSbotCharg130 = T6ttoffWW.addMassPlane([[x , x - 130.0, y]]*2)
    T6ttoffWW_DeltaSbotCharg130.addSource('efficiencyMap',"orig/T6ttWW_EM/cms_sus_13_012_T6ttWW_DeltaSbotCharg130_EM_MAPS/MA5_EM_T6ttWW_DeltaSbotCharg130_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                             
    T6ttoffWW_DeltaSbotCharg165 = T6ttoffWW.addMassPlane([[x , x - 165.0, y]]*2)
    T6ttoffWW_DeltaSbotCharg165.addSource('efficiencyMap',"orig/T6ttWW_EM/cms_sus_13_012_T6ttWW_DeltaSbotCharg165_EM_MAPS/MA5_EM_T6ttWW_DeltaSbotCharg165_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                             
    T6ttWWoff_DeltaChargNeu5 = T6ttWW.addMassPlane([[x , y + 5.0 , y]]*2)
    T6ttWWoff_DeltaChargNeu5.addSource('efficiencyMap',"orig/T6ttWW_EM/cms_sus_13_012_T6ttWW_DeltaChargNeu5_EM_MAPS/MA5_EM_T6ttWW_DeltaChargNeu5_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                           
    T6ttWWoff_DeltaChargNeu40 = T6ttWW.addMassPlane([[x , y + 40.0 , y]]*2)
    T6ttWWoff_DeltaChargNeu40.addSource('efficiencyMap',"orig/T6ttWW_EM/cms_sus_13_012_T6ttWW_DeltaChargNeu40_EM_MAPS/MA5_EM_T6ttWW_DeltaChargNeu40_%s.dat" % dataset_n, "txt")
    

# ********************************************************                                                                                                                          
    T6ZZtt = dataset.addTxName('T6ZZtt')
    T6ZZtt.checked = ''
    T6ZZtt.constraint ="[[['Z'],['t']],[['Z'],['t']]]"
    T6ZZtt.conditionDescription ="None"
    T6ZZtt.condition ="None"
    T6ZZtt.massConstraint = [['dm >= 90.0','dm >= 175.']]*2
    T6ZZtt.dataUrl = None
    T6ZZtt.source = 'SModelS'
    
    T6ZZofftt = dataset.addTxName('T6ZZofftt')
    T6ZZofftt.checked = ''
    T6ZZofftt.constraint ="2.1*[[['jet','jet'],['t']],[['jet','jet'],['t']]]"
    T6ZZofftt.conditionDescription ="None"
    T6ZZofftt.condition ="None"
    T6ZZofftt.massConstraint = [['dm <= 90.0','dm >= 175.']]*2
    T6ZZofftt.dataUrl = None
    T6ZZofftt.source = 'SModelS'

    T6ZZttoff = dataset.addTxName('T6ZZttoff')
    T6ZZttoff.checked = ''
    T6ZZttoff.constraint ="[[['Z'],['b','W']],[['Z'],['b','W']]]"
    T6ZZttoff.conditionDescription ="None"
    T6ZZttoff.condition ="None"
    T6ZZttoff.massConstraint = [['dm <= 90.0','dm <= 175.']]*2
    T6ZZttoff.dataUrl = None
    T6ZZttoff.source = 'SModelS'

    #+++++++ next mass plane block ++++++++++++++                                                                                                                         
    T6ZZtt_x010 = T6ZZtt.addMassPlane([[x, 0.10*x + (1-0.10)*y , y]]*2)
    T6ZZtt_x010.addSource('efficiencyMap',"orig/T6ZZtt_EM/cms_sus_13_012_T6ZZtt_x0.1_EM_MAPS/MA5_EM_T6ZZtt_x0.1_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                  
    T6ZZofftt_x095 = T6ZZofftt.addMassPlane([[x, 0.95*x + (1-0.95)*y , y]]*2)
    T6ZZofftt_x095.addSource('efficiencyMap',"orig/T6ZZtt_EM/cms_sus_13_012_T6ZZtt_x095_EM_MAPS/MA5_EM_T6ZZtt_x0.95_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                  
    T6ZZofftt_x050 = T6ZZofftt.addMassPlane([[x, 0.50*x + (1-0.50)*y , y]]*2)
    T6ZZofftt_x050.addSource('efficiencyMap',"orig/T6ZZtt_EM/cms_sus_13_012_T6ZZtt_x050_EM_MAPS/MA5_EM_T6ZZtt_x050_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                  
    T6ZZtt_x05 = T6ZZtt.addMassPlane([[x, 0.50*x + (1-0.50)*y , y]]*2)
    T6ZZtt_x05.addSource('efficiencyMap',"orig/T6ZZtt_EM/cms_sus_13_012_T6ZZtt_x0.5_EM_MAPS/MA5_EM_T6ZZtt_x0.5_%s.dat" % dataset_n, "txt")

    #+++++++ next mass plane block ++++++++++++++                                                                                                                            
    T6ZZofftt_DeltaStop2Stop5 = T6ZZofftt.addMassPlane([[x, x-5.0 , y]]*2)
    T6ZZofftt_DeltaStop2Stop5.addSource('efficiencyMap',"orig/T6ZZtt_EM/cms_sus_13_012_T6ZZtt_DeltaStop2Stop5_EM_MAPS/MA5_EM_T6ZZtt_DeltaStop2Stop5_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                 
    T6ZZofftt_DeltaStop2Stop45 = T6ZZofftt.addMassPlane([[x, x-45.0 , y]]*2)
    T6ZZofftt_DeltaStop2Stop45.addSource('efficiencyMap',"orig/T6ZZtt_EM/cms_sus_13_012_T6ZZtt_DeltaStop2Stop45_EM_MAPS/MA5_EM_T6ZZtt_DeltaStop2Stop45_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                 
    T6ZZofftt_DeltaStop2Stop80 = T6ZZofftt.addMassPlane([[x, x-80.0 , y]]*2)
    T6ZZofftt_DeltaStop2Stop80.addSource('efficiencyMap',"orig/T6ZZtt_EM/cms_sus_13_012_T6ZZtt_DeltaStop2Stop80_EM_MAPS/MA5_EM_T6ZZtt_DeltaStop2Stop80_%s.dat" % dataset_n, "txt")

    #+++++++ next mass plane block ++++++++++++++                                                                                                                                
    T6ZZttoff_DeltaStopNeu90 = T6ZZttoff.addMassPlane([[x, y + 90.0 , y]]*2)
    T6ZZttoff_DeltaStopNeu90.addSource('efficiencyMap',"orig/T6ZZtt_EM/cms_sus_13_012_T6ZZtt_DeltaStopNeu90_EM_MAPS/MA5_EM_T6ZZtt_DeltaStopNeu90_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                             
    T6ZZttoff_DeltaStopNeu135 = T6ZZttoff.addMassPlane([[x, y + 135.0 , y]]*2)
    T6ZZttoff_DeltaStopNeu135.addSource('efficiencyMap',"orig/T6ZZtt_EM/cms_sus_13_012_T6ZZtt_DeltaStopNeu135_EM_MAPS/MA5_EM_T6ZZtt_DeltaStopNeu135_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                                  
    T6ZZttoff_DeltaStopNeu160 = T6ZZttoff.addMassPlane([[x, y + 160.0 , y]]*2)
    T6ZZttoff_DeltaStopNeu160.addSource('efficiencyMap',"orig/T6ZZtt_EM/cms_sus_13_012_T6ZZtt_DeltaStopNeu160_EM_MAPS/MA5_EM_T6ZZtt_DeltaStopNeu160_%s.dat" % dataset_n, "txt")
    """
# **********************************************************

    T6ttZZ = dataset.addTxName('T6ttZZ')
    T6ttZZ.checked = ''
    T6ttZZ.constraint ="[[['t'],['Z']],[['t'],['Z']]]"
    T6ttZZ.conditionDescription ="None"
    T6ttZZ.condition ="None"
    T6ttZZ.massConstraint = [['dm >= 175.0','dm >= 90.']]*2
    T6ttZZ.dataUrl = None
    T6ttZZ.source = 'SModelS'

    T6ttZZoff = dataset.addTxName('T6ttZZoff')
    T6ttZZoff.checked = ''
    T6ttZZoff.constraint ="2.1*[[['t'],['jet','jet']],[['t'],['jet','jet']]]"
    T6ttZZoff.conditionDescription ="None"
    T6ttZZoff.condition ="None"
    T6ttZZoff.massConstraint = [['dm >= 175.0','dm <=90.']]*2
    T6ttZZoff.dataUrl = None
    T6ttZZoff.source = 'SModelS'

    #T6ttoffZZ = dataset.addTxName('T6ttoffZZ')
    #T6ttoffZZ.checked = ''
    #T6ttoffZZ.constraint ="[[['b','W'],['Z']],[['b','W'],['Z']]]"
    #T6ttoffZZ.conditionDescription ="None"
    #T6ttoffZZ.condition ="None"
    #T6ttoffZZ.massConstraint = [['dm <= 175.0','dm >= 90.']]*2
    #T6ttoffZZ.dataUrl = None
    #T6ttoffZZ.source = 'SModelS'

    #+++++++ next mass plane block ++++++++++++++                                                                                                                             
    T6ttZZ_x050 = T6ttZZ.addMassPlane([[x, 0.50*x + (1-0.50)*y , y]]*2)    
    T6ttZZ_x050.addSource('efficiencyMap',"orig/T6ttZZ_EM/cms_sus_13_012_T6ttZZ_x05_EM_MAPS/MA5_EM_T6ttZZ_x05_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                               
    T6ttZZ_x070 = T6ttZZ.addMassPlane([[x, 0.70*x + (1-0.70)*y , y]]*2)
    T6ttZZ_x070.addSource('efficiencyMap',"orig/T6ttZZ_EM/cms_sus_13_012_T6ttZZ_x070_EM_MAPS/MA5_EM_T6ttZZ_x070_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6ttZZ_DeltaStopNeu200 = T6ttZZ.addMassPlane([[x, x-200.0 , y]]*2)
    T6ttZZ_DeltaStopNeu200.addSource('efficiencyMap',"orig/T6ttZZ_EM/cms_sus_13_012_T6ttZZ_DeltaStopNeu200_EM_MAPS/MA5_EM_T6ttZZ_DeltaStopNeu200_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6ttZZoff_DeltaNeu2LSP5 = T6ttZZoff.addMassPlane([[x, y+5.0 , y]]*2)
    T6ttZZoff_DeltaNeu2LSP5.addSource('efficiencyMap',"orig/T6ttZZ_EM/cms_sus_13_012_T6ttZZ_DeltaNeu2LSP5_EM_MAPS/MA5_EM_T6ttZZ_DeltaNeu2LSP5_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                               
    T6ttZZoff_DeltaNeu2LSP45 = T6ttZZoff.addMassPlane([[x, y+45.0 , y]]*2)
    T6ttZZoff_DeltaNeu2LSP45.addSource('efficiencyMap',"orig/T6ttZZ_EM/cms_sus_13_012_T6ttZZ_DeltaNeu2LSP45_EM_MAPS/MA5_EM_T6ttZZ_DeltaNeu2LSP45_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++
    T6ttZZoff_DeltaNeu2LSP85 = T6ttZZoff.addMassPlane([[x, y+85.0 , y]]*2)
    T6ttZZoff_DeltaNeu2LSP85.addSource('efficiencyMap',"orig/T6ttZZ_EM/cms_sus_13_012_T6ttZZ_DeltaNeu2LSP85_EM_MAPS/MA5_EM_T6ttZZ_DeltaNeu2LSP85_%s.dat" % dataset_n, "txt")
    #+++++++ next mass plane block ++++++++++++++                                                                                                                               
    T6ttZZ_DeltaNeu2LSP95 = T6ttZZ.addMassPlane([[x, y+95.0 , y]]*2)
    T6ttZZ_DeltaNeu2LSP95.addSource('efficiencyMap',"orig/T6ttZZ_EM/cms_sus_13_012_T6ttZZ_DeltaNeu2LSP95_EM_MAPS/MA5_EM_T6ttZZ_DeltaNeu2LSP95_%s.dat" % dataset_n, "txt")


    
    
    
datasets = {"SR_3NJet6_500HT800_450MHT600": ( 454, 418, 66 ), 
            "SR_3NJet6_1250HT1500_450MHTinf": ( 23, 17.6, 4.1 ),
            "SR_3NJet6_1250HT1500_300MHT450": ( 38, 42.8, 9.5 ),
            "SR_6NJet8_1000HT1250_200MHT300": ( 67, 70, 16 ),
            "SR_6NJet8_800HT1000_300MHT450": ( 35, 28.6, 6.9 ),
            "SR_6NJet8_500HT800_300MHT450": ( 62, 52, 12 ),
            "SR_3NJet6_800HT1000_200MHT300": ( 808, 777, 107 ),
            "SR_3NJet6_500HT800_200MHT300": ( 6159, 6088, 665 ),
            "SR_6NJet8_1250HT1500_450MHTinf": ( 2, 0.5, 2.6 ),
            "SR_3NJet6_1000HT1250_300MHT450": ( 129, 137, 20 ),
            "SR_8NJetinf_1250HT1500_200MHTinf": ( 5, 7.1, 3.8 ),
            "SR_8NJetinf_1000HT1250_200MHTinf": ( 8, 5.6, 2.3 ),
            "SR_3NJet6_800HT1000_300MHT450": ( 305, 330, 40 ),
            "SR_8NJetinf_800HT1000_200MHTinf": ( 9, 8.3, 3.4 ),
            "SR_6NJet8_1500HTinf_300MHTinf": ( 3, 7.9, 3.6 ),
            "SR_3NJet6_800HT1000_600MHTinf": ( 52, 54.8, 9.7 ),
            "SR_6NJet8_500HT800_450MHTinf": ( 9, 0.8, 3.3 ),
            "SR_6NJet8_1000HT1250_300MHT450": ( 20, 21.6, 5.8 ),
            "SR_6NJet8_800HT1000_450MHTinf": ( 4, 6.0, 2.8 ),
            "SR_3NJet6_1500HTinf_200MHT300": ( 94, 86, 17 ),
            "SR_3NJet6_500HT800_600MHTinf": ( 62, 57.4, 11.2 ),
            "SR_3NJet6_1250HT1500_200MHT300": ( 98, 109, 18 ),
            "SR_3NJet6_1000HT1250_600MHTinf": ( 32, 22.8, 5.2 ),
            "SR_3NJet6_1000HT1250_200MHT300": ( 335, 305, 41 ),
            "SR_8NJetinf_1500HTinf_200MHTinf": ( 2, 3.3, 4.7 ),
            "SR_6NJet8_1500HTinf_200MHT300": ( 18, 21.1, 8.1 ),
            "SR_6NJet8_1250HT1500_200MHT300": ( 24, 28., 8.2 ),
            "SR_3NJet6_800HT1000_450MHT600": ( 124, 108, 15 ),
            "SR_6NJet8_800HT1000_200MHT300": ( 111, 124, 29 ),
            "SR_6NJet8_500HT800_200MHT300": ( 266, 290, 65 ),
            "SR_3NJet6_500HT800_300MHT450": ( 2305, 2278, 266 ),
            "SR_8NJetinf_500HT800_200MHTinf": ( 8, 4.8, 2.3 ),
            "SR_3NJet6_1000HT1250_450MHT600": ( 34, 32.3, 6.1 ),
            "SR_6NJet8_1000HT1250_450MHTinf": ( 4, 2.2, 3.8 ),
            "SR_6NJet8_1250HT1500_300MHT450": ( 5, 9.4, 3.6 ),
            "SR_3NJet6_1500HTinf_300MHTinf": ( 39, 29.7, 5.8 )
}

for name, numbers in datasets.items():
    #+++++++ dataset block ++++++++++++++
    dataset = DataSetInput( name )
    name = name.replace('SR_','')
    dataset.setInfo(dataType = 'efficiencyMap', dataId = name, 
            observedN = numbers[0], expectedBG = numbers[1], bgError = numbers[2] )
    add ( dataset )

databaseCreator.create()
