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
type = str)
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = str)
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
info = MetaInfoInput('CMS-SUS-13-007')
info.url ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13007'
info.sqrts = 8
info.lumi = 19.3
info.prettyName = '1 lepton + >= 2 b-jets + Etmiss'
info.private = False
info.arxiv = 'http://arxiv.org/abs/arXiv:1311.4937'
info.contact ='Loukas Gouskos <loukas.gouskos@cern.ch>, Markus Stoye <Markus.Stoye@cern.ch>'
info.publication = 'http://www.sciencedirect.com/science/article/pii/S037026931400255X'
info.comment = 'Only two mass planes for T5tttt; implemented Delta Phi method'
info.implementedBy = 'Federico Ambrogi'




#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mu350Nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mu350Nb2", observedN = 2, expectedBG = 1.4 , bgError = 1.17, upperLimit = '2.764E-01*fb', expectedUpperLimit = '2.139E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Muons_350_Nb2.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Muons_350_Nb2.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_mu-350-sig-Nb2')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mu350Nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mu350Nb3", observedN = 0, expectedBG = 0.6 , bgError = 0.58, upperLimit = '1.550E-01*fb', expectedUpperLimit = '1.548E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Muons_350_Nb3.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Muons_350_Nb3.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_mu-350-sig-Nb3')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mu250Nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mu250Nb3", observedN = 1e-05, expectedBG = 1.9 , bgError = 0.89, upperLimit = '1.556E-01*fb', expectedUpperLimit = '2.080E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Muons_250_Nb3.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Muons_250_Nb3.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_mu-250-sig-Nb3')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mu250Nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mu250Nb2", observedN = 9, expectedBG = 6 , bgError = 2.37, upperLimit = '5.748E-01*fb', expectedUpperLimit = '4.164E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Muons_250_Nb2.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Muons_250_Nb2.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_mu-250-sig-Nb2')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("el450Nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "el450Nb2", observedN = 0, expectedBG = 0 , bgError = 0.44, upperLimit = '1.550E-01*fb', expectedUpperLimit = '1.550E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Elect_450_Nb2.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Elect_450_Nb2.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_el-450-sig-Nb2')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("el450Nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "el450Nb3", observedN = 0, expectedBG = 0 , bgError = 0.1, upperLimit = '1.547E-01*fb', expectedUpperLimit = '1.547E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Elect_450_Nb3.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Elect_450_Nb3.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_el-450-sig-Nb3')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mu450Nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mu450Nb3", observedN = 0, expectedBG = 0 , bgError = 0.22, upperLimit = '1.549E-01*fb', expectedUpperLimit = '1.549E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Muons_450_Nb3.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Muons_450_Nb3.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_mu-450-sig-Nb3')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mu450Nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mu450Nb2", observedN = 0, expectedBG = 0 , bgError = 0.72, upperLimit = '1.553E-01*fb', expectedUpperLimit = '1.553E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Muons_450_Nb2.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Muons_250_Nb2.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_mu-450-sig-Nb2')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("el350Nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "el350Nb3", observedN = 0, expectedBG = 0.9 , bgError = 0.8, upperLimit = '1.547E-01*fb', expectedUpperLimit = '1.552E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Elect_350_Nb3.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Elect_350_Nb3.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_el-350-sig-Nb3')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("el350Nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "el350Nb2", observedN = 2, expectedBG = 2.7 , bgError = 2.06, upperLimit = '2.637E-01*fb', expectedUpperLimit = '2.639E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Elect_350_Nb2.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Elect_350_Nb2.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_el-350-sig-Nb2')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("el250Nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "el250Nb2", observedN = 9, expectedBG = 3.8 , bgError = 1.89, upperLimit = '6.484E-01*fb', expectedUpperLimit = '3.006E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Elect_250_Nb2.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Elect_250_Nb2.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_el-250-sig-Nb2')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("el250Nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "el250Nb3", observedN = 4, expectedBG = 1.9 , bgError = 0.98, upperLimit = '3.923E-01*fb', expectedUpperLimit = '2.083E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "DPhi_Method_T1tttt_Elect_250_Nb3.pdf"
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/DPhi_Method_T1tttt_Elect_250_Nb3.pdf"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_excl.txt', 'txt' )
T1tttt_1.addSource('efficiencyMap', 'orig/EffxAcc_Smoothed_DPhi_Paper_Dec2013.root', 'root', objectName = 'EffxAcc-ISR-T1tttt-EffxAcc-ISR-MG_el-250-sig-Nb3')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/edit/CMSPublic/EffxAcc?topicparent=CMSPublic.PhysicsResultsSUS13007;nowysiwyg=1'
T1ttttoff.addMassPlane(T1tttt_1)

databaseCreator.create()
