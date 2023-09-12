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
info = MetaInfoInput('CMS-SUS-12-024')
info.sqrts = '8.0'
info.private = False
info.lumi = '19.4'
info.publication = 'http://www.sciencedirect.com/science/article/pii/S0370269313005339'
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS12024'
info.arxiv = 'http://arxiv.org/abs/1305.2390'
info.contact = 'Keith Ulmer <keith.ulmer@cern.ch>, Josh Thompson <joshua.thompson@cern.ch>, Alessandro Gaz <alessandro.gaz@cern.ch>'
info.prettyName = '0 leptons + >= 3 (1b-)jets + Etmiss'




#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET3_HT3_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET3_HT3_nb3", observedN = 6, expectedBG = 5.9 , bgError = 1.9, upperLimit = '3.9702E-01*fb', expectedUpperLimit = '3.5007E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET3_HT3_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET3_HT3_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET4_HT3_nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET4_HT3_nb2", observedN = 19, expectedBG = 20.7 , bgError = 3.2, upperLimit = '5.7300E-01*fb', expectedUpperLimit = '6.1033E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET4_HT3_nb2")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET4_HT3_nb2")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET2_HT2_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET2_HT2_nb3", observedN = 182, expectedBG = 179 , bgError = 13, upperLimit = '2.0613E+00*fb', expectedUpperLimit = '1.9515E+00*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET2_HT2_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET2_HT2_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET2_HT1_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET2_HT1_nb3", observedN = 161, expectedBG = 157 , bgError = 13, upperLimit = '2.0431E+00*fb', expectedUpperLimit = '1.8875E+00*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET2_HT1_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET2_HT1_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET4_HT2_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET4_HT2_nb3", observedN = 8, expectedBG = 8.4 , bgError = 2.1, upperLimit = '4.2514E-01*fb', expectedUpperLimit = '4.2430E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET4_HT2_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET4_HT2_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET4_HT4_nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET4_HT4_nb2", observedN = 19, expectedBG = 19.0 , bgError = 3.2, upperLimit = '6.2027E-01*fb', expectedUpperLimit = '6.2027E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET4_HT4_nb2")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET4_HT4_nb2")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET4_HT4_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET4_HT4_nb3", observedN = 4, expectedBG = 2.1 , bgError = 1.1, upperLimit = '3.8500E-01*fb', expectedUpperLimit = '2.6167E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET4_HT4_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET4_HT4_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET3_HT4_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET3_HT4_nb3", observedN = 4, expectedBG = 2.9 , bgError = 1.3, upperLimit = '3.6165E-01*fb', expectedUpperLimit = '2.5204E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET3_HT4_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET3_HT4_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET2_HT3_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET2_HT3_nb3", observedN = 18, expectedBG = 23.2 , bgError = 3.8, upperLimit = '5.0635E-01*fb', expectedUpperLimit = '6.8517E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET2_HT3_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET2_HT3_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET4_HT2_nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET4_HT2_nb2", observedN = 66, expectedBG = 70.5 , bgError = 6.3, upperLimit = '9.7444E-01*fb', expectedUpperLimit = '1.1161E+00*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET4_HT2_nb2")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET4_HT2_nb2")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET3_HT1_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET3_HT1_nb3", observedN = 15, expectedBG = 15.5 , bgError = 3.0, upperLimit = '5.5542E-01*fb', expectedUpperLimit = '5.5753E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET3_HT1_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET3_HT1_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET3_HT2_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET3_HT2_nb3", observedN = 36, expectedBG = 32.1 , bgError = 4.3, upperLimit = '9.6536E-01*fb', expectedUpperLimit = '7.8674E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET3_HT2_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET3_HT2_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET2_HT4_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET2_HT4_nb3", observedN = 14, expectedBG = 12.3 , bgError = 2.7, upperLimit = '6.0499E-01*fb', expectedUpperLimit = '5.1236E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET2_HT4_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET2_HT4_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET4_HT3_nb3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET4_HT3_nb3", observedN = 2, expectedBG = 2.0 , bgError = 1.0, upperLimit = '2.6247E-01*fb', expectedUpperLimit = '2.6247E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="AL"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap', "orig/efficiency_T1bbbb_multi.root", "root", objectName = "heff_MET4_HT3_nb3")
T1bbbb.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1bbbb.root", "root", objectName = "graph_smoothed_Obs_T1bbbb")
T1bbbb.figure = "Fig 7 (left)"
T1bbbb.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1bbbb_exclusions_corrected.pdf"
T1bbbb.dataUrl  = None
# +++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked ="AL"
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.condition ="None"
T1ttttoff.conditionDescription ="None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"orig/efficiency_T1tttt_multi.root", "root", objectName = "heff_MET4_HT3_nb3")
T1tttt.addSource('obsExclusion',"orig/graph_smoothed_Obs_T1tttt.root", "root", objectName = "graph_smoothed_Obs_T1tttt")
T1tttt.figure = "Fig 7 (right)"
T1tttt.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024/T1tttt_exclusions_corrected.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12024"
T1ttttoff.addMassPlane(T1tttt)

databaseCreator.create()
