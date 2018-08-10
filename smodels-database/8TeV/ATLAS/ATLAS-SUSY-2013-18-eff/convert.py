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
info = MetaInfoInput('ATLAS-SUSY-2013-18')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/'
info.sqrts = '8.0*TeV'
info.lumi = 20.1
info.prettyName = '0-1 leptons + >= 3 b-jets + Etmiss'
info.private =False
info.arxiv = 'http://arxiv.org/abs/1407.0600'
info.publication ='http://link.springer.com/article/10.1007/JHEP10(2014)024'
info.comment ='Using combined Exclusion line. Superseeded CONF-2013-061(from FastLim) has more topologies'
info.supersedes ='ATLAS-CONF-2013-061'



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR-0l-4j-B")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR-0l-4j-B", observedN = 3, expectedBG = 1.3 , bgError = 0.9, upperLimit = '3.3184E-01*fb', expectedUpperLimit = '2.0669E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.constraint = "[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb_1.figure  = "Figure 5b and Figure 6b"
T1bbbb_1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_05b.png"
T1bbbb_1.addSource('obsExclusion', 'orig/T1bbbb_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'txt')
T1bbbb_1.addSource('efficiencyMap','orig/EffMap_T1bbbb_SR-0l-4j-B.txt', 'txt')
T1bbbb_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304457/all"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR-0l-4j-C")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR-0l-4j-C", observedN = 1, expectedBG = 1.6 , bgError = 0.7, upperLimit = '2.0051E-01*fb', expectedUpperLimit = '2.0054E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.constraint = "[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb_1.figure  = "Figure 5c and Figure 6c"
T1bbbb_1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_05c.png"
T1bbbb_1.addSource('obsExclusion', 'orig/T1bbbb_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'txt')
T1bbbb_1.addSource('efficiencyMap','orig/EffMap_T1bbbb_SR-0l-4j-C.txt', 'txt')
T1bbbb_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304457/all"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR-0l-4j-A")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR-0l-4j-A", observedN = 2, expectedBG = 1.6 , bgError = 0.9, upperLimit = '2.6218E-01*fb', expectedUpperLimit = '2.0194E-01*fb')
#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.constraint = "[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb_1.figure  = "Figure 5a and Figure 6a"
T1bbbb_1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_05a.png"
T1bbbb_1.addSource('obsExclusion', 'orig/T1bbbb_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'txt')
T1bbbb_1.addSource('efficiencyMap','orig/EffMap_T1bbbb_SR-0l-4j-A.txt', 'txt')
T1bbbb_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304457/all"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR-0l-7j-C")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR-0l-7j-C", observedN = 1, expectedBG = 0.9 , bgError = 1, upperLimit = '2.0994E-01*fb', expectedUpperLimit = '1.4848E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "Figure 14c and Figure 15c"
T1tttt_1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_14c.png"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'txt')
T1tttt_1.addSource('efficiencyMap','orig/EffMap_T1tttt_SR-0l-7j-C.txt', 'txt')
T1tttt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304457/all"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR-0l-7j-B")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR-0l-7j-B", observedN = 3, expectedBG = 3.2 , bgError = 1.6, upperLimit = '2.9336E-01*fb', expectedUpperLimit = '2.9436E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "Figure 14b and Figure 15b"
T1tttt_1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_14b.png"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'txt')
T1tttt_1.addSource('efficiencyMap','orig/EffMap_T1tttt_SR-0l-7j-B.txt', 'txt')
T1tttt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304457/all"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR-0l-7j-A")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR-0l-7j-A", observedN = 21, expectedBG = 21.2 , bgError = 4.6, upperLimit = '6.9135E-01*fb', expectedUpperLimit = '6.9189E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "Figure 13b"
T1tttt_1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_14a.png"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'txt')
T1tttt_1.addSource('efficiencyMap','orig/EffMap_T1tttt_SR-0l-7j-A.txt', 'txt')
T1tttt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304457/d43/input http://hepdata.cedar.ac.uk/view/ins1304457/d44/input"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR-1l-6j-A")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR-1l-6j-A", observedN = 7, expectedBG = 13.5 , bgError = 3.2, upperLimit = '3.1818E-01*fb', expectedUpperLimit = '5.2681E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "Figure 14d and Figure 15d"
T1tttt_1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_14d.png"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'txt')
T1tttt_1.addSource('efficiencyMap','orig/EffMap_T1tttt_SR-1l-6j-A.txt', 'txt')
T1tttt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304457/all"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR-1l-6j-C")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR-1l-6j-C", observedN = 0, expectedBG = 2.3 , bgError = 0.7, upperLimit = '1.5160E-01*fb', expectedUpperLimit = '2.4480E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "Figure 14f and Figure 15f"
T1tttt_1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_14f.png"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'txt')
T1tttt_1.addSource('efficiencyMap','orig/EffMap_T1tttt_SR-1l-6j-C.txt', 'txt')
T1tttt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304457/all"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR-1l-6j-B")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR-1l-6j-B", observedN = 0, expectedBG = 6.1 , bgError = 1.8, upperLimit = '1.3940E-01*fb', expectedUpperLimit = '3.7710E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure  = "Figure 14e and Figure 15e"
T1tttt_1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-18/figaux_14e.png"
T1tttt_1.addSource('obsExclusion', 'orig/T1tttt_ATLAS-SUSY-2013-18-Obs_Excl.dat', 'txt')
T1tttt_1.addSource('efficiencyMap','orig/EffMap_T1tttt_SR-1l-6j-B.txt', 'txt')
T1tttt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304457/all"

databaseCreator.create()
