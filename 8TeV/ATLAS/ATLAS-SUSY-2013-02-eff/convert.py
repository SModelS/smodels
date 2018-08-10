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
info = MetaInfoInput('ATLAS-SUSY-2013-02')
info.comment = 'T5WWLSP060 and T6WWLSP060 originally have xvalue on y-axes, changed by us to M2'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication = 'http://link.springer.com/article/10.1007/JHEP09%282014%29176'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/'
info.arxiv = 'http://arxiv.org/abs/1405.7875'
info.prettyName = '0 leptons + 2-6 jets + Etmiss'
info.supersedes = 'ATLAS-CONF-2013-047'



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR6jl")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR6jl", observedN = 121, expectedBG = 111 , bgError = 11, upperLimit = '1.9230E+00*fb', expectedUpperLimit = '1.5312E+00*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_26b"
T2.addSource('efficiencyMap',"orig/T2_SR6jl.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_26b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d139"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_42b"
T1.addSource('efficiencyMap',"orig/T1_SR6jl.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_42b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d185"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR6jm")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR6jm", observedN = 39, expectedBG = 33 , bgError = 6, upperLimit = '1.1173E+00*fb', expectedUpperLimit = '8.6116E-01*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_27b"
T2.addSource('efficiencyMap',"orig/T2_SR6jm.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_27b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d142"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_43b"
T1.addSource('efficiencyMap',"orig/T1_SR6jm.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_43b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d188"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR6jt+")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR6jt+", observedN = 6, expectedBG = 4.9 , bgError = 1.6, upperLimit = '3.9922E-01*fb', expectedUpperLimit = '3.0218E-01*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_29b"
T2.addSource('efficiencyMap',"orig/T2_SR6jt+.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_29b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d148"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_45b"
T1.addSource('efficiencyMap',"orig/T1_SR6jt+.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_45b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d194"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jt", observedN = 133, expectedBG = 125 , bgError = 10, upperLimit = '1.8181E+00*fb', expectedUpperLimit = '1.5124E+00*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_17b"
T2.addSource('efficiencyMap',"orig/T2_SR2jt.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_17b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d112"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_33b"
T1.addSource('efficiencyMap',"orig/T1_SR2jt.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_33b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d158"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR5j")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR5j", observedN = 121, expectedBG = 126 , bgError = 13, upperLimit = '1.5429E+00*fb', expectedUpperLimit = '1.7138E+00*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_25b"
T2.addSource('efficiencyMap',"orig/T2_SR5j.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_25b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d136"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_41b"
T1.addSource('efficiencyMap',"orig/T1_SR5j.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_41b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d182"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jW")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jW", observedN = 0, expectedBG = 2.3 , bgError = 1.4, upperLimit = '1.4709E-01*fb', expectedUpperLimit = '2.5070E-01*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_18b"
T2.addSource('efficiencyMap',"orig/T2_SR2jW.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_18b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d115"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_34b"
T1.addSource('efficiencyMap',"orig/T1_SR2jW.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_34b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d161"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR4jW")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jW", observedN = 16, expectedBG = 14 , bgError = 4, upperLimit = '6.7961E-01*fb', expectedUpperLimit = '5.9339E-01*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_20b"
T2.addSource('efficiencyMap',"orig/T2_SR4jW.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_20b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d121"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_36b"
T1.addSource('efficiencyMap',"orig/T1_SR4jW.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_36b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d167"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR4jt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jt", observedN = 0, expectedBG = 2.5 , bgError = 1.0, upperLimit = '1.4949E-01*fb', expectedUpperLimit = '2.4033E-01*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_24b"
T2.addSource('efficiencyMap',"orig/T2_SR4jt.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_24b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d133"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_40b"
T1.addSource('efficiencyMap',"orig/T1_SR4jt.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_40b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d179"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jl")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jl", observedN = 12315, expectedBG = 13000 , bgError = 1000, upperLimit = '7.7800E+01*fb', expectedUpperLimit = '9.7112E+01*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_15b"
T2.addSource('efficiencyMap',"orig/T2_SR2jl.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_15b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d106"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_31b"
T1.addSource('efficiencyMap',"orig/T1_SR2jl.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_31b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d152"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jm")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jm", observedN = 715, expectedBG = 760 , bgError = 50, upperLimit = '4.2419E+00*fb', expectedUpperLimit = '5.5524E+00*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_16b"
T2.addSource('efficiencyMap',"orig/T2_SR2jm.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_16b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d109"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_32b"
T1.addSource('efficiencyMap',"orig/T1_SR2jm.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_32b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d155"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR4jl-")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jl-", observedN = 2169, expectedBG = 2120 , bgError = 110, upperLimit = '1.3292E+01*fb', expectedUpperLimit = '1.1561E+01*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_21b"
T2.addSource('efficiencyMap',"orig/T2_SR4jl-.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_21b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d124"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_37b"
T1.addSource('efficiencyMap',"orig/T1_SR4jl-.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_37b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d170"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR4jl")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jl", observedN = 608, expectedBG = 630 , bgError = 50, upperLimit = '4.7487E+00*fb', expectedUpperLimit = '5.4345E+00*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_22b"
T2.addSource('efficiencyMap',"orig/T2_SR4jl.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_22b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d127"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_38b"
T1.addSource('efficiencyMap',"orig/T1_SR4jl.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_38b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d173"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR4jm")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jm", observedN = 24, expectedBG = 37 , bgError = 6, upperLimit = '5.0301E-01*fb', expectedUpperLimit = '8.8617E-01*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_23b"
T2.addSource('efficiencyMap',"orig/T2_SR4jm.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_23b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d130"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_39b"
T1.addSource('efficiencyMap',"orig/T1_SR4jm.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_39b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d176"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR3j")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR3j", observedN = 7, expectedBG = 5 , bgError = 1.2, upperLimit = '4.3344E-01*fb', expectedUpperLimit = '3.3172E-01*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_19b"
T2.addSource('efficiencyMap',"orig/T2_SR3j.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_19b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d118"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_35b"
T1.addSource('efficiencyMap',"orig/T1_SR3j.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_35b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d164"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR6jt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR6jt", observedN = 5, expectedBG = 5.2 , bgError = 1.4, upperLimit = '3.3159E-01*fb', expectedUpperLimit = '3.3330E-01*fb')
#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane([[x,y]]*2)
T2.figure  = "figaux_28b"
T2.addSource('efficiencyMap',"orig/T2_SR6jt.dat","txt")
T2.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
T2.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_28b.png"
T2.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d145"
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_44b"
T1.addSource('efficiencyMap',"orig/T1_SR6jt.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_44b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d191"

databaseCreator.create()
