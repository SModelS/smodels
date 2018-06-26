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
info = MetaInfoInput('CMS-PAS-EXO-16-036')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/EXO-16-036/index.html'
info.sqrts = 13
info.lumi = 12.9
info.prettyName ='hscp search'
info.private = False
info.contact ='Andre Lessa <lessa.a.p@gmail.com>; Jan Heisig <heisig@physik.rwth-aachen.de>; SModelS'
info.comment ='Search for long-lived charged particles implemented using a correction factor for the 8 TeV efficiencies.'

#++Define list of datasets++
datasetNames = ['c000','c100','c200','c300']
observedNs = [5,1,0,0] #Extracted from a fit of the data
expectedBGs = [2.63,0.3377,0.0258,0.0045] #Extracted from a log fit of CMS BG
bgErrors = [0.53,0.127,0.0036,0.001] #Estimated from the error between a log fit and a linear fit of the CMS BG
obsUpperLimits = ['0.625*fb','0.346*fb','0.232*fb','0.233*fb']
expUpperLimits = ['0.36*fb','0.232*fb','0.233*fb','0.232*fb']
for i,name in enumerate(datasetNames):
#+++++++ dataset block ++++++++++++++
    dataset = DataSetInput(name)
    dataset.setInfo(dataType = 'efficiencyMap', dataId = name, 
                    observedN=observedNs[i], expectedBG=expectedBGs[i], bgError=bgErrors[i],
					upperLimit = obsUpperLimits[i], expectedUpperLimit = expUpperLimits[i])

    #+++++++ txnames ++++++++++++++++++++
    #+++++++ next txName block ++++++++++++++
    HSCPM1 = dataset.addTxName('THSCPM1')
    HSCPM1.checked =''
    HSCPM1.constraint = "[[],[]]"
    HSCPM1.condition =None
    HSCPM1.finalState = ['HSCP','HSCP']
    HSCPM1.massConstraints = None
    HSCPM1.dataUrl = None
    HSCPM1.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++
    plane = HSCPM1.addMassPlane([[x],[x]])
    plane.setSources(dataLabels= ['efficiencyMap'],dataFiles=['orig/effmap_M1_chargino_cons_mre'+name+'_clean.txt'], dataFormats=['txt'])
#    plane.addSource(dataLabels='obsExclusion',dataFiles='orig/CMS-PAS-EXO-16-036_Figure_003-b.dat', dataFormats='txt', unit='pb')
    #+++++++ next txName block ++++++++++++++
    HSCPM3 = dataset.addTxName('THSCPM3')
    HSCPM3.checked =''
    HSCPM3.constraint = "[[['*']],[['*']]]"  ##Here '*' represents any (single) even particle
    HSCPM3.condition =None
    HSCPM3.conditionDescription =None
    HSCPM3.finalState = ['HSCP','HSCP']
    HSCPM3.massConstraints = None
    HSCPM3.dataUrl = None
    HSCPM3.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++
    plane = HSCPM3.addMassPlane([[x,y],[x,y]])
    plane.setSources(dataLabels= ['efficiencyMap'],dataFiles=['orig/effmap_M3_chargino_cons_mre'+name+'_clean.txt'], dataFormats=['txt'])
    #+++++++ next txName block ++++++++++++++
    HSCPM5 = dataset.addTxName('THSCPM5')
    HSCPM5.checked =''
    HSCPM5.constraint = "[[['*'],['*']],[['*'],['*']]]" ##Here '*' represents any (single) even particle
    HSCPM5.condition =None
    HSCPM5.conditionDescription =None
    HSCPM5.finalState = ['HSCP','HSCP']
    HSCPM5.massConstraints = None
    HSCPM5.dataUrl = None
    HSCPM5.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++
    plane = HSCPM5.addMassPlane([[x,y,z]]*2)
    plane.setSources(dataLabels= ['efficiencyMap'],dataFiles=['orig/effmap_M5_stau_cons_mre'+name+'_clean.txt'], dataFormats=['txt'])
    #+++++++ next txName block ++++++++++++++
    HSCPM7 = dataset.addTxName('THSCPM7')
    HSCPM7.checked =''
    HSCPM7.constraint = "[[['*']],[['*'],['*']]]" ##Here '*' represents any (single) even particle
    HSCPM7.condition =None
    HSCPM7.conditionDescription =None
    HSCPM7.finalState = ['HSCP','HSCP']
    HSCPM7.massConstraints = None
    HSCPM7.dataUrl = None
    HSCPM7.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++
    plane = HSCPM7.addMassPlane([[x,z],[x,y,z]])
    plane.setSources(dataLabels= ['efficiencyMap'],dataFiles=['orig/effmap_M7_stau_cons_mre'+name+'_clean.txt'], dataFormats=['txt'])
    #+++++++ next txName block ++++++++++++++
    HSCPM8 = dataset.addTxName('THSCPM8')
    HSCPM8.checked =''
    HSCPM8.constraint = "[[['*','*']],[['*','*']]]" ##Here '*','*' represents any pair of even particles
    HSCPM8.condition =None
    HSCPM8.conditionDescription =None
    HSCPM8.finalState = ['HSCP','HSCP']
    HSCPM8.massConstraints = None
    HSCPM8.dataUrl = None
    HSCPM8.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++
    plane = HSCPM8.addMassPlane([[x,y],[x,y]])
    plane.setSources(dataLabels= ['efficiencyMap'],dataFiles=['orig/effmap_M8_stau_cons_mre'+name+'_clean.txt'], dataFormats=['txt'])
    #+++++++ next txName block ++++++++++++++
    HSCPM2 = dataset.addTxName('THSCPM2')
    HSCPM2.checked =''
    HSCPM2.constraint = "[[*],[]]" ##Here [*] represents one branch with any list of vertices
    HSCPM2.condition =None
    HSCPM2.conditionDescription =None
    HSCPM2.finalState = ['MET','HSCP']
    HSCPM2.massConstraints = None
    HSCPM2.dataUrl = None
    HSCPM2.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++
    plane = HSCPM2.addMassPlane([['*'],[x]]) ##Here ['*'] represents a mass array with any length
    plane.setSources(dataLabels= ['efficiencyMap'],dataFiles=['orig/effmap_M2_chargino_cons_mre'+name+'_clean.txt'], dataFormats=['txt'])
    #+++++++ next txName block ++++++++++++++
    HSCPM4 = dataset.addTxName('THSCPM4')
    HSCPM4.checked =''
    HSCPM4.constraint = "[[*],[['*']]]" ##Here [*] represents one branch with any list of vertices
    HSCPM4.condition =None
    HSCPM4.conditionDescription =None
    HSCPM4.finalState = ['MET','HSCP']
    HSCPM4.massConstraints = None
    HSCPM4.dataUrl = None
    HSCPM4.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++
    plane = HSCPM4.addMassPlane([['*'],[x,y]]) ##Here ['*'] represents a mass array with any length
    plane.setSources(dataLabels= ['efficiencyMap'],dataFiles=['orig/effmap_M4_chargino_cons_mre'+name+'_clean.txt'], dataFormats=['txt'])
    #+++++++ next txName block ++++++++++++++
    HSCPM6 = dataset.addTxName('THSCPM6')
    HSCPM6.checked =''
    HSCPM6.constraint = "[[*],[['*'],[ '*']]]" ##Here [*] represents one branch with any list of vertices and '*' any single particle
    HSCPM6.condition =None
    HSCPM6.conditionDescription =None
    HSCPM6.finalState = ['MET','HSCP']
    HSCPM6.massConstraints = None
    HSCPM6.dataUrl = None
    HSCPM6.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++
    plane = HSCPM6.addMassPlane([['*'],[x,y,z]]) ##Here ['*'] represents a mass array with any length
    plane.setSources(dataLabels= ['efficiencyMap'],dataFiles=['orig/effmap_M6_stau_cons_mre'+name+'_clean.txt'], dataFormats=['txt'])


databaseCreator.create()

