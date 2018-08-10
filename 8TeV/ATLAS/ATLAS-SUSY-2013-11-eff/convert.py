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

DataSetInput.ntoys=args.ntoys

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2013-11')
info.comment = 'TSlepSlep,TChiWW and TChipChimSlepSnu efficiency maps created by the SModelS collaboration using MadAnalysis5'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication = 'http://link.springer.com/article/10.1007/JHEP05(2014)071'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-11/'
info.arxiv = 'http://arxiv.org/abs/1403.5294'
info.contact = 'SModelS'
info.prettyName = '2 leptons (e,mu) + Etmiss'
info.supersedes = 'ATLAS-CONF-2013-049'



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("WWb-SF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "WWb-SF", observedN = 26, expectedBG = 30.2 , bgError = 3.5)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_WWb-SF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_WWb-SF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_WWb-SF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_WWb-SF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_WWb-SF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_WWb-SF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_WWb-SF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_WWb-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_WWb-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_WWb-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_WWb-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_WWb-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_WWb-SF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("WWa-SF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "WWa-SF", observedN = 73, expectedBG = 86.5 , bgError = 7.4)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_WWa-SF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_WWa-SF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_WWa-SF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_WWa-SF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_WWa-SF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_WWa-SF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_WWa-SF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_WWa-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_WWa-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_WWa-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_WWa-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_WWa-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_WWa-SF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("WWb-DF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "WWb-DF", observedN = 17, expectedBG = 18.1 , bgError = 2.6)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_WWb-DF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_WWb-DF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_WWb-DF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_WWb-DF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_WWb-DF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_WWb-DF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_WWb-DF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_WWb-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_WWb-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_WWb-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_WWb-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_WWb-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_WWb-DF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mT2-90-DF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mT2-90-DF", observedN = 21, expectedBG = 23.3 , bgError = 3.7)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_mT2-90-DF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_mT2-90-DF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_mT2-90-DF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_mT2-90-DF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_mT2-90-DF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_mT2-90-DF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_mT2-90-DF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_mT2-90-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_mT2-90-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_mT2-90-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_mT2-90-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_mT2-90-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_mT2-90-DF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mT2-120-SF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mT2-120-SF", observedN = 5, expectedBG = 8.9 , bgError = 2.1)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_mT2-120-SF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_mT2-120-SF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_mT2-120-SF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_mT2-120-SF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_mT2-120-SF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_mT2-120-SF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_mT2-120-SF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_mT2-120-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_mT2-120-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_mT2-120-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_mT2-120-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_mT2-120-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_mT2-120-SF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("WWc-SF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "WWc-SF", observedN = 10, expectedBG = 20.3 , bgError = 3.5)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_WWc-SF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_WWc-SF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_WWc-SF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_WWc-SF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_WWc-SF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_WWc-SF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_WWc-SF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_WWc-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_WWc-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_WWc-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_WWc-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_WWc-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_WWc-SF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("WWa-DF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "WWa-DF", observedN = 70, expectedBG = 73.6 , bgError = 7.9)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_WWa-DF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_WWa-DF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_WWa-DF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_WWa-DF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_WWa-DF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_WWa-DF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_WWa-DF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_WWa-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_WWa-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_WWa-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_WWa-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_WWa-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_WWa-DF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mT2-90-SF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mT2-90-SF", observedN = 33, expectedBG = 38.2 , bgError = 5.1)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_mT2-90-SF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_mT2-90-SF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_mT2-90-SF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_mT2-90-SF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_mT2-90-SF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_mT2-90-SF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_mT2-90-SF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_mT2-90-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_mT2-90-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_mT2-90-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_mT2-90-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_mT2-90-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_mT2-90-SF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("WWc-DF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "WWc-DF", observedN = 11, expectedBG = 9.0 , bgError = 2.2)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_WWc-DF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_WWc-DF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_WWc-DF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_WWc-DF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_WWc-DF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_WWc-DF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_WWc-DF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_WWc-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_WWc-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_WWc-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_WWc-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_WWc-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_WWc-DF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mT2-120-DF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mT2-120-DF", observedN = 5, expectedBG = 3.6 , bgError = 1.2)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_mT2-120-DF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_mT2-120-DF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_mT2-120-DF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_mT2-120-DF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_mT2-120-DF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_mT2-120-DF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_mT2-120-DF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_mT2-120-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_mT2-120-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_mT2-120-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_mT2-120-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_mT2-120-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_mT2-120-DF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("Zjets")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "Zjets", observedN = 1, expectedBG = 1.4 , bgError = 0.6)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_Zjets.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_Zjets.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_Zjets.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_Zjets.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_Zjets.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_Zjets.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_Zjets.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_Zjets.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_Zjets.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_Zjets.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_Zjets.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_Zjets.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_Zjets.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mT2-150-SF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mT2-150-SF", observedN = 3, expectedBG = 3.2 , bgError = 0.7)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_mT2-150-SF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_mT2-150-SF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_mT2-150-SF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_mT2-150-SF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_mT2-150-SF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_mT2-150-SF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_mT2-150-SF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_mT2-150-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_mT2-150-SF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_mT2-150-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_mT2-150-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_mT2-150-SF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_mT2-150-SF.dat", "txt")


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mT2-150-DF")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mT2-150-DF", observedN = 2, expectedBG = 1.0 , bgError = 0.5)
#+++++++ next txName block ++++++++++++++
#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane([[x,y]]*2)
TSlepSlep_1.addSource('obsExclusion', "orig/exclusion_TSlepSlep.txt", "txt")
TSlepSlep_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TSlepSlep_1_EM_MAPS/MA5_EM_TSlepSlep_1_mT2-150-DF.dat', 'txt')
TSlepSlep_1.dataUrl = None
TChiWW = dataset.addTxName('TChiWW')
TChiWW.constraint ="[[['W+']],[['W-']]]"
TChiWW.massConstraint = None
TChiWW.conditionDescription = None
TChiWW.condition = None
TChiWW.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWW_1 = TChiWW.addMassPlane([[x,y]]*2)
TChiWW_1.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChiWW_1_EM_MAPS/MA5_EM_TChiWW_1_mT2-150-DF.dat', 'txt')
TChiWW_1.addSource('obsExclusion', "orig/exclusion_TChiWW.txt", "txt")
TChiWW_1.dataUrl = None
TChiWWoff = dataset.addTxName('TChiWWoff')
TChiWWoff.constraint = "22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.conditionDescription = "[[['l+','nu']],[['l-','nu']]] > 2* [[['mu+','nu']],[['l-','nu']]]"
TChiWWoff.condition = "Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['mu+','nu']],[['l-','nu']]]); Cgtr([[['l+','nu']],[['l-','nu']]],2.*[[['l+','nu']],[['mu-','nu']]])"
constraint ="22.2*[[['l+','nu']],[['l-','nu']]]"
TChiWWoff.massConstraint = [['dm <= 76.']]*2
TChiWWoff.source = 'SModelS'
#+++++++ next mass plane block ++++++++++++++
TChiWWoff.addMassPlane( TChiWW_1 )
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.constraint = "[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription = "[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition = "Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = 'SModelS'
TChipChimSlepSnu.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_x005 = TChipChimSlepSnu.addMassPlane([[x,.05*x+.95*y,y]]*2)
TChipChimSlepSnu_x005.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x005_EM_MAPS/MA5_EM_TChipChimSlepSnu_x005_mT2-150-DF.dat', 'txt')
TChipChimSlepSnu_x025 = TChipChimSlepSnu.addMassPlane([[x,.25*x+.75*y,y]]*2)
TChipChimSlepSnu_x025.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x025_EM_MAPS/MA5_EM_TChipChimSlepSnu_x025_mT2-150-DF.dat', 'txt')
TChipChimSlepSnu_x05 = TChipChimSlepSnu.addMassPlane([[x,.50*x+.50*y,y]]*2)
TChipChimSlepSnu_x05.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x05_EM_MAPS/MA5_EM_TChipChimSlepSnu_x05_mT2-150-DF.dat', 'txt')
TChipChimSlepSnu_x05.addSource('obsExclusion', "orig/TChipChimSlepSnu_Obs_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionP1', "orig/TChipChimSlepSnu_Obs+1_x05.txt", "txt")
TChipChimSlepSnu_x05.addSource('obsExclusionM1', "orig/TChipChimSlepSnu_Obs-1_x05.txt", "txt")
TChipChimSlepSnu_x075 = TChipChimSlepSnu.addMassPlane([[x,.75*x+.25*y,y]]*2)
TChipChimSlepSnu_x075.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x075_EM_MAPS/MA5_EM_TChipChimSlepSnu_x075_mT2-150-DF.dat', 'txt')
TChipChimSlepSnu_x095 = TChipChimSlepSnu.addMassPlane([[x,.95*x+.05*y,y]]*2)
TChipChimSlepSnu_x095.addSource('efficiencyMap','orig/atlas_susy_2013_11_TChipChimSlepSnu_x095_EM_MAPS/MA5_EM_TChipChimSlepSnu_x095_mT2-150-DF.dat', 'txt')
TChipChimSlepSnu_DELTASleptonNeutralino5 = TChipChimSlepSnu.addMassPlane([[x,y+5.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino5_mT2-150-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino10 = TChipChimSlepSnu.addMassPlane([[x,y+10.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino10_mT2-150-DF.dat", "txt")
TChipChimSlepSnu_DELTASleptonNeutralino15 = TChipChimSlepSnu.addMassPlane([[x,y+15.0,y]]*2)
TChipChimSlepSnu_DELTASleptonNeutralino15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTASleptonNeutralino15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTASleptonNeutralino15_mT2-150-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton5 = TChipChimSlepSnu.addMassPlane([[x,x-5.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton5.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton5_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton5_mT2-150-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton10 = TChipChimSlepSnu.addMassPlane([[x,x-10.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton10.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton10_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton10_mT2-150-DF.dat", "txt")
TChipChimSlepSnu_DELTACharginoSlepton15 = TChipChimSlepSnu.addMassPlane([[x,x-15.0,y]]*2)
TChipChimSlepSnu_DELTACharginoSlepton15.addSource('efficiencyMap',"orig/atlas_susy_2013_11_TChipChimSlepSnu_DELTACharginoSlepton15_EM_MAPS/MA5_EM_TChipChimSlepSnu_DELTACharginoSlepton15_mT2-150-DF.dat", "txt")

databaseCreator.create()
