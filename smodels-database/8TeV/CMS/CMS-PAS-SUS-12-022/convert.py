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
info = MetaInfoInput('CMS-PAS-SUS-12-022')
info.sqrts = '8.0'
info.private = False
info.lumi = '9.2'
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS12022'
info.supersededBy = 'CMS-SUS-13-006 (more data)'
info.contact = 'Benjamin Hooberman <hooberman@gmail.com>, Lesya Shchutska <lesya.shchutska@cern.ch>'
info.prettyName = 'multi-lepton + Etmiss final states (exactly 3 or 4 leptons, 2 SS leptons, 2 SFOS leptons + 2 jets, etc)'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiWZ = dataset.addTxName('TChiWZ')
TChiWZ.checked ="A"
TChiWZ.constraint ="[[['W']],[['Z']]]"
TChiWZ.conditionDescription ="None"
TChiWZ.condition ="None"
TChiWZ.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiWZ = TChiWZ.addMassPlane(2*[[x, y]])
TChiWZ.dataUrl = 'Not defined'
TChiWZ.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiWZ.root', 'orig/exclusion_TChiWZ.root'],
                 dataFormats= ['canvas', 'canvas'],objectNames= ['interpret', 'interpret'],indices= [8, 2],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
TChiChipmStauStau = dataset.addTxName('TChiChipmStauStau')
TChiChipmStauStau.checked ="A"
TChiChipmStauStau.constraint ="[[['ta'],['ta']],[['nu'],['ta']]]"
TChiChipmStauStau.conditionDescription ="None"
TChiChipmStauStau.condition ="None"
TChiChipmStauStau.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiChipmStauStau050 = TChiChipmStauStau.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChiChipmStauStau050.figure = "Fig. 17"
TChiChipmStauStau050.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12022/exclusion_TChiStauSnu_0_5.pdf"
TChiChipmStauStau050.dataUrl = 'Not defined'
TChiChipmStauStau050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiStauSnu_0_5.root', 'orig/exclusion_TChiStauSnu_0_5.root'],
                 dataFormats= ['canvas', 'canvas'],objectNames= ['interpret', 'interpret'],indices= [8, 2],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.checked ="A"
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TSlepSlep = TSlepSlep.addMassPlane(2*[[x, y]])
TSlepSlep.figure = "Fig. 21"
TSlepSlep.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12022/exclusion_TSlepSlep.pdf"
TSlepSlep.dataUrl = 'Not defined'
TSlepSlep.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TSlepSlep.root', 'orig/exclusion_TSlepSlep.root'],
                 dataFormats= ['canvas', 'canvas'],objectNames= ['interpret', 'interpret'],indices= [8, 2],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
TChiChipmSlepStau = dataset.addTxName('TChiChipmSlepStau')
TChiChipmSlepStau.checked ="A"
TChiChipmSlepStau.constraint ="[[['L'],['L']],[['nu'],['ta']]]"
TChiChipmSlepStau.conditionDescription ="[[['L'],['L']],[['nu'],['ta']]] > 2.7*[[['ta'],['ta']],[['nu'],['ta']]],[[['L'],['L']],[['nu'],['ta']]] > 2.7*[[['e'],['e']],[['nu'],['ta']]]"
TChiChipmSlepStau.condition ="Cgtr([[['L'],['L']],[['nu'],['ta']]],3.*[[['ta'],['ta']],[['nu'],['ta']]]);Cgtr([[['L'],['L']],[['nu'],['ta']]],3.*[[['e'],['e']],[['nu'],['ta']]])"
TChiChipmSlepStau.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepStau095 = TChiChipmSlepStau.addMassPlane(2*[[x, x*0.95+(1.-0.95)*y, y]])
TChiChipmSlepStau095.figure = "Fig. 16a"
TChiChipmSlepStau095.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12022/exclusion_TChiSlepSnu_2a_0_95.pdf"
TChiChipmSlepStau095.dataUrl = 'Not defined'
TChiChipmSlepStau095.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2a_0_95.root', 'orig/exclusion_TChiSlepSnu_2a_0_95.root'],
                 dataFormats= ['canvas', 'canvas'],objectNames= ['interpret', 'interpret'],indices= [8, 2],units= [None, 'fb'])
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepStau050 = TChiChipmSlepStau.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChiChipmSlepStau050.figure = "Fig. 16a"
TChiChipmSlepStau050.figureUrl =  "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12022/exclusion_TChiSlepSnu_2a_0_5.pdf"
TChiChipmSlepStau050.dataUrl = 'Not defined'
TChiChipmSlepStau050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2a_0_5.root', 'orig/exclusion_TChiSlepSnu_2a_0_5.root'],
                 dataFormats= ['canvas', 'canvas'],objectNames= ['interpret', 'interpret'],indices= [8, 2],units= [None, 'fb'])
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepStau005 = TChiChipmSlepStau.addMassPlane(2*[[x, x*0.05+(1.-0.05)*y, y]])
TChiChipmSlepStau005.figure = "Fig. 16a"
TChiChipmSlepStau005.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12022/exclusion_TChiSlepSnu_2a_0_05.pdf"
TChiChipmSlepStau005.dataUrl = 'Not defined'
TChiChipmSlepStau005.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2a_0_05.root', 'orig/exclusion_TChiSlepSnu_2a_0_05.root'],
                 dataFormats= ['canvas', 'canvas'],objectNames= ['interpret', 'interpret'],indices= [8, 2],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
TChiChipmSlepL = dataset.addTxName('TChiChipmSlepL')
TChiChipmSlepL.checked ="A"
TChiChipmSlepL.constraint ="2.*([[['L'],['L']],[['L'],['nu']]] + [[['L'],['L']],[['nu'],['L']]])"
TChiChipmSlepL.conditionDescription ="[[['L'],['L']],[['L'],['nu']]] ~ [[['L'],['L']],[['nu'],['L']]],[[['L'],['L']],[['nu'],['L']]] > 2.7*[[['ta'],['ta']],[['nu'],['L']]],[[['L'],['L']],[['L'],['nu']]] > 2.7*[[['ta'],['ta']],[['L'],['nu']]],[[['L'],['L']],[['nu'],['L']]] > 2.7*[[['L'],['L']],[['nu'],['ta']]],[[['L'],['L']],[['L'],['nu']]] > 2.7*[[['L'],['L']],[['ta'],['nu']]],[[['L'],['L']],[['nu'],['L']]] > 2.7*[[['e'],['e']],[['nu'],['L']]],[[['L'],['L']],[['L'],['nu']]] > 2.7*[[['e'],['e']],[['L'],['nu']]],[[['L'],['L']],[['nu'],['L']]] > 2.7*[[['L'],['L']],[['nu'],['e']]],[[['L'],['L']],[['L'],['nu']]] > 2.7*[[['L'],['L']],[['e'],['nu']]]"
TChiChipmSlepL.condition ="Csim([[['L'],['L']],[['L'],['nu']]],[[['L'],['L']],[['nu'],['L']]]);Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['ta'],['ta']],[['nu'],['L']]]);Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['ta'],['ta']],[['L'],['nu']]]);Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['L'],['L']],[['nu'],['ta']]]);Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['L'],['L']],[['ta'],['nu']]]);Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['e'],['e']],[['nu'],['L']]]);Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['e'],['e']],[['L'],['nu']]]);Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['L'],['L']],[['nu'],['e']]]);Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['L'],['L']],[['e'],['nu']]])"
TChiChipmSlepL.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepL005 = TChiChipmSlepL.addMassPlane(2*[[x, x*0.05+(1.-0.05)*y, y]])
TChiChipmSlepL005.figure = "Fig. 15a" 
TChiChipmSlepL005.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12022/exclusion_TChiSlepSnu_2i_0_05.pdf"
TChiChipmSlepL005.dataUrl = 'Not defined'
TChiChipmSlepL005.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2i_0_05.root', 'orig/exclusion_TChiSlepSnu_2i_0_05.root'],
                 dataFormats= ['canvas', 'canvas'],objectNames= ['interpret', 'interpret'],indices= [8, 2],units= [None, 'fb'])
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepL095 = TChiChipmSlepL.addMassPlane(2*[[x, x*0.95+(1.-0.95)*y, y]])
TChiChipmSlepL095.figure = "Fig. 15b" 
TChiChipmSlepL095.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12022/exclusion_TChiSlepSnu_2i_0_95.pdf"
TChiChipmSlepL095.dataUrl = 'Not defined'
TChiChipmSlepL095.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2i_0_95.root', 'orig/exclusion_TChiSlepSnu_2i_0_95.root'],
                 dataFormats= ['canvas', 'canvas'],objectNames= ['interpret', 'interpret'],indices= [8, 2],units= [None, 'fb'])
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepL050 = TChiChipmSlepL.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChiChipmSlepL050.figure = "Fig. 14"
TChiChipmSlepL050.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12022/exclusion_TChiSlepSnu_2i_0_5.pdf"
TChiChipmSlepL050.dataUrl = 'Not defined'
TChiChipmSlepL050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2i_0_5.root', 'orig/exclusion_TChiSlepSnu_2i_0_5.root'],
                 dataFormats= ['canvas', 'canvas'],objectNames= ['interpret', 'interpret'],indices= [8, 2],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.checked ="A"
TChipChimSlepSnu.constraint ="[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription ="[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]],[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition ="Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]);Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]);Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]);Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]);Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]);Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]);Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]);Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]);Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]);Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]);Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]);Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]);Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]);Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]);Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]);Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]);Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"
TChipChimSlepSnu.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu050 = TChipChimSlepSnu.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
TChipChimSlepSnu050.figure = "Fig. 20"
TChipChimSlepSnu050.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS12022/exclusion_TChipmSlepSnu.pdf"
TChipChimSlepSnu050.dataUrl = 'Not defined'
TChipChimSlepSnu050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChipmSlepSnu.root', 'orig/exclusion_TChipmSlepSnu.root'],
                 dataFormats= ['canvas', 'canvas'],objectNames= ['interpret', 'interpret'],indices= [8, 2],units= [None, 'fb'])



databaseCreator.create()
