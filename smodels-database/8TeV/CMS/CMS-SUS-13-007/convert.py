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
info.supersedes =''
info.implementedBy = 'Federico Ambrogi'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure = 'combLimit_T1tttt_b.pdf'
T1tttt_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/combLimit_T1tttt_b.pdf'
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/limits_model_A.txt'
T1tttt_1.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['./orig/T1tttt_excl.txt', './orig/T1tttt.txt'],
                 dataFormats= ['txt', 'txt'],objectNames= [None, 'None'],units= [None, 'fb'])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++++ next txName block ++++++++++++++
T5tttt = dataset.addTxName('T5tttt')
T5tttt.constraint = "[[['t'],['t']],[['t'],['t']]]"
T5tttt.conditionDescription = None
T5tttt.condition = None
T5tttt.massConstraint = None
T5tttt.source = 'CMS'
T5ttttoff = dataset.addTxName('T5ttttoff')
T5ttttoff.constraint = "[[['b','W'],['b','W']],[['b','W'],['b','W']]]"
T5ttttoff.conditionDescription = None
T5ttttoff.condition = None
T5ttttoff.massConstraint = [['dm <= 169.0', 'dm <= 169.0'], ['dm <= 169.0', 'dm <= 169.0']]
T5ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T5tttt_1 = T5tttt.addMassPlane(2*[[1000., x, y]])
T5tttt_1.figure = 'combLimit_T1t1t_b.pdf'
T5tttt_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/combLimit_T1t1t_b.png'
T5tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/limits_model_B.txt'
T5tttt_1.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['./orig/T5tttt_mg1TeV_excl.txt', './orig/T5tttt_mg1TeV.txt'],
                 dataFormats= ['txt', 'txt'],objectNames= [None, 'None'],units= [None, 'fb'])
T5ttttoff.addMassPlane(T5tttt_1)
#+++++++ next mass plane block ++++++++++++++
T5tttt_2 = T5tttt.addMassPlane(2*[[x, y, 50.]])
T5tttt_2.figure = 'combLimit_T5tttt_b.pdf'
T5tttt_2.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/combLimit_T5tttt_b.png'
T5tttt_2.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/limits_model_C.txt'
T5tttt_2.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['./orig/T5tttt_mLSP50GeV_excl.txt', './orig/T5tttt_mLSP50GeV.txt'],
                 dataFormats= ['txt', 'txt'],objectNames= [None, 'None'],units= [None, 'fb'])
T5ttttoff.addMassPlane(T5tttt_2)



databaseCreator.create()
