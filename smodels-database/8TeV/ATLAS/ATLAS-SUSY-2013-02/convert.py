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
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication = 'http://link.springer.com/article/10.1007/JHEP09%282014%29176'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/'
info.arxiv = 'http://arxiv.org/abs/1405.7875'
info.prettyName = '0 leptons + 2-6 jets + Etmiss'
info.supersedes = 'ATLAS-CONF-2013-047'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TGQ = dataset.addTxName('TGQ')
TGQ.constraint ="[[['jet']],[['jet','jet']]]"
TGQ.conditionDescription ="None"
TGQ.condition ="None"
TGQ.source = 'ATLAS'
TGQ.round_to = 6 ## round to 6 digits to make PCA work.
#+++++++ next mass plane block ++++++++++++++
TGQ0 = TGQ.addMassPlane([[0.96*x, y],[x,y]])
TGQ0.figure = "fig(aux) 9b"
TGQ0.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_09b.png"
TGQ0.dataUrl = 'Not defined'
TGQ0.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TGQ.txt', 'orig/limit_TGQ.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
T6WW = dataset.addTxName('T6WW')
T6WW.constraint ="[[['jet'],['W']],[['jet'],['W']]]"
T6WW.conditionDescription ="None"
T6WW.condition = "None"
T6WW.massConstraint = None
T6WW.source = 'ATLAS'
T6WWoff = dataset.addTxName('T6WWoff')
T6WWoff.constraint = "2.23 * [[['jet'],['jet','jet']],[['jet'],['jet','jet']]]"
T6WWoff.conditionDescription = "None"
T6WWoff.condition ="None"
T6WWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6WWoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T6WWLSP060 = T6WW.addMassPlane(2*[[x, y*(x-60.0)+60.0, 60.0]])
T6WWLSP060.figure = "fig(aux) 11c"                                                                                              
T6WWLSP060.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_11c.png"                      
T6WWLSP060.dataUrl = 'Not defined'
T6WWLSP060.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6WWLSP060.txt', 'orig/limit_T6WWLSP060.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])
T6WWoff.addMassPlane(T6WWLSP060)
#+++++++ next mass plane block ++++++++++++++
T6WW050 = T6WW.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
T6WW050.figure = "fig(aux) 11d"
T6WW050.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_11d.png"
T6WW050.dataUrl = 'Not defined'
T6WW050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6WW050.txt', 'orig/limit_T6WW050.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])
T6WWoff.addMassPlane(T6WW050)

#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane(2*[[x, y]])
T2.figure = "Fig(aux) 8b"
T2.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_09c.png"
T2.dataUrl = 'Not defined'
T2.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T2.txt', 'orig/limit_T2.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane(2*[[x, y]])
T1.figure = "Fig(aux) 9a"
T1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_09a.png"
T1.dataUrl = 'Not defined'
T1.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T1.txt', 'orig/limit_T1.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription ="None"
T5WW.condition ="None"
T5WW.massConstraint = None
T5WW.source = 'ATLAS'
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint = "2.23 * [[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription = "None"
T5WWoff.condition = "None"
T5WWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T5WWoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T5WWLSP060 = T5WW.addMassPlane(2*[[x, y*(x-60.)+60.0, 60.0]])
T5WWLSP060.figure = "fig(aux) 10a"
T5WWLSP060.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_11a.png"
T5WWLSP060.dataUrl = 'Not defined'
T5WWLSP060.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T5WWLSP060.txt', 'orig/limit_T5WWLSP060.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])
T5WWoff.addMassPlane(T5WWLSP060)
#+++++++ next mass plane block ++++++++++++++
T5WW050 = T5WW.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
T5WW050.figure = "fig(aux) 11b"
T5WW050.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_11b.png"
T5WW050.dataUrl = 'Not defined'
T5WW050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T5WW050.txt', 'orig/limit_T5WW050.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])
T5WWoff.addMassPlane(T5WW050)

#+++++++ next txName block ++++++++++++++
T5tctc = dataset.addTxName('T5tctc')
T5tctc.constraint ="[[['t'],['jet']],[['t'],['jet']]]"
T5tctc.conditionDescription ="None"
T5tctc.condition ="None"
T5tctc.massConstraint = None
T5tctc.source = 'ATLAS'
T5tctcoff = dataset.addTxName('T5tctcoff')
T5tctcoff.constraint ="[[['W','b'],['jet']],[['W','b'],['jet']]]"
T5tctcoff.massConstraint = [['dm <= 169.0', 'dm >= 0.0'], ['dm <= 169.0', 'dm >= 0.0']]
T5tctcoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T5tctcD020 = T5tctc.addMassPlane(2*[[x, y, y-20.0]])
T5tctcD020.figure = "fig(aux) 13"
T5tctcD020.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_13.png"
T5tctcD020.dataUrl = 'Not defined'
T5tctcD020.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T5tctcD020.txt', 'orig/limit_T5tctcD020.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])
T5tctcoff.addMassPlane(T5tctcD020)



databaseCreator.create()
