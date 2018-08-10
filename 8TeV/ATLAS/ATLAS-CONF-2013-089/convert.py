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
info = MetaInfoInput('ATLAS-CONF-2013-089')
info.comment = 'There are more topologies (with more than three masses). Superseded by ATLAS-SUSY-2013-20.'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-089/; https://cds.cern.ch/record/1595272'
info.prettyName = '2 leptons (e,mu) + Etmiss'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T6WW = dataset.addTxName('T6WW')
T6WW.checked ="VM"
T6WW.constraint ="[[['jet'],['W']],[['jet'],['W']]]"
T6WW.conditionDescription ="None"
T6WW.condition ="None"
T6WW.source = "ATLAS"
T6WW.massConstraint = None
T6WWoff = dataset.addTxName('T6WWoff')
T6WWoff.constraint = "22*[[['jet'],['l','nu']],[['jet'],['l','nu']]]"
T6WWoff.conditionDescription = "[[['jet'],['mu','nu']],[['jet'],['mu','nu']]] > [[['jet'], ['e','nu']],[['jet'],['e','nu']]],[[['jet'],['e','nu']],[['jet'],['mu','nu']]] > 2*[[['jet'], ['e','nu']],[['jet'],['e','nu']]]"
T6WWoff.condition = "Cgtr([[['jet'],['mu','nu']],[['jet'],['mu','nu']]], [[['jet'],['e','nu']],[['jet'], ['e','nu']]]);Cgtr([[['jet'],['e','nu']],[['jet'],['mu','nu']]], 2.*[[['jet'],['e','nu']],[['jet'], ['e','nu']]])"
T6WWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6WWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6WW050 = T6WW.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
T6WW050.figure = 'Fig. 23a'
T6WW050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-089/fig_23a.png'
T6WW050.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-089/fig_23a_PRELIMINARY.data'
T6WW050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6WW050.txt', 'orig/T6WW050.data'],
                 dataFormats= ['txt', 'txt'])
T6WWoff.addMassPlane(T6WW050)
#+++++++ next mass plane block ++++++++++++++
T6WWLSP060 = T6WW.addMassPlane(2*[[x, y*(x-60.0)+60.0, 60.0]])
T6WWLSP060.figure = 'Fig. 23b'
T6WWLSP060.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-089/fig_23b.png'
T6WWLSP060.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-089/fig_23b_PRELIMINARY.data'
T6WWLSP060.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6WWLSP060.txt', 'orig/T6WWLSP060.data'],
                 dataFormats= ['txt', 'txt'])
T6WWoff.addMassPlane(T6WWLSP060)

#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.checked ="VM"
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription ="None"
T5WW.condition ="None"
T5WW.source = "ATLAS"
T5WW.massConstraint = None
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="20.2*[[['jet','jet'],['l','nu']],[['jet','jet'],['l','nu']]]"
T5WWoff.conditionDescription = "[[['jet','jet'],['mu','nu']],[['jet','jet'],['mu','nu']]] > [[['jet','jet'],['e','nu']],[['jet','jet'],['e','nu']]], [[['jet','jet'],['e','nu']],[['jet','jet'],['mu','nu']]] > 2.* [[['jet','jet'],['e','nu']],[['jet','jet'],['e','nu']]]"
T5WWoff.condition = "Cgtr([[['jet','jet'],['mu','nu']],[['jet','jet'],['mu','nu']]],[[['jet','jet'], ['e','nu']],[['jet','jet'],['e','nu']]]);Cgtr([[['jet','jet'],['e','nu']],[['jet','jet'],['mu','nu']]],2.*[[['jet','jet'], ['e','nu']],[['jet','jet'],['e','nu']]])"
T5WWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T5WWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T5WWLSP060 = T5WW.addMassPlane(2*[[x, y*x-y*60.0+60.0, 60.0]])
T5WWLSP060.figure = 'Fig (aux). 22b'
T5WWLSP060.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-089/fig_22b.png'
T5WWLSP060.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-089/fig_22b_PRELIMINARY.data'
T5WWLSP060.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T5WWLSP060.txt', 'orig/T5WWLSP060.data'],
                 dataFormats= ['txt', 'txt'])
T5WWoff.addMassPlane(T5WWLSP060)
#+++++++ next mass plane block ++++++++++++++
T5WW050 = T5WW.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
T5WW050.figure = 'Fig. 22a'
T5WW050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-089/fig_22a.png'
T5WW050.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-089/fig_22a_PRELIMINARY.data'
T5WW050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T5WW050.txt', 'orig/T5WW050.data'],
                 dataFormats= ['txt', 'txt'])
T5WWoff.addMassPlane(T5WW050)



databaseCreator.create()
