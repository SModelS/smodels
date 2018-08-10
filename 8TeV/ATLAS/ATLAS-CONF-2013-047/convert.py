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
info = MetaInfoInput('ATLAS-CONF-2013-047')
info.comment = 'T5WWLSP060 and T6WWLSP060 originally have xvalue on y-axes, changed by us to M2'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/; http://cds.cern.ch/record/1547563'
info.supersededBy = 'ATLAS-SUSY-2013-02'
info.prettyName = '0 leptons + 2-6 jets + Etmiss'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TGQ = dataset.addTxName('TGQ')
TGQ.checked ="VM"
TGQ.constraint ="[[['jet']],[['jet','jet']]]"
TGQ.conditionDescription ="None"
TGQ.condition ="None"
TGQ.source = 'ATLAS'
TGQ.round_to = 6 ## round to 6 digits to make PCA work.
#+++++++ next mass plane block ++++++++++++++
TGQ0 = TGQ.addMassPlane([[0.96*x, y], [x,y]])
TGQ0.figure = 'Fig. 19b'
TGQ0.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_19b.png'
TGQ0.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_19b_PRELIMINARY.data'
TGQ0.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/TGQ_exc.dat', 'orig/TGQ.dat'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
T6WW = dataset.addTxName('T6WW')
T6WW.checked ="VM"
T6WW.constraint ="[[['jet'],['W']],[['jet'],['W']]]"
T6WW.conditionDescription ="None"
T6WW.condition ="None"
T6WW.massConstraint = None
T6WW.source = 'ATLAS'
T6WWoff = dataset.addTxName('T6WWoff')
T6WWoff.constraint ="2.23 * [[['jet'],['jet','jet']],[['jet'],['jet','jet']]]"
T6WWoff.conditionDescription = "None"
T6WWoff.condition = "None"
T6WWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6WWoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T6WWLSP060 = T6WW.addMassPlane(2*[[x, y, 60.0]])
T6WWLSP060.figure = 'Fig. 21c'
T6WWLSP060.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_21c.png'
T6WWLSP060.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_21c_PRELIMINARY.data'
T6WWLSP060.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T6WWLSP060_exc.dat', 'orig/T6WWLSP060.dat'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])
T6WWoff.addMassPlane(T6WWLSP060)
#+++++++ next mass plane block ++++++++++++++
T6WW050 = T6WW.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
T6WW050.figure = 'Fig. 21d'
T6WW050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_21d.png'
T6WW050.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_21d_PRELIMINARY.data'
T6WW050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T6WW050_exc.dat', 'orig/T6WW050.dat'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])
T6WWoff.addMassPlane(T6WW050)

#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.checked ="VM"
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane(2*[[x, y]])
T2.figure = 'Fig. 19c'
T2.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_19c.png'
T2.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_19c_PRELIMINARY.data'
T2.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2_exc.dat', 'orig/T2.dat'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.checked ="VM"
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane(2*[[x, y]])
T1.figure = 'Fig. 19a'
T1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_19a.png'
T1.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_19a_PRELIMINARY.data'
T1.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T1_exc.dat', 'orig/T1.dat'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.checked ="VM"
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
T5WWLSP060 = T5WW.addMassPlane(2*[[x, y, 60.0]])
T5WWLSP060.figure = 'Fig. 21a'
T5WWLSP060.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_21a.png'
T5WWLSP060.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_21b_PRELIMINARY.data'
T5WWLSP060.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T5WWLSP60_exc.dat', 'orig/T5WWLSP060.dat'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])
T5WWoff.addMassPlane(T5WWLSP060)
#+++++++ next mass plane block ++++++++++++++
T5WW050 = T5WW.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
T5WW050.figure = 'Fig. 21b'
T5WW050.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_21b.png'
T5WW050.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_21a_PRELIMINARY.data'
T5WW050.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T5WW050_exc.dat', 'orig/T5WW050.dat'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])
T5WWoff.addMassPlane(T5WW050)

#+++++++ next txName block ++++++++++++++
T5tctc = dataset.addTxName('T5tctc')
T5tctc.checked ="VM"
T5tctc.constraint ="[[['t'],['jet']],[['t'],['jet']]]"
T5tctc.conditionDescription ="None"
T5tctc.condition ="None"
T5tctc.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T5tctcD020 = T5tctc.addMassPlane(2*[[x, y, y-20.0]])
T5tctcD020.figure = "fig 24"
T5tctcD020.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_24.png"
T5tctcD020.dataUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_24_PRELIMINARY.data"
T5tctcD020.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T5tctc_exc.dat', 'orig/T5tctc.dat'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])



databaseCreator.create()
