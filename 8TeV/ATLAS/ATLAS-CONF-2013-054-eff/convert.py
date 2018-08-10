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



#+++++++ Datasets info ++++++++++++++
datasetsInfo = [
{'dataId': '8j50 flavor >=2 b-jets', 'dataType': 'efficiencyMap', 'observedN': '44', 'expectedUpperLimit': '1.24*fb', 'folder': 'data-cut2', 'expectedBG': '50.0', 'bgError': '10.0', 'upperLimit': '1.1*fb'},
{'dataId': '8j50 flavor 1 b-jets', 'dataType': 'efficiencyMap', 'observedN': '44', 'expectedUpperLimit': '1.21*fb', 'folder': 'data-cut1', 'expectedBG': '40.0', 'bgError': '10.0', 'upperLimit': '1.1*fb'},
{'dataId': '9j50 flavor >=2 b-jets', 'dataType': 'efficiencyMap', 'observedN': '7', 'expectedUpperLimit': '0.448*fb', 'folder': 'data-cut5', 'expectedBG': '8.0', 'bgError': '2.7', 'upperLimit': '0.37*fb'},
{'dataId': '8j80 flavor >=2 b-jets', 'dataType': 'efficiencyMap', 'observedN': '3', 'expectedUpperLimit': '0.298*fb', 'folder': 'data-cut12', 'expectedBG': '3.3', 'bgError': '2.2', 'upperLimit': '0.31*fb'},
{'dataId': '9j50 MJ 420', 'dataType': 'efficiencyMap', 'observedN': '9', 'expectedUpperLimit': '0.601*fb', 'folder': 'data-cut16', 'expectedBG': '11.0', 'bgError': '5.0', 'upperLimit': '0.5*fb'},
{'dataId': '8j80 flavor 0 b-jets', 'dataType': 'efficiencyMap', 'observedN': '2', 'expectedUpperLimit': '0.147*fb', 'folder': 'data-cut10', 'expectedBG': '0.9', 'bgError': '0.6', 'upperLimit': '0.24*fb'},
{'dataId': '8j50 MJ 420', 'dataType': 'efficiencyMap', 'observedN': '37', 'expectedUpperLimit': '1.54*fb', 'folder': 'data-cut14', 'expectedBG': '45.0', 'bgError': '14.0', 'upperLimit': '1*fb'},
{'dataId': '7j80 flavor >=2 b-jets', 'dataType': 'efficiencyMap', 'observedN': '13', 'expectedUpperLimit': '1.08*fb', 'folder': 'data-cut9', 'expectedBG': '25.0', 'bgError': '10.0', 'upperLimit': '0.6*fb'},
{'dataId': '10j50 MJ 420', 'dataType': 'efficiencyMap', 'observedN': '1', 'expectedUpperLimit': '0.254*fb', 'folder': 'data-cut18', 'expectedBG': '2.2', 'bgError': '2.0', 'upperLimit': '0.2*fb'},
{'dataId': '8j80 flavor 1 b-jets', 'dataType': 'efficiencyMap', 'observedN': '1', 'expectedUpperLimit': '0.2*fb', 'folder': 'data-cut11', 'expectedBG': '1.5', 'bgError': '0.9', 'upperLimit': '0.17*fb'},
{'dataId': '10j50 flavor', 'dataType': 'efficiencyMap', 'observedN': '3', 'expectedUpperLimit': '0.2*fb', 'folder': 'data-cut6', 'expectedBG': '1.4', 'bgError': '0.3', 'upperLimit': '0.29*fb'},
{'dataId': '7j80 flavor 1 b-jets', 'dataType': 'efficiencyMap', 'observedN': '17', 'expectedUpperLimit': '0.743*fb', 'folder': 'data-cut8', 'expectedBG': '17.0', 'bgError': '6.0', 'upperLimit': '0.8*fb'},
{'dataId': '8j50 MJ 340', 'dataType': 'efficiencyMap', 'observedN': '69', 'expectedUpperLimit': '2.08*fb', 'folder': 'data-cut13', 'expectedBG': '75.0', 'bgError': '19.0', 'upperLimit': '1.7*fb'},
{'dataId': '10j50 MJ 340', 'dataType': 'efficiencyMap', 'observedN': '1', 'expectedUpperLimit': '0.302*fb', 'folder': 'data-cut17', 'expectedBG': '3.2', 'bgError': '3.5', 'upperLimit': '0.2*fb'},
{'dataId': '8j50 flavor 0 b-jets', 'dataType': 'efficiencyMap', 'observedN': '40', 'expectedUpperLimit': '0.762*fb', 'folder': 'data-cut0', 'expectedBG': '35.0', 'bgError': '4.0', 'upperLimit': '0.97*fb'},
{'dataId': '9j50 flavor 0 b-jets', 'dataType': 'efficiencyMap', 'observedN': '5', 'expectedUpperLimit': '0.265*fb', 'folder': 'data-cut3', 'expectedBG': '3.3', 'bgError': '0.7', 'upperLimit': '0.34*fb'},
{'dataId': '9j50 flavor 1 b-jets', 'dataType': 'efficiencyMap', 'observedN': '8', 'expectedUpperLimit': '0.369*fb', 'folder': 'data-cut4', 'expectedBG': '6.1', 'bgError': '1.7', 'upperLimit': '0.43*fb'},
{'dataId': '9j50 MJ 340', 'dataType': 'efficiencyMap', 'observedN': '13', 'expectedUpperLimit': '0.797*fb', 'folder': 'data-cut15', 'expectedBG': '17.0', 'bgError': '7.0', 'upperLimit': '0.5*fb'},
{'dataId': '7j80 flavor 0 b-jets', 'dataType': 'efficiencyMap', 'observedN': '12', 'expectedUpperLimit': '0.463*fb', 'folder': 'data-cut7', 'expectedBG': '11.0', 'bgError': '2.2', 'upperLimit': '0.5*fb'}
]

#+++++++ Constraints info ++++++++++++++
constraintsDict = {
"T2gg" : "[[['jet']],[['jet']]]",
"T1" : "[[['jet','jet']],[['jet','jet']]]",
"T1bbbb" : "[[['b','b']],[['b','b']]]",
"T5btbt" : "[[['b'],['t']],[['b'],['t']]]",
"T1bbqq" : "[[['b','b']],[['jet','jet']]]",
"T5bbbt" : "[[['b'],['b']],[['b'],['t']]]",
"T1bbbt" : "[[['b','b']],[['b','t']]]",
"TGQbtq" : "[[['b','t']],[['jet']]]",
"TGQ" : "[[['jet']],[['jet','jet']]]",
"TGQqtt" : "[[['jet']],[['t+','t-']]]",
"T2tt" : "[[['t']],[['t']]]",
"T5tttt" : "[[['t+'],['t-']],[['t+'],['t-']]]+[[['t-'],['t+']],[['t-'],['t+']]]",
"T1btqq" : "[[['b','t']],[['jet','jet']]]",
"T1btbt" : "[[['b','t']],[['b','t']]]",
"T1tttt" : "[[['t+','t-']],[['t+','t-']]]",
"T2" : "[[['jet']],[['jet']]]",
"T5tbtb" : "[[['t'],['b']],[['t'],['b']]]",
"T5tbtt" : "[[['t'],['b']],[['t+'],['t-']]]+[[['t'],['b']],[['t-'],['t+']]]",
"TGQbbq" : "[[['b','b']],[['jet']]]",
"T2bt" : "[[['b']],[['t']]]",
"T1bbtt" : "[[['b','b']],[['t+','t-']]]",
"T5bbbb" : "[[['b'],['b']],[['b'],['b']]]",
"T1qqtt" : "[[['jet','jet']],[['t+','t-']]]",
"T1bttt" : "[[['b','t']],[['t+','t-']]]",
"T2bb" : "[[['b']],[['b']]]"
}

#+++++++ Figure info ++++++++++++++
figure={}
figure["T1"]='fig_07a.pdf'
figureUrl={}
figureUrl["T1"]='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_07a.pdf'


#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-CONF-2013-054')
info.comment = 'created from fastlim-1.0;  superseded by ATLAS-SUSY-2013-04, fastlim maps in ATLAS-CONF-2013-054 contain different topologies than ATLAS-SUSY-2013-04(MA5 maps). Topologies included here are, T1, T1bbbb, T1bbbt, T1bbqq, T1bbtt, T1btbt, T1btqq, T1bttt, T1qqtt, T1tttt, T2, T2bb, T2bt, T2tt, T5bbbb, T5bbbt, T5btbt, T5tbtb, T5tbtt, T5tttt, TGQ, TGQbbq, TGQbtq, TGQQtt.'
info.lastUpdate = '2016/08/16'
info.implementedBy = 'WW'
info.sqrts = '8*TeV'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-054/'
info.contact = 'fastlim'
info.lumi = '20.3/fb'
info.prettyName = '0 leptons + >= 7-10 jets + Etmiss'


for dt in datasetsInfo:
    folder = dt.pop('folder')
    dataset = DataSetInput(folder)        
    dataset.setInfo(**dt)
    origDir = os.path.basename(folder)+"-orig/"
    for i in os.listdir(origDir):
        if i[-5:]!=".effi": continue
        txname=i[:-5]
        tmp = dataset.addTxName(txname)
        tmp.constraint = constraintsDict[txname]
        tmp.conditionDescription=None
        tmp.condition = None
        tmp.source = 'Fastlim-v1.0'
        tmp.dataUrl = None
        if i[:2] in [ "T5", "T6" ]:
            tmp_1 = tmp.addMassPlane([[x,y,z]]*2)
            tmp_1.addSource('efficiencyMap', origDir+'%s.effi' % txname, 'effi')
        else:
            tmp_1 = tmp.addMassPlane([[x,y]]*2)
            tmp_1.addSource('efficiencyMap', origDir+'%s.effi' % txname, 'effi')
        if os.path.exists(origDir+'%s_excl.dat' % txname):
            tmp_1.addSource('obsExclusion', origDir+'%s_excl.dat' % txname, 'txt')
        if txname in figure:
            tmp_1.figure = figure[txname]
        if txname in figureUrl:
            tmp_1.figureUrl = figureUrl[txname]
databaseCreator.create()

