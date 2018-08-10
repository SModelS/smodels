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
{'dataId': 'SL2b_low', 'dataType': 'efficiencyMap', 'observedN': '24', 'expectedUpperLimit': '0.688*fb', 'folder': 'data-cut2', 'expectedBG': '24.1', 'bgError': '4.1', 'upperLimit': '0.66*fb'},
{'dataId': 'SL1b_high', 'dataType': 'efficiencyMap', 'observedN': '6', 'expectedUpperLimit': '0.309*fb', 'folder': 'data-cut1', 'expectedBG': '4.0', 'bgError': '1.1', 'upperLimit': '0.39*fb'},
{'dataId': 'SL5j', 'dataType': 'efficiencyMap', 'observedN': '9', 'expectedUpperLimit': '0.552*fb', 'folder': 'data-cut5', 'expectedBG': '14.8', 'bgError': '3.7', 'upperLimit': '0.35*fb'},
{'dataId': 'incHL6j_m', 'dataType': 'efficiencyMap', 'observedN': '0', 'expectedUpperLimit': '0.195*fb', 'folder': 'data-cut12', 'expectedBG': '1.7', 'bgError': '0.6', 'upperLimit': '0.15*fb'},
{'dataId': 'incHL5j_m', 'dataType': 'efficiencyMap', 'observedN': '2', 'expectedUpperLimit': '0.237*fb', 'folder': 'data-cut10', 'expectedBG': '2.5', 'bgError': '0.8', 'upperLimit': '0.22*fb'},
{'dataId': 'incHL5j_e', 'dataType': 'efficiencyMap', 'observedN': '4', 'expectedUpperLimit': '0.269*fb', 'folder': 'data-cut9', 'expectedBG': '3.6', 'bgError': '1.0', 'upperLimit': '0.3*fb'},
{'dataId': 'incHL6j_e', 'dataType': 'efficiencyMap', 'observedN': '2', 'expectedUpperLimit': '0.241*fb', 'folder': 'data-cut11', 'expectedBG': '2.0', 'bgError': '0.7', 'upperLimit': '0.23*fb'},
{'dataId': 'SL2m', 'dataType': 'efficiencyMap', 'observedN': '7', 'expectedUpperLimit': '0.204*fb', 'folder': 'data-cut6', 'expectedBG': '1.6', 'bgError': '1.0', 'upperLimit': '0.57*fb'},
{'dataId': 'incHL3j_m', 'dataType': 'efficiencyMap', 'observedN': '5', 'expectedUpperLimit': '0.235*fb', 'folder': 'data-cut8', 'expectedBG': '2.7', 'bgError': '0.9', 'upperLimit': '0.38*fb'},
{'dataId': 'SL1b_low', 'dataType': 'efficiencyMap', 'observedN': '8', 'expectedUpperLimit': '0.353*fb', 'folder': 'data-cut0', 'expectedBG': '6.1', 'bgError': '1.4', 'upperLimit': '0.43*fb'},
{'dataId': 'SL2b_high', 'dataType': 'efficiencyMap', 'observedN': '3', 'expectedUpperLimit': '0.276*fb', 'folder': 'data-cut3', 'expectedBG': '3.6', 'bgError': '1.4', 'upperLimit': '0.26*fb'},
{'dataId': 'SL3j', 'dataType': 'efficiencyMap', 'observedN': '7', 'expectedUpperLimit': '0.335*fb', 'folder': 'data-cut4', 'expectedBG': '5.6', 'bgError': '1.6', 'upperLimit': '0.4*fb'},
{'dataId': 'incHL3j_e', 'dataType': 'efficiencyMap', 'observedN': '4', 'expectedUpperLimit': '0.26*fb', 'folder': 'data-cut7', 'expectedBG': '3.9', 'bgError': '1.0', 'upperLimit': '0.3*fb'}
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
info = MetaInfoInput('ATLAS-CONF-2013-062')
info.comment = 'created from fastlim-1.0'
info.lastUpdate = '2016/08/16'
info.implementedBy = 'WW'
info.sqrts = '8*TeV'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-062/'
info.contact = 'fastlim'
info.lumi = '20.3/fb'
info.prettyName = '1 lepton + jets + ETmiss'


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

