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
info = MetaInfoInput('CMS-SUS-13-011')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13011'
info.sqrts = 8
info.lumi = 19.5
info.prettyName = '1 lepton + >= 4 (1b-)jets + Etmiss'
info.private = False
info.arxiv = 'arXiv:1308.1586v2'
info.contact = 'Benjamin Hooberman <hooberman@gmail.com>, Mariarosaria DAlfonso <dalfonso@mail.cern.ch>'
info.publication = 'http://link.springer.com/article/10.1140%2Fepjc%2Fs10052-013-2677-2'
info.comment = ''
info.supersedes ='CMS-PAS-SUS-12-023'
info.supersededBy =''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked =''
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.source = "CMS"
T2tt.massConstraint = None
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked =''
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0'], ['dm <= 169.0']]
T2ttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure = '20a'
T2tt_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/fig20a.pdf'
T2tt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13011'
T2tt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13011'
T2tt_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased.root'
T2tt_1.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased.C'
T2tt_1.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2tt_exclusion.dat', 'orig/topneutralino_cutbased.root'],
                 dataFormats= ['txt', 'root'],objectNames= ['None', 'xsec_upperlimit'])
T2ttoff.addMassPlane(T2tt_1)

#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.checked =''
T6bbWW.constraint ="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription = "None"
T6bbWW.condition ="None"
T6bbWW.source = "CMS"
T6bbWW.massConstraint = None
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.checked =''
T6bbWWoff.constraint ="2.3*([[['b'],['L','nu']],[['b'],['jet','jet']]])"
T6bbWWoff.conditionDescription ="None"
T6bbWWoff.condition ="None"
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T6bbWW_1 = T6bbWW.addMassPlane(2*[[x, x*0.5+(1.-0.5)*y, y]])
T6bbWW_1.figure = '10c'
T6bbWW_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/fig10c.pdf'
T6bbWW_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/'
T6bbWW_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/bottomchargino_x50_cutbased.root'
T6bbWW_1.dataUrl =''
T6bbWW_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/bottomchargino_x50_cutbased.root'
T6bbWW_1.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/Obs_BDT_x05.txt', 'orig/bottomchargino_x50_BDT.root'],
                 dataFormats= ['txt', 'root'],objectNames= ['None', 'xsec_upperlimit'])
T6bbWWoff.addMassPlane(T6bbWW_1)
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane(2*[[x, x*0.25+(1.-0.25)*y, y]])
T6bbWW_2.figure ='10b'
T6bbWW_2.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/fig10b.pdf'
T6bbWW_2.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/'
T6bbWW_2.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/bottomchargino_x25_BDT.root'
T6bbWW_2.dataUrl =''
T6bbWW_2.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/bottomchargino_x25_BDT.root'
T6bbWW_2.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/Obs_BDT_x025.txt', 'orig/bottomchargino_x25_BDT.root'],
                 dataFormats= ['txt', 'root'],objectNames= ['None', 'xsec_upperlimit'])
T6bbWWoff.addMassPlane(T6bbWW_2)
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane(2*[[x, x*0.75+(1.-0.75)*y, y]])
T6bbWW_3.figure = '10d'
T6bbWW_3.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/fig10d.pdf'
T6bbWW_3.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/'
T6bbWW_3.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/bottomchargino_x75_BDT.root'
T6bbWW_3.dataUrl =''
T6bbWW_3.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/bottomchargino_x75_BDT.root'
T6bbWW_3.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/Obs_BDT_x075.txt', 'orig/bottomchargino_x75_BDT.root'],
                 dataFormats= ['txt', 'root'],objectNames= ['None', 'xsec_upperlimit'])
T6bbWWoff.addMassPlane(T6bbWW_3)



databaseCreator.create()
