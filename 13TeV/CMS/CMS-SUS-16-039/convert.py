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
info = MetaInfoInput('CMS-SUS-16-039')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'Multilepton EWK searches'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1709.05406'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Submitted to J. High Energy Phys.'
info.comment = 'Moriond 2017. Negative values in TChiChipmSlepL and TChiChipmSlepStau root files. Added TChiChipmStauStau.'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++

TChiWZ=dataset.addTxName('TChiWZ')
TChiWZ.checked=''
TChiWZ.constraint="[[['W']],[['Z']]]"
TChiWZ.condition=None
TChiWZ.conditionDescription = None
TChiWZ.source="CMS"

#offshell txName block

TChiWZoff=dataset.addTxName('TChiWZoff')
TChiWZoff.checked=''
TChiWZoff.constraint="71.*([[['mu+','mu-']],[['l','nu']]] + [[['e+','e-']],[['l','nu']]])"
TChiWZoff.condition = "cGtr([[['mu+','mu-']],[['l','nu']]],[[['e+','e-']],[['l','nu']]])"
TChiWZoff.massConstraint = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription=None
TChiWZoff.source="CMS"

#++++++next mass plane block+++++++++

TChiWZ_1 = TChiWZ.addMassPlane(2*[[x,y]])
TChiWZ_1.figure='Fig. 18-a'
TChiWZ_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_018-a.png'
TChiWZ_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_018-a.root'
TChiWZ_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_018-a.root'
TChiWZ_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_018-a.root'
TChiWZ_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-039_Figure_018-a.root','orig/CMS-SUS-16-039_Figure_018-a.root','orig/CMS-SUS-16-039_Figure_018-a.root','orig/CMS-SUS-16-039_Figure_018-a.root','orig/CMS-SUS-16-039_Figure_018-a.root','orig/CMS-SUS-16-039_Figure_018-a.root','orig/CMS-SUS-16-039_Figure_018-a.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['TChiWZ;1','TChiWZ;1','TChiWZ;1','TChiWZ;1','TChiWZ;1','TChiWZ;1','TChiWZ;1'],
                    indices= [4, 6, 5, 7, 9, 8, 2],units=[None,None,None,None,None,None,'pb'])

TChiWZoff.addMassPlane(TChiWZ_1)

#++++++next txName block+++++++++++++++

TChiWH=dataset.addTxName('TChiWH')
TChiWH.checked=''
TChiWH.constraint="[[['W']],[['higgs']]]"
TChiWH.condition=None
TChiWH.conditionDescription=None
TChiWH.source="CMS"


#++++++next mass plane block++++++++
TChiWH_1 = TChiWH.addMassPlane(2*[[x,y]])
TChiWH_1.figure='Fig. 18-b'
TChiWH_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_018-b.png'
TChiWH_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_018-b.root'
TChiWH_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_018-b.root'
TChiWH_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_018-b.root'
TChiWH_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-039_Figure_018-b.root','orig/CMS-SUS-16-039_Figure_018-b.root','orig/CMS-SUS-16-039_Figure_018-b.root','orig/CMS-SUS-16-039_Figure_018-b.root','orig/CMS-SUS-16-039_Figure_018-b.root','orig/CMS-SUS-16-039_Figure_018-b.root','orig/CMS-SUS-16-039_Figure_018-b.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1'],
                    indices= [4, 6, 5, 7, 9, 8, 2],units=[None,None,None,None,None,None,'pb'])

#+++++next txName block+++++++++++++++
TChiChipmStauStau=dataset.addTxName('TChiChipmStauStau')
TChiChipmStauStau.checked=''
TChiChipmStauStau.constraint="[[['ta+'],['ta-']],[['nu'],['ta']]]+[[['ta-'],['ta+']],[['nu'],['ta']]]"
TChiChipmStauStau.condition=None
TChiChipmStauStau.conditionDescription=None
TChiChipmStauStau.source = "CMS"
TChiChipmStauStau.massConstraint=None

#++++++next txName block+++++++++++++++

TChiChipmStauStau_1=TChiChipmStauStau.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
TChiChipmStauStau_1.figure='Fig. 17'
TChiChipmStauStau_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_017.png'
TChiChipmStauStau_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_017.root'
TChiChipmStauStau_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_017.root'
TChiChipmStauStau_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_017.root'
TChiChipmStauStau_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-039_Figure_017.root','orig/CMS-SUS-16-039_Figure_017.root','orig/CMS-SUS-16-039_Figure_017.root','orig/CMS-SUS-16-039_Figure_017.root','orig/CMS-SUS-16-039_Figure_017.root','orig/CMS-SUS-16-039_Figure_017.root','orig/CMS-SUS-16-039_Figure_017.root'], dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],
                        objectNames=['TChiSlepSnu_TD;1','TChiSlepSnu_TD;1','TChiSlepSnu_TD;1','TChiSlepSnu_TD;1','TChiSlepSnu_TD;1','TChiSlepSnu_TD;1','TChiSlepSnu_TD;1'],
                    indices= [4,6, 5, 7, 9, 8, 2],units=[None,None,None,None,None,None,'pb'])


#+++++next txName block+++++++++++++++
TChiChipmSlepStau=dataset.addTxName('TChiChipmSlepStau')
TChiChipmSlepStau.checked=''
TChiChipmSlepStau.constraint="[[['L'],['L']],[['nu'],['ta']]]"
TChiChipmSlepStau.conditionDescription=None
TChiChipmSlepStau.condition =  "cGtr([[['L'],['L']],[['nu'],['ta']]],3.*[[['ta'],['ta']],[['nu'],['ta']]]);cGtr(3.*[[['mu'],['mu']],[['nu'],['ta']]],[[['L'],['L']],[['nu'],['ta']]])"
TChiChipmSlepStau.source = "CMS"
TChiChipmSlepStau.massConstraint=None

#++++++next txName block+++++++++++++++
TChiChipmSlepStau_1=TChiChipmSlepStau.addMassPlane(2*[[x,0.05*x+0.95*y,y]])
TChiChipmSlepStau_1.figure='Fig. 16-a'
TChiChipmSlepStau_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-a.png'
TChiChipmSlepStau_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-a.root'
TChiChipmSlepStau_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-a.root'
TChiChipmSlepStau_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-a.root'
TChiChipmSlepStau_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-039_Figure_016-a.root','orig/CMS-SUS-16-039_Figure_016-a.root','orig/CMS-SUS-16-039_Figure_016-a.root','orig/CMS-SUS-16-039_Figure_016-a.root','orig/CMS-SUS-16-039_Figure_016-a.root','orig/CMS-SUS-16-039_Figure_016-a.root','orig/CMS-SUS-16-039_Figure_016-a.root'], dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],
                        objectNames=['TChiSlepSnu_TE_05;1','TChiSlepSnu_TE_05;1','TChiSlepSnu_TE_05;1','TChiSlepSnu_TE_05;1','TChiSlepSnu_TE_05;1','TChiSlepSnu_TE_05;1',
'TChiSlepSnu_TE_05;1'],indices= [4,6, 5, 7, 9, 8, 2],units=[None,None,None,None,None,None,'pb'])


#++++++next txName block+++++++++++++++
TChiChipmSlepStau_2=TChiChipmSlepStau.addMassPlane(2*[[x,0.95*x+0.05*y,y]])
TChiChipmSlepStau_2.figure='Fig. 16-c'
TChiChipmSlepStau_2.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-c.png'
TChiChipmSlepStau_2.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-c.png'
TChiChipmSlepStau_2.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-c.png'
TChiChipmSlepStau_2.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-c.png'
TChiChipmSlepStau_2.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-039_Figure_016-c.root','orig/CMS-SUS-16-039_Figure_016-c.root','orig/CMS-SUS-16-039_Figure_016-c.root','orig/CMS-SUS-16-039_Figure_016-c.root','orig/CMS-SUS-16-039_Figure_016-c.root','orig/CMS-SUS-16-039_Figure_016-c.root','orig/CMS-SUS-16-039_Figure_016-c.root'], dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],
                        objectNames=['TChiSlepSnu_TE_95;1','TChiSlepSnu_TE_95;1','TChiSlepSnu_TE_95;1','TChiSlepSnu_TE_95;1','TChiSlepSnu_TE_95;1','TChiSlepSnu_TE_95;1',
'TChiSlepSnu_TE_95;1'],
                    indices= [4,6, 5, 7, 9, 8, 2],units=[None,None,None,None,None,None,'pb'])

#++++++next txName block+++++++++++++++
TChiChipmSlepStau_3=TChiChipmSlepStau.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
TChiChipmSlepStau_3.figure='Fig. 16-b'
TChiChipmSlepStau_3.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-b.png'
TChiChipmSlepStau_3.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-b.root'
TChiChipmSlepStau_3.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-b.root'
TChiChipmSlepStau_3.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_016-b.root'
TChiChipmSlepStau_3.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-039_Figure_016-b.root','orig/CMS-SUS-16-039_Figure_016-b.root','orig/CMS-SUS-16-039_Figure_016-b.root','orig/CMS-SUS-16-039_Figure_016-b.root','orig/CMS-SUS-16-039_Figure_016-b.root','orig/CMS-SUS-16-039_Figure_016-b.root','orig/CMS-SUS-16-039_Figure_016-b.root'], dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],
                        objectNames=['TChiSlepSnu_TE;1','TChiSlepSnu_TE;1','TChiSlepSnu_TE;1','TChiSlepSnu_TE;1','TChiSlepSnu_TE;1','TChiSlepSnu_TE;1','TChiSlepSnu_TE;1'],
                    indices= [4,6, 5, 7, 9, 8, 2],units=[None,None,None,None,None,None,'pb'])

#++++++next txName block+++++++++++++++

TChiChipmSlepL=dataset.addTxName('TChiChipmSlepL')
TChiChipmSlepL.checked=''
TChiChipmSlepL.constraint="2*(([[['e+'],['e-']],[['L'],['nu']]]+[[['e+'],['e-']],[['nu'],['L']]]+[[['e-'],['e+']],[['L'],['nu']]]+[[['e-'],['e+']],[['nu'],['L']]])+([[['mu+'],['mu-']],[['L'],['nu']]]+[[['mu+'],['mu-']],[['nu'],['L']]]+[[['mu-'],['mu+']],[['L'],['nu']]]+[[['mu-'],['mu+']],[['nu'],['L']]])+([[['ta+'],['ta-']],[['L'],['nu']]]+[[['ta+'],['ta-']],[['nu'],['L']]]+[[['ta-'],['ta+']],[['L'],['nu']]]+[[['ta-'],['ta+']],[['nu'],['L']]]))"
TChiChipmSlepL.condition="cSim(([[['e+'],['e-']],[['L'],['nu']]]+[[['mu+'],['mu-']],[['L'],['nu']]]+[[['ta+'],['ta-']],[['L'],['nu']]]+[[['e-'],['e+']],[['L'],['nu']]]+[[['mu-'],['mu+']],[['L'],['nu']]]+[[['ta-'],['ta+']],[['L'],['nu']]]),([[['e+'],['e-']],[['nu'],['L']]]+[[['e-'],['e+']],[['nu'],['L']]]+[[['mu+'],['mu-']],[['nu'],['L']]]+[[['mu-'],['mu+']],[['nu'],['L']]]+[[['ta+'],['ta-']],[['nu'],['L']]]+[[['ta-'],['ta+']],[['nu'],['L']]]));cGtr(3*(([[['mu+'],['mu-']],[['L'],['nu']]]+[[['mu+'],['mu-']],[['nu'],['L']]]+[[['mu-'],['mu+']],[['L'],['nu']]]+[[['mu-'],['mu+']],[['nu'],['L']]])),([[['e+'],['e-']],[['L'],['nu']]]+[[['e+'],['e-']],[['nu'],['L']]]+[[['e-'],['e+']],[['L'],['nu']]]+[[['e-'],['e+']],[['nu'],['L']]]+[[['mu+'],['mu-']],[['L'],['nu']]]+[[['mu+'],['mu-']],[['nu'],['L']]]+[[['mu-'],['mu+']],[['L'],['nu']]]+[[['mu-'],['mu+']],[['nu'],['L']]]+[[['ta+'],['ta-']],[['L'],['nu']]]+[[['ta+'],['ta-']],[['nu'],['L']]]));cGtr(([[['e+'],['e-']],[['L'],['nu']]]+[[['e+'],['e-']],[['nu'],['L']]]+[[['e-'],['e+']],[['L'],['nu']]]+[[['e-'],['e+']],[['nu'],['L']]]+[[['mu+'],['mu-']],[['L'],['nu']]]+[[['mu+'],['mu-']],[['nu'],['L']]]+[[['mu-'],['mu+']],[['L'],['nu']]]+[[['mu-'],['mu+']],[['nu'],['L']]]+[[['ta+'],['ta-']],[['L'],['nu']]]+[[['ta+'],['ta-']],[['nu'],['L']]]+[[['ta-'],['ta+']],[['L'],['nu']]]+[[['ta-'],['ta+']],[['nu'],['L']]]),3*([[['ta+'],['ta-']],[['L'],['nu']]]+[[['ta+'],['ta-']],[['nu'],['L']]]+[[['ta-'],['ta+']],[['L'],['nu']]]+[[['ta-'],['ta+']],[['nu'],['L']]]))"
TChiChipmSlepL.conditionDescription=None
TChiChipmSlepL.source="CMS"
TChiChipmSlepL.massConstraint=None
TChiChipmSlepL.round_to = 6

#++++++ mass plane block++++++++
TChiChipmSlepL_1 = TChiChipmSlepL.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
TChiChipmSlepL_1.figure='Fig. 14'
TChiChipmSlepL_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_014.png'
TChiChipmSlepL_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_014.root'
TChiChipmSlepL_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_014.root'
TChiChipmSlepL_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_014.root'
TChiChipmSlepL_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-039_Figure_014.root','orig/CMS-SUS-16-039_Figure_014.root','orig/CMS-SUS-16-039_Figure_014.root','orig/CMS-SUS-16-039_Figure_014.root','orig/CMS-SUS-16-039_Figure_014.root','orig/CMS-SUS-16-039_Figure_014.root','orig/CMS-SUS-16-039_Figure_014.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],
                    objectNames=['TChiSlepSnu_FD;1','TChiSlepSnu_FD;1','TChiSlepSnu_FD;1','TChiSlepSnu_FD;1','TChiSlepSnu_FD;1','TChiSlepSnu_FD;1','TChiSlepSnu_FD;1'],
                    indices=[4,6, 5, 7, 9, 8, 2] ,units=[None,None,None,None,None,None,'pb'])


#++++++next mass plane block++++++++
TChiChipmSlepL_2 = TChiChipmSlepL.addMassPlane(2*[[x,0.05*x+0.95*y,y]])
TChiChipmSlepL_2.figure='Fig. 15-a'
TChiChipmSlepL_2.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_015-a.png'
TChiChipmSlepL_2.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_015-a.root'
TChiChipmSlepL_2.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_015-a.root'
TChiChipmSlepL_2.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_015-a.root'
TChiChipmSlepL_2.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-039_Figure_015-a.root','orig/CMS-SUS-16-039_Figure_015-a.root','orig/CMS-SUS-16-039_Figure_015-a.root','orig/CMS-SUS-16-039_Figure_015-a.root','orig/CMS-SUS-16-039_Figure_015-a.root','orig/CMS-SUS-16-039_Figure_015-a.root','orig/CMS-SUS-16-039_Figure_015-a.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['TChiSlepSnu_FD_05;1','TChiSlepSnu_FD_05;1','TChiSlepSnu_FD_05;1','TChiSlepSnu_FD_05;1','TChiSlepSnu_FD_05;1','TChiSlepSnu_FD_05;1','TChiSlepSnu_FD_05;1'],
                    indices= [4,6, 5, 7, 9, 8, 2],units=[None,None,None,None,None,None,'pb'])

#++++++next mass plane block++++++++
TChiChipmSlepL_3 = TChiChipmSlepL.addMassPlane(2*[[x,0.95*x+0.05*y,y]])
TChiChipmSlepL_3.figure='Fig. 15-b'
TChiChipmSlepL_3.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_015-b.png'
TChiChipmSlepL_3.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_015-b.root'
TChiChipmSlepL_3.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_015-b.root'
TChiChipmSlepL_3.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-039/CMS-SUS-16-039_Figure_015-b.root'
TChiChipmSlepL_3.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-16-039_Figure_015-b.root','orig/CMS-SUS-16-039_Figure_015-b.root','orig/CMS-SUS-16-039_Figure_015-b.root','orig/CMS-SUS-16-039_Figure_015-b.root','orig/CMS-SUS-16-039_Figure_015-b.root','orig/CMS-SUS-16-039_Figure_015-b.root','orig/CMS-SUS-16-039_Figure_015-b.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],
                    objectNames=['TChiSlepSnu_FD_95;1','TChiSlepSnu_FD_95;1','TChiSlepSnu_FD_95;1','TChiSlepSnu_FD_95;1','TChiSlepSnu_FD_95;1','TChiSlepSnu_FD_95;1','TChiSlepSnu_FD_95;1'],
                    indices=[4,6, 5, 7, 9, 8, 2] ,units=[None,None,None,None,None,None,'pb'])

databaseCreator.create()
