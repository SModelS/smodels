#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt, txname.txt, twiki.txt and sms.py.

.. moduleauthor: Wolfgang Waltenberger

"""
import sys
import os
import argparse
import types

argparser = argparse.ArgumentParser(description = \
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath', \
help = 'path to the package smodels_utils',\
type = types.StringType)
args = argparser.parse_args()

if args.utilsPath:
    utilsPath = args.utilsPath
else:
    databaseRoot = '../../../../'
    sys.path.append(os.path.abspath(databaseRoot))
    from utilsPath import utilsPath
    utilsPath = databaseRoot + utilsPath

sys.path.append(os.path.abspath(utilsPath))
from smodels_utils.dataPreparation.inputObjects import TxNameInput, MetaInfoInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.origPlotObjects import x, y


#+++++++ global info block ++++++++++++++

import os, glob
dir=os.getcwd()
print dir
pos1=dir.find("ATLAS/")+6                                                                                                       
pos=dir.find("-ANA")
expid = dir[pos1:pos]  
print "expid=",expid
signalregion=dir[pos+1:]
print "signalregion=",signalregion

info = MetaInfoInput(expid)
info.signalRegion = signalregion
info.url ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/%s/' % expid 
info.sqrts = 8
info.prettyName = ''
info.private = True
info.arxiv = ''
info.contact ='fastlim'
info.publication = ''
info.lumi = 1.0 
info.comment = 'created from fastlim-1.0'
info.supersededBy = ''
info.implementedBy = ''


constraints =  { "T2tt": "[[['t+']],[['t-']]]", "T2bb": "[[['b']],[['b']]]", 
                 "T2": "[[['jet']],[['jet']]]",
                 "T2bt": "[[['b']],[['t']]]", "T1tttt": "[[['t+','t-']],[['t+','t-']]]",
                 "T5tttt": "[[['t+'],['t-']],[['t+'],['t-']]]",
                 "T5bbbb": "[[['b'],['b']],[['b'],['b']]]",
                 "T5bbbt": "[[['b'],['b']],[['b'],['t']]]",
                 "T1bbtt": "[[['b','b']],[['t','t']]]", "T1btbt": "[[['b','t']],[['b','t']]]",
                 "T1bbqq": "[[['b','b']],[['jet','jet']]]", "T1bbbb": "[[['b','b']],[['b','b']]]",
                 "T1bbbt": "[[['b','b']],[['b','t']]]", "T5btbt": "[[['b'],['t']],[['b'],['t']]]",
                 "T5tbtb": "[[['t'],['b']],[['t'],['b']]]", "T5tbtt": "[[['t'],['b']],[['t'],['t']]]",
                 "TGQqtt": "[[['jet']],[['t+','t-']]]", "TGQ": "[[['jet']],[['jet','jet']]]",
                 "TGQbtq": "[[['b','t']],[['jet']]]", "TGQbbq": "[[['b','b']],[['jet']]]",
                 "T1btqq": "[[['b','t']],[['jet','jet']]]", "T1qqtt": "[[['jet','jet']],[['t','t']]]",
                 "T1bttt": "[[['b','t']],[['t','t']]]", "T1": "[[['jet','jet']],[['jet','jet']]]" }

#+++++++ next txName block ++++++++++++++
#T1 = TxName('T1')
#T1.on.constraint = constraints["T1"]
#T1.on.conditionDescription = None
#T1.on.condition = None
#
##+++++++ next mass plane block ++++++++++++++
#T1_1 = T1.addMassPlane(motherMass = x , lspMass = y )
##----figure----
#T1_1.figure = 'fig_07a.pdf'
#T1_1.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_07a.pdf'
## ----limit source----
### T1_1.obsUpperLimit.setSource( './orig/T1.effi', 'txt', objectName = None, index = None )
#T1_1.efficiencyMap.setSource( './orig/T1.effi', 'txt', objectName = None, index = None )
## T1_1.obsUpperLimit.unit = 'fb'
## T1_1.expUpperLimit.setSource( path, type, objectName = None, index = None )
## ----exclusion source----
#T1_1.obsExclusion.setSource( './orig/T1_exc.dat', 'txt', objectName = None, index = None )
##T1_1.obsExclusionM1.setSource( path, type, objectName = None, index = None )
#T1_1.obsExclusionP1.setSource( path, type, objectName = None, index = None )
#T1_1.expExclusion.setSource( path, type, objectName = None, index = None )
#T1_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T1_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
#----global url settings ----
#T1_1.dataUrl = 
#T1_1.histoDataUrl = 
#----limit url settings ----
# T1_1.obsUpperLimit.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13007/limits_model_A.txt'
#T1_1.expectedlimit.dataUrl =
#----exclusion url settings ----
#T1_1.exclusionDataUrl =
#T1_1.exclusion.dataUrl =
#T1_1.exclusionM1.dataUrl =
#T1_1.exclusionP1.dataUrl =
#T1_1.expectedExclusion.dataUrl =
#T1_1.expectedExclusionM1.dataUrl =
#T1_1.expectedExclusionP1.dataUrl =

figure={}
figure["T1"]='fig_07a.pdf'

figureUrl={}
figureUrl["T1"]='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-047/fig_07a.pdf'

for i in os.listdir("orig/"):
    if i[-5:]!=".effi": continue
    txname=i[:-5]
    print txname 
    tmp= TxNameInput ( txname )
    tmp.on.constraint = constraints[txname]
    tmp.on.conditionDescription=None
    tmp.on.condition = None
    if i[:2] in [ "T5", "T6" ]:
        tmp.globalEfficiencyMap.setSource ( './orig/%s.effi' % txname, 'effi', objectName = None, index = None )
        continue
    if txname=="TGQ":
        tmp_1 = tmp.addMassPlane (motherMass = x , lspMass = y )
    else:
        tmp_1 = tmp.addMassPlane (motherMass = x , lspMass = y )
    tmp_1.efficiencyMap.setSource( './orig/%s.effi' % txname, 'effi', objectName = None, index = None ) 
    if os.path.exists ( './orig/%s_excl.dat' % txname ):
        tmp_1.obsExclusion.setSource( './orig/%s_excl.dat' % txname, 'txt', objectName = None, index = None )
    if txname in figure:
        tmp_1.figure = figure[txname]
    if txname in figureUrl:
        tmp_1.figureUrl = figureUrl[txname]

def translate ( filename ):
    ff=open(filename)
    lines=ff.readlines()[3:]
    ff.close()
    newfilename=filename.replace("./orig/","").replace(".effi",".hlp")
    w=open(newfilename,"w")
    for line in lines:
        a=line.split()
        f=[ float(a[0]), float(a[1]), float(a[2]), float(a[3]) ]
        line="[[[%.1f*GeV, %.1f*GeV, %.1f*GeV], [%.1f*GeV, %.1f*GeV, %.1f*GeV]], %f]," % \
              ( f[0], f[1], f[2], f[0], f[1], f[2], f[3] )
        w.write ( line +"\n" )
    w.close()


#T5tttt = TxName ( "T5tttt" )
#T5tttt.on.constraint = "[[['t+'],['t-']],[['t+'],['t-']]]"
#T5tttt.on.conditionDescription = None
#T5tttt.on.condition = None
## T5tttt.globalEfficiencyMap.setSource ( './orig/T5tttt.effi', 'effi', objectName = None, index = None )
#translate ( './orig/T5tttt.effi' )

# ----limit source----
## T1_1.obsUpperLimit.setSource( './orig/T1.effi', 'txt', objectName = None, index = None )
#        tmp_1.efficiencyMap.setSource( './orig/T1.effi', 'txt', objectName = None, index = None )
# T1_1.obsUpperLimit.unit = 'fb'
# T1_1.expUpperLimit.setSource( path, type, objectName = None, index = None )
# ----exclusion source----
        

#T2tt = TxName('T2tt')
#T2tt.on.constraint = "[[['t']],[['t']]]"
#T2tt.on.conditionDescription = None
#T2tt.on.condition = None
#T2tt_1 = T2tt.addMassPlane(motherMass = x , lspMass = y )
#T2tt_1.efficiencyMap.setSource( './orig/T2tt.effi', 'effi', objectName = None, index = None )

databaseCreator.infoFileDirectory="./"
databaseCreator.create()

import os
os.unlink ("info.txt")
