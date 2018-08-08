#!/usr/bin/env python3

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
type = str)
argparser.add_argument ('-smodelsPath', '--smodelsPath',
help = 'path to the package smodels_utils',\
type = str)
argparser.add_argument ('-no', '--noUpdate',
help = 'do not update the lastUpdate field.',\
action= "store_true" )
argparser.add_argument ('-r', '--resetValidation',
help = 'reset the validation flag',\
action= "store_true" )
argparser.add_argument('-v', '--verbose',
help='specifying the level of verbosity (error, warning, info, debug)',
default = 'warning', type = str)
argparser.add_argument('-n', '--ntoys',
help='number of Monte Carlo toys',
default = 50000, type = int)

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
from smodels_utils.dataPreparation.commandlineArgs import setEnv
setEnv ( args )
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.datasetCreation import DatasetsFromLatex,createAggregationList
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z

# DataSetInput.ntoys = 50

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-16-050-agg')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/index.html'
info.sqrts = '13.0*TeV'
info.lumi = 35.9
info.prettyName = 'hadronic top tagging'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1719.11188'
info.comment = ''
info.implementedBy = 'HARG'
info.contact = 'SModelS'
info.supersedes = "CMS-PAS-SUS-16-050"
## we dont yet perform combinations
# info.combinableWith = "CMS-PAS-SUS-16-052, CMS-PAS-SUS-16-024, CMS-PAS-SUS-16-022"

max_datasets = None


## all 84 signal regions take around 16s on commodore??

## corr>.4: 60 aggregate regions, agreements 98/97/97. Takes around 30s(?) per point on wnouc.
# aggregation = [[1, 6, 10, 11, 22, 38], [2], [3], [4], [5], [7, 12, 13, 16, 17, 18, 21, 28, 29, 30, 32, 33, 34, 36, 37, 58, 68, 69, 70, 82], [8], [9], [14], [15], [19], [20], [23], [24], [25], [26], [27], [31], [35], [39], [40], [41], [42], [43], [44], [45], [46], [47], [48], [49], [50], [51], [52], [53], [54], [55], [56], [57], [59], [60], [61], [62], [63], [64], [65], [66], [67], [71], [72], [73], [74], [75], [76], [77], [78], [79], [80], [81], [83], [84]]

## corr>.4: remove 0-points (23,25,41,77), 56 aggregate regions, agreements 98/97/97. Takes around 10s per point on higgs. overexcludes for one case (T2tt?)
aggregation = [[1, 6, 10, 11, 22, 38], [2], [3], [4], [5], [7, 12, 13, 16, 17, 18, 21, 28, 29, 30, 32, 33, 34, 36, 37, 58, 68, 69, 70, 82], [8], [9], [14], [15], [19], [20], [24], [26], [27], [31], [35], [39], [40], [42], [43], [44], [45], [46], [47], [48], [49], [50], [51], [52], [53], [54], [55], [56], [57], [59], [60], [61], [62], [63], [64], [65], [66], [67], [71], [72], [73], [74], [75], [76], [78], [79], [80], [81], [83], [84]]

## corr>.4: remove low-points (5,8,23,24,25,35,41,45,48,52,59,62,66,67,74,77,78,80,84), 41 aggregate regions, agreements 91/97/94. Takes around 5.5s per point on commodore.
# aggregation = [[1, 6, 10, 11, 22, 38], [2], [3], [4], [7, 12, 13, 16, 17, 18, 21, 28, 29, 30, 32, 33, 34, 36, 37, 58, 68, 69, 70, 82], [9], [14], [15], [19], [20], [26], [27], [31], [39], [40], [42], [43], [44], [46], [47], [49], [50], [51], [53], [54], [55], [56], [57], [60], [61], [63], [64], [65], [71], [72], [73], [75], [76], [79], [81], [83]]

## corr>.4: move low-points (5,8,23,24,25,35,41,45,48,52,59,62,66,67,74,77,78,80,84) to its own region, 42 aggregate regions, agreements 97/90/93. Takes around 7s per point on higgs.
# aggregation = [[1, 6, 10, 11, 22, 38], [2], [3], [4], [7, 12, 13, 16, 17, 18, 21, 28, 29, 30, 32, 33, 34, 36, 37, 58, 68, 69, 70, 82], [9], [14], [15], [19], [20], [26], [27], [31], [39], [40], [42], [43], [44], [46], [47], [49], [50], [51], [53], [54], [55], [56], [57], [60], [61], [63], [64], [65], [71], [72], [73], [75], [76], [79], [81], [83], [5,8,23,24,25,35,41,45,48,52,59,62,66,67,74,77,78,80,84] ]

## corr>.3: 
# aggregation = [[1, 2, 6, 10, 11, 15, 22, 23, 38, 53, 54, 59, 63], [3, 4, 9, 51], [5, 7, 12, 13, 16, 17, 18, 21, 24, 28, 29, 30, 32, 33, 34, 35, 36, 37, 57, 58, 68, 69, 70, 82], [14], [19], [20], [26], [27], [31], [39], [40], [42], [43], [44], [46], [47], [49], [50], [55], [56], [60], [61], [64], [65], [71], [72], [73], [75], [76], [79], [81], [83]]
## low pointers [1, 2, 5, 6, 7, 8, 11, 12, 13, 23, 24, 25, 29, 33, 34, 35, 36, 37, 38, 41, 45, 48, 52, 59, 62, 66, 67, 69, 74, 77, 78, 80, 84]

info.createCovarianceMatrix ( "orig/CMS-SUS-16-050_Figure-aux_001.root", \
        "total_covar", addOrder=True, max_datasets = max_datasets, aggregate=aggregation )

creator = DatasetsFromLatex ( "orig/tables.tex", max_datasets = max_datasets,
         ds_name = "t#1b#2MT2#3MET#4", aggregate=aggregation )

for ctr,dataset in enumerate(creator):
    #+++++++ next txName block ++++++++++++++
    T2tt = dataset.addTxName('T2tt')
    T2tt.checked =''
    T2tt.constraint = "[[['t']],[['t']]]"
    T2tt.conditionDescription = None
    T2tt.condition =None
    T2tt.source ='CMS'
    #+++++++ next mass plane block ++++++++++++++
    T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
    objNames = 'acc_bin%d' % ctr
    if aggregation != None:
        objNames = [ 'acc_bin%d' % (i-1) for i in aggregation[ctr] ]
    #----exclusion source----
    T2tt_1.addSource( 'efficiencyMap', 'orig/acc_maps_T2tt.root', 'root',
                      objectName = objNames )
    T2tt_1.addSource( 'obsExclusion', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ObsLim' )
    T2tt_1.addSource( 'obsExclusionP1', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ObsLimSup' )
    T2tt_1.addSource( 'obsExclusionM1', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ObsLimSdn' )
    T2tt_1.addSource( 'expExclusion', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ExpLim' )
    T2tt_1.addSource( 'expExclusionP1', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ExpLimSup' )
    T2tt_1.addSource( 'expExclusionM1', 'orig/CMS-SUS-16-050_Figure_008.root', 'root', objectName = 'ExpLimSdn' )
    T2tt_1.dataUrl = "cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/acc_maps_T2tt.root"
    T2tt_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure_008.png"

    #+++++++ next txName block ++++++++++++++
    T1tttt = dataset.addTxName('T1tttt')
    T1tttt.checked =''
    T1tttt.constraint = "[[['t','t']],[['t','t']]]"
    T1tttt.conditionDescription = None
    T1tttt.condition =None
    T1tttt.source ='CMS'
    #+++++++ next mass plane block ++++++++++++++
    T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
    #----exclusion source----
    T1tttt_1.addSource( 'efficiencyMap', 'orig/acc_maps_T1tttt.root', 'root',
                      objectName = objNames )
    T1tttt_1.addSource( 'obsExclusion', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ObsLim' )
    T1tttt_1.addSource( 'obsExclusionP1', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ObsLimSup' )
    T1tttt_1.addSource( 'obsExclusionM1', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ObsLimSdn' )
    T1tttt_1.addSource( 'expExclusion', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ExpLim' )
    T1tttt_1.addSource( 'expExclusionP1', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ExpLimSup' )
    T1tttt_1.addSource( 'expExclusionM1', 'orig/CMS-SUS-16-050_Figure_009-a.root', 'root', objectName = 'ExpLimSdn' )
    T1tttt_1.dataUrl = "cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/acc_maps_T1tttt.root"
    T1tttt_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure_009-d.png"

    #+++++++ next txName block ++++++++++++++
    T5tctc = dataset.addTxName('T5tctc')
    T5tctc.checked =''
    T5tctc.constraint = "[[['t'],['c']],[['t'],['c']]]"
    T5tctc.conditionDescription = None
    T5tctc.condition =None
    T5tctc.source ='CMS'
    #+++++++ next mass plane block ++++++++++++++
    T5tctc_1 = T5tctc.addMassPlane([[x,y+20,y]]*2)
    #----exclusion source----
    T5tctc_1.addSource( 'efficiencyMap', 'orig/acc_maps_T5ttcc.root', 'root',
                      objectName = objNames )
    T5tctc_1.addSource( 'obsExclusion', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ObsLim' )
    T5tctc_1.addSource( 'obsExclusionP1', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ObsLimSup' )
    T5tctc_1.addSource( 'obsExclusionM1', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ObsLimSdn' )
    T5tctc_1.addSource( 'expExclusion', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ExpLim' )
    T5tctc_1.addSource( 'expExclusionP1', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ExpLimSup' )
    T5tctc_1.addSource( 'expExclusionM1', 'orig/CMS-SUS-16-050_Figure_009-d.root', 'root', objectName = 'ExpLimSdn' )
    T5tctc_1.dataUrl = "cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/acc_maps_T5ttcc.root"
    T5tctc_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-050/CMS-SUS-16-050_Figure-aux_005.png"

databaseCreator.create()
