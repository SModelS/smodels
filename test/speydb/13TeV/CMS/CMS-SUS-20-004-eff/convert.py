#!/usr/bin/env python3

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import types
import math
from smodels_utils.dataPreparation.argParser import getParserArgs
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.datasetCreation import aggregateDataSets, \
         createAggregationOrder
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z
getParserArgs()

#databaseCreator.ncpus = 1

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-20-004')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-004/'
info.sqrts = '13.0*TeV'
info.lumi = 137
info.prettyName  = '2 h(b b), EWK'
info.publication = 'JHEP 05 (2022) 014'
info.publicationDOI = 'https://doi.org/10.1007/JHEP05(2022)014'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.arxiv = 'https://arxiv.org/abs/2201.04206'
info.comment = "Datasets include 3rd momenta -> first usage of SLv2 (thanks to Bill Ford). Need to use prefit BG numbers; postfit ones give bad validation. For TChiHH we did not add the gravitino results (TChiHH-G model in the paper), as that would give us an unrealistic interpolation for mother masses between 800 and 1200 GeV."
info.private = False
info.implementedBy = 'WW'
#info.supersedes = "CMS-PAS-SUS-16-050"

aggregates = None
# aggregates = [[0, 1, 5, 9, 10, 14, 21, 22, 30, 37, 38, 41, 48, 52, 53, 58, 59, 62, 65], [2, 3, 8, 50, 60], [4, 6, 7, 11, 12, 15, 16, 17, 18, 20, 23, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 47, 55, 56, 57, 67, 68, 69, 81], [13, 19, 40, 43, 44, 51, 64, 66, 71], [24, 39], [42], [45], [46], [49], [54], [61], [63], [70, 75], [72, 77], [73], [74, 76], [78], [79], [80], [82], [83]]

datasets = []

dsnames = ['SR1', 'SR2', 'SR3', 'SR4', 'SR5', 'SR6', 'SR7', 'SR8','SR9', 'SR10', 'SR11', 'SR12', 'SR13', 'SR14', 'SR15', 'SR16', 'SR17', 'SR18', 'SR19', 'SR20', 'SR21', 'SR22']

nobs = [138., 91., 14., 3., 54., 38., 4., 0., 8., 2., 4., 0., 1., 3., 1., 1., 42., 6., 1., 4., 0., 0.]
nbg = [160.57, 90.439, 11.512, 2.7508, 53.529, 28.319, 2.5856, 2.6029, 5.0628, 2.1711, 0.063528, 0.88532, 2.6812, 1.2599, 0.41911, 0.67201, 33.567, 7.3285, 1.6546, 4.0289, 0.88189, 0.19966] # prefit
bgerr = [14, 9.7, 3.4, 2.3, 8.8, 5.6, 1.5, 2.4, 1.6, 0.79, 0.11, 1.42, 1.06, 0.62, 0.61, 1.10, 6.1, 2.0, 1.04, 1.5, 0.40, 0.21] # prefit
third_momenta = [456.69, 216.349, 20.4459, 11.1834, 223.787, 66.9344, 2.73344, 12.2239, 2.24935, 0.300774, 0.00170956, 3.08295, 0.772409, 0.173894, 0.24212, 1.43719, 96.3718, 3.88287, 0.900496, 2.30716, 0.0461215, 0.00956027] # prefit
postfit = False
if postfit:
    nbg = [149.7, 91.5, 12.8, 2.8, 54.1, 33.2, 3.2, 1.27, 5.9, 2.31, 0.72, 0.52, 2.58, 1.62, 1.16, 0.78, 37.0, 7.2, 1.50, 4.0, 0.74, 0.14] # postfit
    bgerr = [8.9, 6.9, 2.6, 1.4, 5.6, 4.2, 1.3, 0.98, 1.4, 0.73, 0.53, 0.65, 0.85, 0.65, 0.87, 0.76, 4.2, 1.5, 0.75, 1.2, 0.29, 0.13] # postfit
    third_momenta = [89.6133207755313, 52.01993867878129, 6.545475410175557, 1.8979767656009354, 37.73442845552516, 19.097367532338634, 1.2721691347255832, 0.7569841640613788, 1.1385623192095085, 0.2223261805522786, 0.1280992100309235, 0.2589021843675942, 0.34535607421498965, 0.17294582720232696, 0.5360461361387237, 0.408862185466996, 23.44710013517256, 1.1442286650030837, 0.2863944651158953, 0.9292135600842435, 0.01702728432752175, 0.0022959935123803834] # postfit

for i in range(22):
    dsname = f'SR{i+1}'
    dataset = DataSetInput( dsname )
    dataset.setInfo(dataType = 'efficiencyMap', dataId = dsname,
           observedN = nobs[i], expectedBG=nbg[i], bgError = bgerr[i],
           comment = dsnames[i], thirdMoment = third_momenta[i] )

    TChiHH = dataset.addTxName('TChiHH')
    TChiHH.validationTarball = "[[x, 0.0], [x, 0.0]]:TChiHHN3.tar.gz; [[x, 1.0], [x, 1.0]]:TChiHHN3.tar.gz; [[x, y], [x, y]]:TChiHHN3.tar.gz"
    TChiHH.checked = ''
    TChiHH.constraint = "[[['higgs']],[['higgs']]]"
    TChiHH.condition = None
    TChiHH.conditionDescription = None
    TChiHH.source = "CMS"
    #T2tt.massConstraint=[['dm>=338.0'],['dm>=338.0']]


    T5HH=dataset.addTxName('T5HH')
    T5HH.checked=''
    T5HH.constraint="[[['jet','jet'],['higgs']],[['jet','jet'],['higgs']]]"
    T5HH.condition=None
    T5HH.conditionDescription = None
    T5HH.source="CMS"

    TChiHH_1 = TChiHH.addMassPlane(2*[[x,y]])
    TChiHH_1.addSource( 'efficiencyMap', f'orig/SR{i+1}.root', 'root', objectName = f'Bin{i+1}SignalEfficiency_N1N2;1')
    TChiHH_1.efficiencyMap._unit = ".3364"
    TChiHH_1.xrange = "[200,500]"
    TChiHH_1.yrange = "[0,200]"
    TChiHH_1.addSource( 'obsExclusion', 'orig/TChiHHobservedmasslimitcurve.csv', 'csv' )
    TChiHH_1.addSource( 'expExclusion', 'orig/TChiHHexpectedmasslimitcurve.csv', 'csv')
    TChiHH_1.addSource( 'expExclusionP1', 'orig/TChiHHexpected+1s.d.masslimitcurve.csv', 'csv')
    TChiHH_1.addSource( 'expExclusionM1', 'orig/TChiHHexpected-1s.d.masslimitcurve.csv', 'csv' )
    TChiHH_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-004/CMS-SUS-20-004_Figure-aux_004-c.root"
    TChiHH_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-004/CMS-SUS-20-004_Figure-aux_004-c.png"
    

    ## adding these lines only makes for a silly interpolation 
    ## to mother masses between 800 and 1200 gev
    if False:
        TChiHH_1 = TChiHH.addMassPlane(2*[[x,0.]])
        TChiHH_1.addSource( 'efficiencyMap', f'orig/Eff,TChiHH-G_bin{i+1}.csv', 'csv')
        TChiHH_1.efficiencyMap._unit = ".3364"
        TChiHH_1.addSource( 'obsExclusion', 'orig/TChiHH-Gobservedlimitcurve.csv', 'csv')
        TChiHH_1.addSource( 'expExclusion', 'orig/TChiHH-Gexpectedlimitcurve.csv', 'csv')
        TChiHH_1.addSource( 'expExclusionP1', 'orig/Xseclimits,TChiHH-G_expP1.csv', 'csv')
        TChiHH_1.addSource( 'expExclusionM1', 'orig/Xseclimits,TChiHH-G_expM1.csv', 'csv')
        TChiHH_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-004/CMS-SUS-20-004_Figure-aux_003-a.root"
        TChiHH_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-004/CMS-SUS-20-004_Figure-aux_003-a.png"
        TChiHH_1.xrange = "[220,520]"
        TChiHH_1.yrange = "[2,220]"
        
        TChiHH_2 = TChiHH.addMassPlane(2*[[x,1.]])
        TChiHH_2.addSource( 'efficiencyMap', f'orig/Eff,TChiHH-G_bin{i+1}.csv', 'csv')
        TChiHH_2.efficiencyMap._unit = ".3364"
        TChiHH_2.addSource( 'obsExclusion', 'orig/TChiHH-Gobservedlimitcurve.csv', 'csv')
        TChiHH_2.addSource( 'expExclusion', 'orig/TChiHH-Gexpectedlimitcurve.csv', 'csv')
        TChiHH_2.addSource( 'expExclusionP1', 'orig/Xseclimits,TChiHH-G_expP1.csv', 'csv')
        TChiHH_2.addSource( 'expExclusionM1', 'orig/Xseclimits,TChiHH-G_expM1.csv', 'csv')
        TChiHH_2.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-004/CMS-SUS-20-004_Figure-aux_003-a.root"
        TChiHH_2.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-004/CMS-SUS-20-004_Figure-aux_003-a.png"
        TChiHH_2.xrange = "[220,520]"
        TChiHH_2.yrange = "[2,220]"

    T5HH_1 = T5HH.addMassPlane(2*[[x,x-50,1.]])
    T5HH_1.addSource( 'efficiencyMap', f'orig/Eff,T5HH_bin{i+1}.csv', 'csv')
    T5HH_1.addSource( 'obsExclusion', 'orig/T5HHobservedlimitcurve.csv', 'csv')
    T5HH_1.addSource( 'expExclusion', 'orig/T5HHexpectedlimitcurve.csv', 'csv')
    T5HH_1.addSource( 'expExclusionP1', 'orig/Xseclimits,T5HH_expP1.csv', 'csv')
    T5HH_1.addSource( 'expExclusionM1', 'orig/Xseclimits,T5HH_expM1.csv', 'csv')
    T5HH_1.dataUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-004/CMS-SUS-20-004_Figure-aux_003-b.root"
    T5HH_1.figureUrl = "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-20-004/CMS-SUS-20-004_Figure-aux_003-b.png"
    T5HH_1.xrange = "[2000,2500]"

    datasets.append ( dataset )

info.createCovarianceMatrix ( "orig/Covariancematrix.csv", "total_covar", matrixIsCorrelations = False, addOrder = True, max_datasets = None, datasets = datasets, aggregate = aggregates )
# info.createCovarianceMatrix ( "orig/Correlationmatrix.csv", "total_covar", matrixIsCorrelations = True, addOrder = True, max_datasets = None, datasets = datasets, aggregate = aggregates )

databaseCreator.create() 
