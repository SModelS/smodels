#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import copy
import csv
from smodels.tools.physicsUnits import fb
from smodels_utils.dataPreparation.argParser import getParserArgs
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z

getParserArgs()

#+++++++ global info block ++++++++++++++
info 			= MetaInfoInput('ATLAS-SUSY-2019-09')
info.url 		= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/'
info.sqrts 		= 13
info.lumi 		= 139
info.prettyName = '3 leptons EW-ino'
info.publication = "https://link.springer.com/content/pdf/10.1140/epjc/s10052-021-09749-7.pdf"
info.private 	= False
info.arxiv 		= 'https://arxiv.org/abs/2106.01676'
info.contact    = 'atlas-phys-susy-conveners@cern.ch'
info.jsonFiles  = "{'bkg_onshell_simplified.json' : ['SRWZ_1', 'SRWZ_2', 'SRWZ_3', 'SRWZ_4', 'SRWZ_5', 'SRWZ_6', 'SRWZ_7', 'SRWZ_8', 'SRWZ_9', 'SRWZ_10', 'SRWZ_11', 'SRWZ_12', 'SRWZ_13', 'SRWZ_14', 'SRWZ_15', 'SRWZ_16', 'SRWZ_17', 'SRWZ_18', 'SRWZ_19', 'SRWZ_20'], 'bkg_offshell_simplified.json' : ['SRhigh_0Jb', 'SRhigh_0Jc', 'SRhigh_0Jd', 'SRhigh_0Je', 'SRhigh_0Jf1', 'SRhigh_0Jf2', 'SRhigh_0Jg1', 'SRhigh_0Jg2', 'SRhigh_nJa', 'SRhigh_nJb', 'SRhigh_nJc', 'SRhigh_nJd', 'SRhigh_nJe', 'SRhigh_nJf', 'SRhigh_nJg', 'SRlow_0Jb', 'SRlow_0Jc', 'SRlow_0Jd', 'SRlow_0Je', 'SRlow_0Jf1', 'SRlow_0Jf2', 'SRlow_0Jg1', 'SRlow_0Jg2', 'SRlow_nJb', 'SRlow_nJc', 'SRlow_nJd', 'SRlow_nJe', 'SRlow_nJf1', 'SRlow_nJf2', 'SRlow_nJg1', 'SRlow_nJg2']}"



############# Retrieve tables content from a csv file exrracted from the tex file foud on the Arxiv paper page ################
def replace (element):
    element = element.replace('\\','')
    element = element.replace('WZi{','WZ_')
    element = element.replace('i{','')
    element = element.replace('{','')
    element = element.replace('}','')
    element = element.replace('z','_0')
    element = element.replace('n','_n')
    element = element.replace('j','J')
    element = element.replace('tcc','')
    element = element.replace('$','')
    element = element.replace('~','')
    element = element.replace('SF','_SF')
    element = element.replace('DF','_DF')
    return element

######## tables = {topology_1 : {SR_1 : {observed : x, expected : y, bgError : z}, SR_2 : {...}, ... }, topology_2 : {...}, ... } ########
tables={}
with open("orig/ANA-SUSY-2019-09-PAPER.csv",newline='') as csvfile :
    reader = csv.reader(csvfile, delimiter=" ") #Every space is a seperator between two elements of the row
    for row in reader : #We go rw by row
        if row != [] : #If the row is not empty

            #The following can be used to retrieve the name in the label of the tex table and stored in the variable topo
            # if row[0][0:7] == '\\label{' :
            #     if row[0][11:18] == 'results':
            #         topo = row[0][19:-1]
            #     else :
            #         topo = row[0][13:-1]


            if row[0] in ['{Regions}','{Regions', 'Region'] : #If the row contains the SR names
                list_SR = [] #We store the SR order in order to write the following yields in the corresponding dictionary
                for element in row :
                    if element != '' and '\SR' in element : #If the element of the row really is a SR name
                        element = replace(element) #We reshape the SR name to match the name in the corresponding .csv file title
                        topo = element[0:element.index('_')] #The name 'topo' is used to easily loop over a given topology (not mandatory). We only keep the topology part of the name
                        list_SR.append(element) #Store the SR in the list of SR names, at the following of the others
                        if topo not in tables : #If it is the first time we encounter this topology
                            tables[topo] = {element : {'observed' : None, 'expected' : None, 'bgError' : None}}
                        else : #If the dictionary for this topoogy already exists
                            tables[topo][element] = {'observed' : None, 'expected' : None, 'bgError' : None}

            if row[0] == 'Observed' : #If the row corresponds to the 'observed' yields
                index_SR = 0 #We want to write the yields in the order given by the previous 'Region' row
                # row = clean_row(row)
                for element in row :
                    if '$' in element : #If the element is a number
                        tables[topo][list_SR[index_SR]]['observed'] = float(replace(element)) #Add the number to the dictionary following the order given by list_SR (the order of the SR obtained according to the previous 'Region' row)
                        index_SR += 1  #If we just wrote a yield, the next one will be for the next SR according to the order of the previous 'Region' row

            if row[0] == 'Fitted': #If the row corresponds to the 'expected' and 'bgError' yields
                index_SR = 0 #We want to write the yields in the order given by the previous 'Region' row
                bgError = False #In the row, the 'expected' yield always comes before the corresponding 'bgError'. The boolean is to know if the element of the 'for' loop is the 'expected' or not
                # row = clean_row(row)
                for element in row :
                    if '$' in element : #If the element is a number
                        if bgError : #If the 'expected' yield has already been written, the next number is the corresponding 'bgError'
                            tables[topo][list_SR[index_SR]]['bgError'] = float(replace(element)) #Same as for the 'Observed' row
                            bgError = False #The boolean goes to False because the next number in the row must be an 'expected' yield
                            index_SR += 1 #The next number will not be for the same SR
                        else : #If it is the first number or the previous one was a 'bgError', this one must now be an 'expected' yield
                            tables[topo][list_SR[index_SR]]['expected'] = float(replace(element)) #Same as for the 'Observed' row
                            bgError = True #The next numer will be a 'bgError' but for the same SR




for topo in tables : #To go through the results
    for SR in tables[topo] :
        # print ('{' + topo + ' : {' + SR + ' : ' + str(tables[topo][SR]) + '}}') #Print the content returned if you want to check
        for event in tables[topo][SR] :
            if tables[topo][SR][event] == None : #Check that every entry of the dictionary has been rewritten
                print ("ERROR : Not all None rewritten")



        ######## TChiWZ ########
        if topo == 'SRWZ' :

            #+++++++ dataset block ++++++++++++++
            locals()[SR] = DataSetInput(SR) #locals() converts its string argument inside [] into a local variable
            locals()[SR].setInfo(dataType = 'efficiencyMap', dataId = SR,
                            observedN = tables[topo][SR]['observed'], expectedBG = tables[topo][SR]['expected'], bgError = tables[topo][SR]['bgError'])


            #+++++++ next txName block ++++++++++++++
            TChiWZ                      = locals()[SR].addTxName('TChiWZ')
            TChiWZ.checked 				= 'no'
            TChiWZ.constraint           = "[[['W']],[['Z']]]"
            TChiWZ.condition            = None
            TChiWZ.massConstraint       = [['dm > 70.0'], ['dm > 80.0']]
            TChiWZ.conditionDescription = None
            TChiWZ.source               = "ATLAS"

            #+++++++ next mass plane block ++++++++++++++
            TChiWZ_1                  = TChiWZ.addMassPlane(2*[[x, y]])
            TChiWZ_1.figure           = 'Extracted from the patchsets'
            TChiWZ_1.figureUrl        = 'n/a'
            TChiWZ_1.exclusionDataUrl = 'https://www.hepdata.net/record/ins1866951?table=Fig%2016a%20WZ%20Exclusion:%20Wino-bino(%2b),%20Obs'
            TChiWZ_1.dataUrl          = 'https://doi.org/10.17182/hepdata.95751.v1/r3'
            TChiWZ_1.setSources(dataLabels = [ 'obsExclusion', 'expExclusion'],
                                dataFiles = [ "orig/Fig16aWZExclusion:Wino-bino(+),onshell_Obs.csv",
                                             # "orig/Fig16aWZExclusion:Wino-bino(+),Obs_Up.csv",
                                             # "orig/Fig16aWZExclusion:Wino-bino(+),Obs_Down.csv",
                                              "orig/Fig16aWZExclusion:Wino-bino(+),onshell_Exp.csv",
                                            #  "orig/Fig16aWZExclusion:Wino-bino(+),Exp_Up.csv",
                                            #  "orig/Fig16aWZExclusion:Wino-bino(+),Exp_Down.csv"
                                            ],
                                units = [ None ]*2,
                                dataFormats = [ 'csv' ]*2 )
            TChiWZ_1.addSource("efficiencyMap", "orig/" + SR + "_efficiency.csv",
                              unit = "", dataFormat = "csv")




        ######## TChiWZoff ########
        if topo == 'SRlow' or topo == 'SRhigh' :

            #+++++++ dataset block ++++++++++++++
            locals()[SR] = DataSetInput(SR) #locals() converts its string argument inside [] into a local variable
            locals()[SR].setInfo(dataType = 'efficiencyMap', dataId = SR,
                            observedN = tables[topo][SR]['observed'], expectedBG = tables[topo][SR]['expected'], bgError = tables[topo][SR]['bgError'])


            #+++++++ next mass plane block ++++++++++++++
            TChiWZoff                      = locals()[SR].addTxName('TChiWZoff')
            TChiWZoff.checked              = 'no'
            TChiWZoff.constraint           = "[[['mu+','mu-']],[['l','nu']]] + [[['e+','e-']],[['l','nu']]]"
            TChiWZoff.condition            = None
            TChiWZoff.massConstraint       = [['dm < 250.'], ['dm < 250.']]
            TChiWZoff.conditionDescription = None
            TChiWZoff.source               = "ATLAS"

            #+++++++ next mass plane block ++++++++++++++ Wino/Bino +
            TChiWZoff_plus                  = TChiWZoff.addMassPlane(2*[[x, x-y]])
            TChiWZoff_plus.figure           = 'Extracted from the patchsets'
            TChiWZoff_plus.figureUrl        = 'n/a'
            TChiWZoff_plus.exclusionDataUrl = 'https://www.hepdata.net/record/ins1866951?table=Fig%2016b%20WZ%20Exclusion:%20Wino-bino(%2b)%20($\Delta%20m$),%20offshell_Obs'
            TChiWZoff_plus.dataUrl          = 'https://doi.org/10.17182/hepdata.95751.v1/r3'
            TChiWZoff_plus.setSources(dataLabels = [ 'obsExclusion', 'expExclusion'],
                                dataFiles = [ "orig/Fig16bWZExclusion:Wino-bino(+)(Deltam),offshell_Obs.csv",
                                              "orig/Fig16bWZExclusion:Wino-bino(+)(Deltam),offshell_Exp.csv"],
                                units = [ None ]*2,
                                dataFormats = [ 'csv' ]*2 )
            TChiWZoff_plus.addSource("efficiencyMap", "orig/" + SR + "_cuts_winobino(+)_efficiency.csv",
                              unit = "", dataFormat = "csv")

            TChiWZoff_plus.efficiencyMap._unit = "71."


            #+++++++ next mass plane block ++++++++++++++ Higgsino
            TChiWZoff_hino                   = TChiWZoff.addMassPlane([[x, x-y], [(2*x-y)/2,x-y]])
            TChiWZoff_hino.figure            = 'Extracted from the patchsets'
            TChiWZoff_hino.validationTarball = "TChiWZoffCd2.tar.gz"
            TChiWZoff_hino.figureUrl         = 'n/a'
            TChiWZoff_hino.exclusionDataUrl  = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/fig_16d.png'
            TChiWZoff_hino.dataUrl           = 'https://doi.org/10.17182/hepdata.95751.v1/r3'
            TChiWZoff_hino.setSources(dataLabels = [ 'obsExclusion', 'expExclusion'],
                                    ############ Use global (offshell+compressed) or onshell only exclusion curves ? ##############
                                dataFiles = [ "orig/Fig16dWZExclusion:Higgsino(Deltam),offshell_Obs.csv",
                                              # "orig/Fig16dWZExclusion:Higgsino(Deltam),Obs_Up.csv",
                                              # "orig/Fig16dWZExclusion:Higgsino(Deltam),Obs_Down.csv",
                                              "orig/Fig16dWZExclusion:Higgsino(Deltam),offshell_Exp.csv",
                                              # "orig/Fig16dWZExclusion:Higgsino(Deltam),Exp_Up.csv",
                                              # "orig/Fig16dWZExclusion:Higgsino(Deltam),Exp_Down.csv"
                                              ],
                                units = [ None ]*2,
                                dataFormats = [ 'csv' ]*2 )
            TChiWZoff_hino.addSource("efficiencyMap", "orig/" + SR + "_cuts_higgsino_efficiency.csv",
                              unit = "", dataFormat = "csv")

            TChiWZoff_hino.efficiencyMap._unit = "71."





####### WH - SRs not exctracted fro the patchsets #######

#+++++++ dataset block ++++++++++++++
SR_WH_0j = DataSetInput('SR_WH_0j') ####################################### No clear description of this more inclusive SR, need to count the SR with n_jet=0 and mll>105 GeV ?
                                    ####################################### Here the SR with mll>105 GeV were not taken into account, even when n_jet=0
SR_WH_0j.setInfo(dataType = 'efficiencyMap', dataId = 'SR_WH_0j',
                observedN = 261, expectedBG = 244, bgError = 14)

#+++++++ next txName block ++++++++++++++
TChiWH                      = SR_WH_0j.addTxName('TChiWH')
TChiWH.checked 				= 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.condition            = None
TChiWH.massConstraint       = None
TChiWH.conditionDescription = None
TChiWH.source               = "ATLAS"

#+++++++ next mass plane block ++++++++++++++
TChiWH_1                  = TChiWH.addMassPlane(2*[[x,y]])
TChiWH_1.figure           = 'Aux. Fig. 11a; Aux. Fig. 11b'
TChiWH_1.figureUrl        = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/figaux_11a.png; https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/figaux_11b.png'
TChiWH_1.exclusionDataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/fig_17.png'
TChiWH_1.dataUrl          = 'https://www.hepdata.net/record/ins1866951?table=AuxFig%2011a%20Acc:%20Onshell%20SR$_{low-m_{ll}-0j}^{Wh}$; https://www.hepdata.net/record/ins1866951?table=AuxFig%2011b%20Eff:%20Onshell%20SR$_{low-m_{ll}-0j}^{Wh}$'
TChiWH_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
                    dataFiles = [ "orig/Fig17WhExclusion,Obs.csv",
                                  "orig/Fig17WhExclusion,Obs_Up.csv",
                                  "orig/Fig17WhExclusion,Obs_Down.csv",
                                  "orig/Fig17WhExclusion,Exp.csv",
                                  "orig/Fig17WhExclusion,Exp_Up.csv",
                                  "orig/Fig17WhExclusion,Exp_Down.csv"],
                    units = [ None ]*6,
                    dataFormats = [ 'csv' ]*6 )
TChiWH_1.addSource("efficiencyMap", ("orig/AuxFig11aAcc:OnshellSR_{low-m_{ll}-0j}^{Wh}.csv", "orig/AuxFig11bEff:OnshellSR_{low-m_{ll}-0j}^{Wh}.csv"),
                  unit = "/100000", dataFormat = "mcsv")




#+++++++ dataset block ++++++++++++++
SR_WH_nj = DataSetInput('SR_WH_nj')
SR_WH_nj.setInfo(dataType = 'efficiencyMap', dataId = 'SR_WH_nj',
                observedN = 488, expectedBG = 492.6, bgError = 26)

#+++++++ next txName block ++++++++++++++
TChiWH                      = SR_WH_nj.addTxName('TChiWH')
TChiWH.checked 				= 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.condition            = None
TChiWH.massConstraint       = None
TChiWH.conditionDescription = None
TChiWH.source               = "ATLAS"

#+++++++ next mass plane block ++++++++++++++
TChiWH_1                  = TChiWH.addMassPlane(2*[[x,y]])
TChiWH_1.figure           = 'Aux. Fig. 11c; Aux. Fig. 11d'
TChiWH_1.figureUrl        = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/figaux_11c.png; https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/figaux_11d.png'
TChiWH_1.exclusionDataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/fig_17.png'
TChiWH_1.dataUrl = 'https://www.hepdata.net/record/ins1866951?table=AuxFig%2011c%20Acc:%20Onshell%20SR$_{low-m_{ll}-nj}^{Wh}$; https://www.hepdata.net/record/ins1866951?table=AuxFig%2011d%20Eff:%20Onshell%20SR$_{low-m_{ll}-nj}^{Wh}$'
TChiWH_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
                    dataFiles = [ "orig/Fig17WhExclusion,Obs.csv",
                                  "orig/Fig17WhExclusion,Obs_Up.csv",
                                  "orig/Fig17WhExclusion,Obs_Down.csv",
                                  "orig/Fig17WhExclusion,Exp.csv",
                                  "orig/Fig17WhExclusion,Exp_Up.csv",
                                  "orig/Fig17WhExclusion,Exp_Down.csv"],
                    units = [ None ]*6,
                    dataFormats = [ 'csv' ]*6 )
TChiWH_1.addSource("efficiencyMap", ("orig/AuxFig11cAcc:OnshellSR_{low-m_{ll}-nj}^{Wh}.csv", "orig/AuxFig11dEff:OnshellSR_{low-m_{ll}-nj}^{Wh}.csv"),
                  unit = "/100000", dataFormat = "mcsv")




#+++++++ dataset block ++++++++++++++
SR_WH_DFOS = DataSetInput('SR_WH_DFOS')
SR_WH_DFOS.setInfo(dataType = 'efficiencyMap', dataId = 'SR_WH_DFOS',
                observedN = 20, expectedBG = 11.5, bgError = 2.4)

#+++++++ next txName block ++++++++++++++
TChiWH                      = SR_WH_DFOS.addTxName('TChiWH')
TChiWH.checked 				= 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.condition            = None
TChiWH.massConstraint       = None
TChiWH.conditionDescription = None
TChiWH.source               = "ATLAS"

#+++++++ next mass plane block ++++++++++++++
TChiWH_1                  = TChiWH.addMassPlane(2*[[x,y]])
TChiWH_1.figure           = 'Aux. Fig. 11e; Aux. Fig. 11f'
TChiWH_1.figureUrl        = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/figaux_11e.png; https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/figaux_11f.png'
TChiWH_1.exclusionDataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-09/fig_17.png'
TChiWH_1.dataUrl = 'https://www.hepdata.net/record/ins1866951?table=AuxFig%2011e%20Acc:%20Onshell%20SR$_{DFOS}^{Wh}$; https://www.hepdata.net/record/ins1866951?table=AuxFig%2011f%20Eff:%20Onshell%20SR$_{DFOS}^{Wh}$'
TChiWH_1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'expExclusionP1', 'expExclusionM1'],
                    dataFiles = [ "orig/Fig17WhExclusion,Obs.csv",
                                  "orig/Fig17WhExclusion,Obs_Up.csv",
                                  "orig/Fig17WhExclusion,Obs_Down.csv",
                                  "orig/Fig17WhExclusion,Exp.csv",
                                  "orig/Fig17WhExclusion,Exp_Up.csv",
                                  "orig/Fig17WhExclusion,Exp_Down.csv"],
                    units = [ None ]*6,
                    dataFormats = [ 'csv' ]*6 )
TChiWH_1.addSource("efficiencyMap", ("orig/AuxFig11eAcc:OnshellSR_{DFOS}^{Wh}.csv", "orig/AuxFig11fEff:OnshellSR_{DFOS}^{Wh}.csv"),
                  unit = "/100000", dataFormat = "mcsv")


databaseCreator.create()
