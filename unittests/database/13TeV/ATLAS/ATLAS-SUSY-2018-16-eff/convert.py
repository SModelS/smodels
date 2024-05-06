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
from smodels.base.physicsUnits import fb
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z

datasets_ewk_excl = {'CR_ewk_VV_high': {'SRname': 'CRVV_MLL_hghmet_cuts', 'observedN': 687, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'CR_ewk_VV_low': {'SRname': 'CRVV_MLL_lowmet_cuts', 'observedN': 721, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'CR_ewk_tau_high': {'SRname': 'CRtau_MLL_hghmet_cuts', 'observedN': 592, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'CR_ewk_tau_low': {'SRname': 'CRtau_MLL_lowmet_cuts', 'observedN': 2247, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'CR_ewk_top_high': {'SRname': 'CRtop_MLL_hghmet_cuts', 'observedN': 4327, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'CR_ewk_top_low': {'SRname': 'CRtop_MLL_lowmet_cuts', 'observedN': 6164, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            # 'SR_ewk_1l1T_a': {'SRname': 'SR_eMLLa_Onelep1track_cuts', 'observedN': 0, 'expectedBG': 0.5, 'bgError': 0.5},
            # 'SR_ewk_1l1T_b': {'SRname': 'SR_eMLLb_Onelep1track_cuts', 'observedN': 8, 'expectedBG': 6.0, 'bgError': 1.9},
            # 'SR_ewk_1l1T_c': {'SRname': 'SR_eMLLc_Onelep1track_cuts', 'observedN': 8, 'expectedBG': 7.6, 'bgError': 2.1},
            # 'SR_ewk_1l1T_d': {'SRname': 'SR_eMLLd_Onelep1track_cuts', 'observedN': 24, 'expectedBG': 20.7, 'bgError': 3.4},
            # 'SR_ewk_1l1T_e': {'SRname': 'SR_eMLLe_Onelep1track_cuts', 'observedN': 24, 'expectedBG': 24., 'bgError': 4.},
            # 'SR_ewk_1l1T_f': {'SRname': 'SR_eMLLf_Onelep1track_cuts', 'observedN': 16, 'expectedBG': 18.1, 'bgError': 3.1},
            'SR_ewk_2l_ee_high_c': {'SRname': 'SRee_eMLLc_hghmet_cuts', 'observedN': 1, 'expectedBG': 0.7, 'bgError': 0.4},
            'SR_ewk_2l_ee_low_c': {'SRname': 'SRee_eMLLc_lowmet_deltaM_high_cuts', 'observedN': 7, 'expectedBG': 5.3, 'bgError': 1.5},
            'SR_ewk_2l_ee_med_c': {'SRname': 'SRee_eMLLc_lowmet_deltaM_low_cuts', 'observedN': 0, 'expectedBG': 0.11, 'bgError': 0.08},
            'SR_ewk_2l_ee_high_d': {'SRname': 'SRee_eMLLd_hghmet_cuts', 'observedN': 16, 'expectedBG': 10.3, 'bgError': 2.5},
            'SR_ewk_2l_ee_low_d': {'SRname': 'SRee_eMLLd_lowmet_deltaM_high_cuts', 'observedN': 11, 'expectedBG': 8.6, 'bgError': 1.8},
            'SR_ewk_2l_ee_med_d': {'SRname': 'SRee_eMLLd_lowmet_deltaM_low_cuts', 'observedN': 4, 'expectedBG': 5.1, 'bgError': 1.6},
            'SR_ewk_2l_ee_high_e': {'SRname': 'SRee_eMLLe_hghmet_cuts', 'observedN': 13, 'expectedBG': 12.1, 'bgError': 2.2},
            'SR_ewk_2l_ee_low_e': {'SRname': 'SRee_eMLLe_lowmet_deltaM_high_cuts', 'observedN': 16, 'expectedBG': 16.7, 'bgError': 2.5},
            'SR_ewk_2l_ee_med_e': {'SRname': 'SRee_eMLLe_lowmet_deltaM_low_cuts', 'observedN': 11, 'expectedBG': 7.3, 'bgError': 1.9},
            'SR_ewk_2l_ee_high_f': {'SRname': 'SRee_eMLLf_hghmet_cuts', 'observedN': 8, 'expectedBG': 10.1, 'bgError': 1.7},
            'SR_ewk_2l_ee_low_f': {'SRname': 'SRee_eMLLf_lowmet_deltaM_high_cuts', 'observedN': 16, 'expectedBG': 15.5, 'bgError': 2.6},
            'SR_ewk_2l_ee_med_f': {'SRname': 'SRee_eMLLf_lowmet_deltaM_low_cuts', 'observedN': 4, 'expectedBG': 2.2, 'bgError': 0.9},
            'SR_ewk_2l_ee_high_g': {'SRname': 'SRee_eMLLg_hghmet_cuts', 'observedN': 8, 'expectedBG': 10.4, 'bgError': 1.7},
            'SR_ewk_2l_ee_low_g': {'SRname': 'SRee_eMLLg_lowmet_deltaM_high_cuts', 'observedN': 10, 'expectedBG': 12.9, 'bgError': 2.1},
            'SR_ewk_2l_ee_high_h': {'SRname': 'SRee_eMLLh_hghmet_cuts', 'observedN': 18, 'expectedBG': 19.3, 'bgError': 2.5},
            'SR_ewk_2l_ee_low_h': {'SRname': 'SRee_eMLLh_lowmet_deltaM_high_cuts', 'observedN': 9, 'expectedBG': 18.8, 'bgError': 2.2},
            'SR_ewk_2l_mm_high_a': {'SRname': 'SRmm_eMLLa_hghmet_cuts', 'observedN': 5, 'expectedBG': 3.4, 'bgError': 1.2},
            'SR_ewk_2l_mm_low_a': {'SRname': 'SRmm_eMLLa_lowmet_deltaM_high_cuts', 'observedN': 9, 'expectedBG': 15.4, 'bgError': 2.4},
            'SR_ewk_2l_mm_med_a': {'SRname': 'SRmm_eMLLa_lowmet_deltaM_low_cuts', 'observedN': 16, 'expectedBG': 14.6, 'bgError': 2.9},
            'SR_ewk_2l_mm_high_b': {'SRname': 'SRmm_eMLLb_hghmet_cuts', 'observedN': 5, 'expectedBG': 3.5, 'bgError': 1.3},
            'SR_ewk_2l_mm_low_b': {'SRname': 'SRmm_eMLLb_lowmet_deltaM_high_cuts', 'observedN': 7, 'expectedBG': 8., 'bgError': 1.7},
            'SR_ewk_2l_mm_med_b': {'SRname': 'SRmm_eMLLb_lowmet_deltaM_low_cuts', 'observedN': 8, 'expectedBG': 6.9, 'bgError': 2.1},
            'SR_ewk_2l_mm_high_c': {'SRname': 'SRmm_eMLLc_hghmet_cuts', 'observedN': 0, 'expectedBG': 3.9, 'bgError': 1.3},
            'SR_ewk_2l_mm_low_c': {'SRname': 'SRmm_eMLLc_lowmet_deltaM_high_cuts', 'observedN': 7, 'expectedBG': 6.5, 'bgError': 1.6},
            'SR_ewk_2l_mm_med_c': {'SRname': 'SRmm_eMLLc_lowmet_deltaM_low_cuts', 'observedN': 6, 'expectedBG': 6.2, 'bgError': 1.9},
            'SR_ewk_2l_mm_high_d': {'SRname': 'SRmm_eMLLd_hghmet_cuts', 'observedN': 9, 'expectedBG': 11., 'bgError': 2.},
            'SR_ewk_2l_mm_low_d': {'SRname': 'SRmm_eMLLd_lowmet_deltaM_high_cuts', 'observedN': 12, 'expectedBG': 11.3, 'bgError': 1.9},
            'SR_ewk_2l_mm_med_d': {'SRname': 'SRmm_eMLLd_lowmet_deltaM_low_cuts', 'observedN': 41, 'expectedBG': 34, 'bgError': 4},
            'SR_ewk_2l_mm_high_e': {'SRname': 'SRmm_eMLLe_hghmet_cuts', 'observedN': 23, 'expectedBG': 17.8, 'bgError': 2.7},
            'SR_ewk_2l_mm_low_e': {'SRname': 'SRmm_eMLLe_lowmet_deltaM_high_cuts', 'observedN': 17, 'expectedBG': 15.6, 'bgError': 2.3},
            'SR_ewk_2l_mm_med_e': {'SRname': 'SRmm_eMLLe_lowmet_deltaM_low_cuts', 'observedN': 59, 'expectedBG': 52, 'bgError': 6},
            'SR_ewk_2l_mm_high_f': {'SRname': 'SRmm_eMLLf_hghmet_cuts', 'observedN': 3, 'expectedBG': 8.3, 'bgError': 1.4},
            'SR_ewk_2l_mm_low_f': {'SRname': 'SRmm_eMLLf_lowmet_deltaM_high_cuts', 'observedN': 18, 'expectedBG': 16.7, 'bgError': 2.3},
            'SR_ewk_2l_mm_med_f': {'SRname': 'SRmm_eMLLf_lowmet_deltaM_low_cuts', 'observedN': 21, 'expectedBG': 18.5, 'bgError': 3.2},
            'SR_ewk_2l_mm_high_g': {'SRname': 'SRmm_eMLLg_hghmet_cuts', 'observedN': 5, 'expectedBG': 10.1, 'bgError': 1.5},
            'SR_ewk_2l_mm_low_g': {'SRname': 'SRmm_eMLLg_lowmet_deltaM_high_cuts', 'observedN': 16, 'expectedBG': 15.3, 'bgError': 2.},
            'SR_ewk_2l_mm_high_h': {'SRname': 'SRmm_eMLLh_hghmet_cuts', 'observedN': 20, 'expectedBG': 19.6, 'bgError': 2.3},
            'SR_ewk_2l_mm_low_h': {'SRname': 'SRmm_eMLLh_lowmet_deltaM_high_cuts', 'observedN': 44, 'expectedBG': 35.9, 'bgError': 3.3},
            }

datasets_ewk_incl = {'SR_ewk_1l1T': {'Topologies': ['TChiZoff','TChiWZoff'], 'letter': ['g','h'], 'observedN': 80, 'expectedBG': 76.9, 'bgError': 6.7},
                    'SR_ewk_2l_low': {'Topologies': ['TChiZoff','TChiWZoff','TChiWWoff'], 'letter': ['c','d'], 'observedN': 199, 'expectedBG': 202.5, 'bgError': 8.2},
                    'SR_ewk_2l_med': {'Topologies': ['TChiZoff','TChiWZoff','TChiWWoff'], 'letter': ['e','f'], 'observedN': 170, 'expectedBG': 146.9, 'bgError': 9.3},
                    'SR_ewk_2l_high': {'Topologies': ['TChiZoff','TChiWZoff','TChiWWoff'], 'letter': ['a','b'], 'observedN': 134, 'expectedBG': 140.5, 'bgError': 7.0}
                    }

datasets_slep = {'CR_slep_VV_high': {'SRname': 'CRVV_MT2_hghmet_cuts', 'observedN': 367, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'CR_slep_VV_low': {'SRname': 'CRVV_MT2_lowmet_cuts', 'observedN': 272, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'CR_slep_tau_high': {'SRname': 'CRtau_MT2_hghmet_cuts', 'observedN': 218, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'CR_slep_tau_low': {'SRname': 'CRtau_MT2_lowmet_cuts', 'observedN': 799, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'CR_slep_top_high': {'SRname': 'CRtop_MT2_hghmet_cuts', 'observedN': 2053, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'CR_slep_top_low': {'SRname': 'CRtop_MT2_lowmet_cuts', 'observedN': 3106, 'expectedBG': 1, 'bgError': 1, 'upperLimit': 10000000*fb, 'expectedUpperLimit': 10000000*fb},
            'SR_slep_2l_ee_high_a': {'SRname': 'SRee_eMT2a_hghmet_cuts', 'observedN': 3, 'expectedBG': 4.0, 'bgError': 1.1},
            'SR_slep_2l_ee_low_a': {'SRname': 'SRee_eMT2a_lowmet_V2_cuts', 'observedN': 8, 'expectedBG': 6.0, 'bgError': 1.4},
            'SR_slep_2l_ee_high_b': {'SRname': 'SRee_eMT2b_hghmet_cuts', 'observedN': 3, 'expectedBG': 3.6, 'bgError': 1.0},
            'SR_slep_2l_ee_low_b': {'SRname': 'SRee_eMT2b_lowmet_V2_cuts', 'observedN': 5, 'expectedBG': 5.3, 'bgError': 2.1},
            'SR_slep_2l_ee_high_c': {'SRname': 'SRee_eMT2c_hghmet_cuts', 'observedN': 9, 'expectedBG': 7.9, 'bgError': 1.9},
            'SR_slep_2l_ee_low_c': {'SRname': 'SRee_eMT2c_lowmet_V2_cuts', 'observedN': 15, 'expectedBG': 11.6, 'bgError': 2.5},
            'SR_slep_2l_ee_high_d': {'SRname': 'SRee_eMT2d_hghmet_cuts', 'observedN': 13, 'expectedBG': 13.2, 'bgError': 2.1},
            'SR_slep_2l_ee_low_d': {'SRname': 'SRee_eMT2d_lowmet_V2_cuts', 'observedN': 19, 'expectedBG': 22.9, 'bgError': 3.3},
            'SR_slep_2l_ee_high_e': {'SRname': 'SRee_eMT2e_hghmet_cuts', 'observedN': 9, 'expectedBG': 8.6, 'bgError': 1.4},
            'SR_slep_2l_ee_low_e': {'SRname': 'SRee_eMT2e_lowmet_V2_cuts', 'observedN': 30, 'expectedBG': 31., 'bgError': 4.},
            'SR_slep_2l_ee_high_f': {'SRname': 'SRee_eMT2f_hghmet_cuts', 'observedN': 6, 'expectedBG': 5.7, 'bgError': 1.0},
            'SR_slep_2l_ee_low_f': {'SRname': 'SRee_eMT2f_lowmet_V2_cuts', 'observedN': 24, 'expectedBG': 23.3, 'bgError': 3.0},
            'SR_slep_2l_ee_high_g': {'SRname': 'SRee_eMT2g_hghmet_cuts', 'observedN': 8, 'expectedBG': 7.0, 'bgError': 1.2},
            'SR_slep_2l_ee_low_g': {'SRname': 'SRee_eMT2g_lowmet_V2_cuts', 'observedN': 32, 'expectedBG': 27.1, 'bgError': 3.1},
            'SR_slep_2l_ee_high_h': {'SRname': 'SRee_eMT2h_hghmet_cuts', 'observedN': 6, 'expectedBG': 6.8, 'bgError': 1.1},
            'SR_slep_2l_ee_low_h': {'SRname': 'SRee_eMT2h_lowmet_V2_cuts', 'observedN': 11, 'expectedBG': 16.8, 'bgError': 2.1},
            'SR_slep_2l_mm_high_a': {'SRname': 'SRmm_eMT2a_hghmet_cuts', 'observedN': 10, 'expectedBG': 11.0, 'bgError': 2.2},
            'SR_slep_2l_mm_low_a': {'SRname': 'SRmm_eMT2a_lowmet_V2_cuts', 'observedN': 3, 'expectedBG': 5.2, 'bgError': 1.1},
            'SR_slep_2l_mm_high_b': {'SRname': 'SRmm_eMT2b_hghmet_cuts', 'observedN': 3, 'expectedBG': 5.8, 'bgError': 1.3},
            'SR_slep_2l_mm_low_b': {'SRname': 'SRmm_eMT2b_lowmet_V2_cuts', 'observedN': 6, 'expectedBG': 4.3, 'bgError': 1.0},
            'SR_slep_2l_mm_high_c': {'SRname': 'SRmm_eMT2c_hghmet_cuts', 'observedN': 11, 'expectedBG': 8.6, 'bgError': 1.6},
            'SR_slep_2l_mm_low_c': {'SRname': 'SRmm_eMT2c_lowmet_V2_cuts', 'observedN': 15, 'expectedBG': 12.8, 'bgError': 1.8},
            'SR_slep_2l_mm_high_d': {'SRname': 'SRmm_eMT2d_hghmet_cuts', 'observedN': 12, 'expectedBG': 14.2, 'bgError': 1.9},
            'SR_slep_2l_mm_low_d': {'SRname': 'SRmm_eMT2d_lowmet_V2_cuts', 'observedN': 23, 'expectedBG': 24.8, 'bgError': 2.6},
            'SR_slep_2l_mm_high_e': {'SRname': 'SRmm_eMT2e_hghmet_cuts', 'observedN': 9, 'expectedBG': 10.0, 'bgError': 1.5},
            'SR_slep_2l_mm_low_e': {'SRname': 'SRmm_eMT2e_lowmet_V2_cuts', 'observedN': 37, 'expectedBG': 38., 'bgError': 5.},
            'SR_slep_2l_mm_high_f': {'SRname': 'SRmm_eMT2f_hghmet_cuts', 'observedN': 11, 'expectedBG': 11.2, 'bgError': 1.6},
            'SR_slep_2l_mm_low_f': {'SRname': 'SRmm_eMT2f_lowmet_V2_cuts', 'observedN': 44, 'expectedBG': 37.8, 'bgError': 3.3},
            'SR_slep_2l_mm_high_g': {'SRname': 'SRmm_eMT2g_hghmet_cuts', 'observedN': 10, 'expectedBG': 11.5, 'bgError': 1.5},
            'SR_slep_2l_mm_low_g': {'SRname': 'SRmm_eMT2g_lowmet_V2_cuts', 'observedN': 41, 'expectedBG': 36.0, 'bgError': 3.4},
            'SR_slep_2l_mm_high_h': {'SRname': 'SRmm_eMT2h_hghmet_cuts', 'observedN': 8, 'expectedBG': 7.8, 'bgError': 1.4},
            'SR_slep_2l_mm_low_h': {'SRname': 'SRmm_eMT2h_lowmet_V2_cuts', 'observedN': 28, 'expectedBG': 28.0, 'bgError': 2.7}
            }

#+++++++ global info block ++++++++++++++
info 			= MetaInfoInput('ATLAS-SUSY-2018-16')
info.url 		= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-16/'
info.sqrts 		= 13
info.lumi 		= 139
info.prettyName = '2 soft l + jets, EWK'
info.publication = 'Phys. Rev. D 101 (2020) 052005'
info.publicationDOI = 'http://doi.org/10.1103/PhysRevD.101.052005'
info.private 	= False
info.arxiv 		= 'http://arxiv.org/abs/1911.12606'
info.contact    = 'atlas-phys-susy-conveners@cern.ch'
info.comment    = 'Be careful when testing models allowing for the production of only electrons or muons against this analysis. The EMs in the slepton CRs are correct only if both leptons can be produced.\nMoreover, but not related, SR_ewk_1l1T, SR_ewk_2l_low, SR_ewk_2l_med and SR_ewk_2l_high are inclusive SRs with different mass plane than the exclusive ones. Fake individual json files were created in order for SModelS to get a result from them even when the combination of SRs is on.'
info.jsonFiles  = "{'orig/EWKinos_bkgonly_SR2lOnly.json': " + str(list(datasets_ewk_excl.keys())) + ", 'orig/Sleptons_bkgonly.json': " + str(list(datasets_slep.keys())) + ", 'orig/SR_ewk_1l1T_bkgonly.json': ['SR_ewk_1l1T'], 'orig/SR_ewk_2l_low_bkgonly.json': ['SR_ewk_2l_low'], 'orig/SR_ewk_2l_med_bkgonly.json': ['SR_ewk_2l_med'], 'orig/SR_ewk_2l_high_bkgonly.json': ['SR_ewk_2l_high']" + "}"
info.includeCRs = True
info.signalUncertainty = 0.22

for SR in datasets_ewk_excl:
    #+++++++ dataset block ++++++++++++++
    dataset = DataSetInput(SR)
    if 'upperLimit' in datasets_ewk_excl[SR].keys():
        dataset.setInfo(dataType = 'efficiencyMap', dataId = SR,
                        observedN = datasets_ewk_excl[SR]['observedN'], expectedBG = datasets_ewk_excl[SR]['expectedBG'], bgError = datasets_ewk_excl[SR]['bgError'], upperLimit = datasets_ewk_excl[SR]['upperLimit'], expectedUpperLimit = datasets_ewk_excl[SR]['expectedUpperLimit'])
    else:
        dataset.setInfo(dataType = 'efficiencyMap', dataId = SR,
                        observedN = datasets_ewk_excl[SR]['observedN'], expectedBG = datasets_ewk_excl[SR]['expectedBG'], bgError = datasets_ewk_excl[SR]['bgError'])

    #+++++++ next mass plane block ++++++++++++++
    TChiWZoff                      = dataset.addTxName('TChiWZoff')
    TChiWZoff.checked              = 'no'
    if "2l_ee" in SR:
        TChiWZoff.constraint       = "[[['e+','e-']],[['jet','jet']]]"
    elif "2l_mm" in SR:
        TChiWZoff.constraint       = "[[['mu+','mu-']],[['jet','jet']]]"
    # elif "1l1T" in SR:
    #     TChiWZoff.constraint       = "[[['mu+','mu-']],[['jet','jet']]] + [[['e+','e-']],[['jet','jet']]]"
    elif "CR" in SR:
        TChiWZoff.constraint       = "[[['mu+','mu-']],[['jet','jet']]] + [[['e+','e-']],[['jet','jet']]]"
    TChiWZoff.condition            = None
    TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
    TChiWZoff.conditionDescription = None
    TChiWZoff.source               = "ATLAS"


    #+++++++ next mass plane block ++++++++++++++
    TChiWZoff_plus                   = TChiWZoff.addMassPlane([[x, x-y], [x, x-y]])
    TChiWZoff_plus.figure            = 'Extracted from the patchsets'
    TChiWZoff_plus.figureUrl         = 'n/a'
    TChiWZoff_plus.exclusionDataUrl  = 'https://doi.org/10.17182/hepdata.91374.v5/t5; https://doi.org/10.17182/hepdata.91374.v5/t6'
    TChiWZoff_plus.dataUrl           = 'n/a'
    TChiWZoff_plus.validationTarball = "TChiWZoff_2018-16.tar.gz"
    TChiWZoff_plus.setSources(dataLabels = [ 'obsExclusion', 'expExclusion'],
                        dataFiles = [ "orig/HEPData-ins1767649-v5-Figure_14c_Observed.csv",
                                      "orig/HEPData-ins1767649-v5-Figure_14c_Expected.csv",
                                      ],
                        units = [ None ]*2,
                        dataFormats = [ 'csv' ]*2 )
    TChiWZoff_plus.addSource("efficiencyMap", "orig/interpolated_" + datasets_ewk_excl[SR]['SRname'] + "_winobino(+)_efficiency.csv",
                      unit = "", dataFormat = "csv")

for SR in datasets_ewk_incl:
    #+++++++ dataset block ++++++++++++++
    dataset = DataSetInput(SR)
    if 'upperLimit' in datasets_ewk_incl[SR].keys():
        dataset.setInfo(dataType = 'efficiencyMap', dataId = SR,
                        observedN = datasets_ewk_incl[SR]['observedN'], expectedBG = datasets_ewk_incl[SR]['expectedBG'], bgError = datasets_ewk_incl[SR]['bgError'], upperLimit = datasets_ewk_incl[SR]['upperLimit'], expectedUpperLimit = datasets_ewk_incl[SR]['expectedUpperLimit'])
    else:
        dataset.setInfo(dataType = 'efficiencyMap', dataId = SR,
                        observedN = datasets_ewk_incl[SR]['observedN'], expectedBG = datasets_ewk_incl[SR]['expectedBG'], bgError = datasets_ewk_incl[SR]['bgError'])

    if 'TChiZoff' in datasets_ewk_incl[SR]['Topologies']:
        #+++++++ next mass plane block ++++++++++++++
        TChiZoff                      = dataset.addTxName('TChiZoff')
        TChiZoff.checked              = 'no'
        TChiZoff.constraint           = "[[['e+','e-']],[]] + [[['mu+','mu-']],[]]"
        TChiZoff.condition            = None
        TChiZoff.massConstraint       = [['dm < 86.0'],[]]
        TChiZoff.conditionDescription = None
        TChiZoff.source               = "ATLAS"


        #+++++++ next mass plane block ++++++++++++++
        TChiZoff_hino                   = TChiZoff.addMassPlane([[x, x-y], [x-y]])
        TChiZoff_hino.figure            = 'Aux. Fig. 29 (a-h)'
        TChiZoff_hino.figureUrl         = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-16/figaux_29a.png'
        TChiZoff_hino.exclusionDataUrl  = 'https://doi.org/10.17182/hepdata.91374.v5/t1; https://doi.org/10.17182/hepdata.91374.v5/t2'
        TChiZoff_hino.dataUrl           = 'https://doi.org/10.17182/hepdata.91374.v5/t33 (t33-t40)'
        # TChiZoff_hino.validationTarball =  "TChi_nonDegenHinoLSP_BRN2toZN1100.tar.gz"
        TChiZoff_hino.setSources(dataLabels = [ 'obsExclusion', 'expExclusion'],
                            dataFiles = [ "orig/HEPData-ins1767649-v5-Figure_14a_Observed.csv",
                                          "orig/HEPData-ins1767649-v5-Figure_14a_Expected.csv",
                                          ],
                            units = [ None ]*2,
                            dataFormats = [ 'csv' ]*2 )
        TChiZoff_hino.addSource("efficiencyMap", (f"orig/interpolated_Figure29{datasets_ewk_incl[SR]['letter'][0]}.csv", f"orig/interpolated_Figure29{datasets_ewk_incl[SR]['letter'][1]}.csv"),
                          unit = "", dataFormat = "mcsv")
        TChiZoff_hino.efficiencyMap._unit = "/672.9" # Acceptance unit is /10000 and BR(Z*->e+e-) + BR(Z*->µ+µ-) = 0.06729

    if 'TChiWZoff' in datasets_ewk_incl[SR]['Topologies']:
        #+++++++ next mass plane block ++++++++++++++
        TChiWZoff                      = dataset.addTxName('TChiWZoff')
        TChiWZoff.checked              = 'no'
        TChiWZoff.constraint           = "[[['e+','e-']],[['jet','jet']]] + [[['mu+','mu-']],[['jet','jet']]]"
        TChiWZoff.condition            = None
        TChiWZoff.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
        TChiWZoff.conditionDescription = None
        TChiWZoff.source               = "ATLAS"


        #+++++++ next mass plane block ++++++++++++++
        TChiWZoff_hino                   = TChiWZoff.addMassPlane([[x, x-y], [x-y/2, x-y]])
        TChiWZoff_hino.figure            = 'Aux. Fig. 30 (a-h); Aux. Fig. 31 (a-h)'
        TChiWZoff_hino.figureUrl         = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-16/figaux_30a.png; https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-16/figaux_31a.png'
        TChiWZoff_hino.exclusionDataUrl  = 'https://doi.org/10.17182/hepdata.91374.v5/t1; https://doi.org/10.17182/hepdata.91374.v5/t2'
        TChiWZoff_hino.dataUrl           = 'https://doi.org/10.17182/hepdata.91374.v5/t41 (t41-t56)'
        TChiWZoff_hino.validationTarball =  "TChi_nonDegenHinoLSP_BRN2toZN1100.tar.gz"
        # TChiWZoff_hino.validationTarball =  "TChiWZoff_2018-16.tar.gz"
        TChiWZoff_hino.setSources(dataLabels = [ 'obsExclusion', 'expExclusion'],
                            dataFiles = [ "orig/HEPData-ins1767649-v5-Figure_14a_Observed.csv",
                                          "orig/HEPData-ins1767649-v5-Figure_14a_Expected.csv",
                                          ],
                            units = [ None ]*2,
                            dataFormats = [ 'csv' ]*2 )
        TChiWZoff_hino.addSource("efficiencyMap", (f"orig/interpolated_meanFigures3031{datasets_ewk_incl[SR]['letter'][0]}.csv", f"orig/interpolated_meanFigures3031{datasets_ewk_incl[SR]['letter'][1]}.csv"),
                          unit = "", dataFormat = "mcsv")
        TChiWZoff_hino.efficiencyMap._unit = "/450.8" # Acceptance unit is /10000 and (BR(Z*->e+e-) + BR(Z*->µ+µ-))*(BR(W*->u db) + BR(W*->c sb)) = 0.06729*0.67000 =  0.04508

    if 'TChiWWoff' in datasets_ewk_incl[SR]['Topologies']:
        #+++++++ next mass plane block ++++++++++++++
        TChiWWoff                      = dataset.addTxName('TChiWWoff')
        TChiWWoff.checked              = 'no'
        TChiWWoff.constraint           = "[[['e+','nu']],[['e-','nu']]] + [[['mu+','nu']],[['mu-','nu']]]"
        TChiWWoff.condition            = None
        TChiWWoff.massConstraint       = [['dm < 76.0']]*2
        TChiWWoff.conditionDescription = None
        TChiWWoff.source               = "ATLAS"


        #+++++++ next mass plane block ++++++++++++++
        TChiWWoff_hino                   = TChiWWoff.addMassPlane([[x-y/2, x-y], [x-y/2, x-y]])
        TChiWWoff_hino.figure            = 'Aux. Fig. 28 (a-f)'
        TChiWWoff_hino.figureUrl         = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-16/figaux_28a.png'
        TChiWWoff_hino.exclusionDataUrl  = 'https://doi.org/10.17182/hepdata.91374.v5/t1; https://doi.org/10.17182/hepdata.91374.v5/t2'
        TChiWWoff_hino.dataUrl           = 'https://doi.org/10.17182/hepdata.91374.v5/t27 (t27-t32)'
        # TChiWWoff_hino.validationTarball =  "TChi_nonDegenHinoLSP_BRN2toZN1100.tar.gz"
        TChiWWoff_hino.setSources(dataLabels = [ 'obsExclusion', 'expExclusion'],
                            dataFiles = [ "orig/HEPData-ins1767649-v5-Figure_14a_Observed.csv",
                                          "orig/HEPData-ins1767649-v5-Figure_14a_Expected.csv",
                                          ],
                            units = [ None ]*2,
                            dataFormats = [ 'csv' ]*2 )
        TChiWWoff_hino.addSource("efficiencyMap", (f"orig/interpolated_Figure28{datasets_ewk_incl[SR]['letter'][0]}.csv", f"orig/interpolated_Figure28{datasets_ewk_incl[SR]['letter'][1]}.csv"),
                          unit = "", dataFormat = "mcsv")
        TChiWWoff_hino.efficiencyMap._unit = "/242.0" # Acceptance unit is /10000 and BR(W*->e+ nu)*BR(W*->e- nu) + BR(W*->µ+ nu)*BR(W*->µ- nu)) = 0.01210 + 0.01210 =  0.02420


for SR in datasets_slep:
    #+++++++ dataset block ++++++++++++++
    dataset = DataSetInput(SR)
    if 'upperLimit' in datasets_slep[SR].keys():
        dataset.setInfo(dataType = 'efficiencyMap', dataId = SR,
                        observedN = datasets_slep[SR]['observedN'], expectedBG = datasets_slep[SR]['expectedBG'], bgError = datasets_slep[SR]['bgError'], upperLimit = datasets_slep[SR]['upperLimit'], expectedUpperLimit = datasets_slep[SR]['expectedUpperLimit'])
    else:
        dataset.setInfo(dataType = 'efficiencyMap', dataId = SR,
                        observedN = datasets_slep[SR]['observedN'], expectedBG = datasets_slep[SR]['expectedBG'], bgError = datasets_slep[SR]['bgError'])

    #+++++++ next mass plane block ++++++++++++++
    TSlepSlep                      = dataset.addTxName('TSlepSlep')
    TSlepSlep.checked              = 'no'
    if "2l_ee" in SR:
        TSlepSlep.constraint       = "[[['e+']],[['e-']]]"
    elif "2l_mm" in SR:
        TSlepSlep.constraint       = "[[['mu+']],[['mu-']]]"
    elif "CR" in SR:
        TSlepSlep.constraint       = "[[['e+']],[['e-']]] + [[['mu+']],[['mu-']]]"
    TSlepSlep.condition            = None
    TSlepSlep.massConstraint       = [['dm < 86.0'], ['dm < 76.0']]
    TSlepSlep.conditionDescription = None
    TSlepSlep.source               = "ATLAS"


    #+++++++ next mass plane block ++++++++++++++
    TSlepSlep                   = TSlepSlep.addMassPlane([[x, x-y], [x, x-y]])
    TSlepSlep.figure            = 'Extracted from the patchsets'
    TSlepSlep.figureUrl         = 'n/a'
    TSlepSlep.exclusionDataUrl  = 'https://doi.org/10.17182/hepdata.91374.v5/t9; https://doi.org/10.17182/hepdata.91374.v5/t10'
    TSlepSlep.dataUrl           = 'n/a'
    TSlepSlep.setSources(dataLabels = [ 'obsExclusion', 'expExclusion'],
                        dataFiles = [ "orig/HEPData-ins1767649-v5-Figure_16a_Observed.csv",
                                      "orig/HEPData-ins1767649-v5-Figure_16a_Expected.csv",
                                      ],
                        units = [ None ]*2,
                        dataFormats = [ 'csv' ]*2 )
    TSlepSlep.addSource("efficiencyMap", "orig/interpolated_" + datasets_slep[SR]['SRname'] + "_sleptons_degenerate_efficiency.csv",
                      unit = "", dataFormat = "csv")

databaseCreator.create()
