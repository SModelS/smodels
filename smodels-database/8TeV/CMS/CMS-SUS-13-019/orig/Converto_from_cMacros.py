import ROOT
import os,sys

'''
files    = ['T1qqqq-SUS13019-final_XSEC.root','T2bb-SUS13019-final_XSEC.root','T2tt-SUS13019-final_XSEC.root',
            'T1bbbb-SUS13019-final_XSEC.root','T1tttt-SUS13019-final_XSEC.root','T2qq-SUS13019-final_XSEC.root',]

eff_name = ['XSec_limit_combined',   'gr_cTH_obs_combined',  'gr_mTH_obs_combined',
            'gr_pTH_obs_combined', 'gr_cTH_exp_combined', 'gr_cTH_psig_combined' , 'gr_cTH_msig_combined'           ]


for FILE in files:
    OUT = ROOT.TFile("Extracted_"+FILE,'NEW')
    OUT = ROOT.TFile("Extracted_"+FILE,'UPDATE')
    for DATA in eff_name:
       print 'This is the eff', DATA
       inFile = ROOT.TFile(FILE)
       canvas = inFile.Get('cXSEC_')
       effi = canvas.GetPrimitive(DATA)
       OUT = ROOT.TFile("Extracted_"+FILE,'UPDATE')
       effi.Write()
       OUT.Close()
'''


files    = ['T5WH-SUS13019-final_XSEC.root']
eff_name = ['h_Limit' , 'gr_Obs', 'gr_ObsPS', 'gr_ObsMS', 'gr_Exp', 'gr_ExpPS', 'gr_ExpMS' ]
for FILE in files:
    OUT = ROOT.TFile("Extracted_"+FILE,'NEW')
    OUT = ROOT.TFile("Extracted_"+FILE,'UPDATE')
    for DATA in eff_name:
        print 'This is the eff', DATA
        inFile = ROOT.TFile(FILE)
        canvas = inFile.Get('cXSEC_')
        effi = canvas.GetPrimitive(DATA)
        OUT = ROOT.TFile("Extracted_"+FILE,'UPDATE')
        effi.Write()
        OUT.Close()