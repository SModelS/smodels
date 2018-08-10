import ROOT
import os,sys

#files    = ['T1bbbb_exclusion_corrected.root']
#data_names = ["hXsec_exp_corr","graph_smoothed_Obs_T1bbbb"]
files    = ["T1tttt_exclusion_corrected.root"]
data_names = ["hXsec_obs_final","graph_smoothed_Obs_T1tttt"]



for File in files:
 for eff in data_names:

    print 'This is the eff', eff
    inFile = ROOT.TFile(File)
    canvas = inFile.Get('cCONT_')
    effi = canvas.GetPrimitive(eff)
    OUT = ROOT.TFile(eff+'.root','recreate')
    effi.Write()

