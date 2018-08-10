import ROOT
import os,sys

#files    = ["T1bbbb_exclusion_corrected.root"]
#data_names = ["graph_smoothed_Obs_T1bbbb"]

files    = ["T1tttt_exclusion_corrected.root"]
data_names = ["graph_smoothed_Obs_T1tttt"]


#First run the cMacro file to root (root[0] .x file.C) and save it to ROOT format
# (file.root). Then add file.root to the files list above and add the object names
#to data_names


for File in files:
 for eff in data_names:

    print 'This is the eff', eff
    inFile = ROOT.TFile(File)
    canvas = inFile.Get('cCONT_')
    effi = canvas.GetPrimitive(eff)
    OUT = ROOT.TFile(eff+'.root','recreate')
    effi.Write()

