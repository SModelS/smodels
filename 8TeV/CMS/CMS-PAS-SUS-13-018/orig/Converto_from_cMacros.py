import ROOT
import os,sys

files    = ['Exracted_T2bb_CMS-PAS-SUS-13-018.root']
data_names = ['combined_expLimit',
            'combined_expExclOneTimesProspino',
            'combined_expExclPlusOneSigmaProspino',
            'combined_expExclMinusOneSigmaProspino',
            'combined_obsExclOneTimesProspino',
            'combined_obsExclPlusSysErrProspino',
            'combined_obsExclMinusSysErrProspino']

output = 'Extracted_T2bb_CMS-PAS-SUS-13-018'



for File in files:
 for eff in data_names:

    print 'This is the eff', eff
    inFile = ROOT.TFile(File)
    canvas = inFile.Get('cCONT_')
    effi = canvas.GetPrimitive(eff)
    OUT = ROOT.TFile(eff+'.root','recreate')
    effi.Write()

