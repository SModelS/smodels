import ROOT

'''
    This script extracts the efficiency maps for the CMS-SUS-13-013 analysis from the ModelA1 root file
    Eff are in percetnage, so it divides correctly the numbers.
    Federico A.
'''

FILE = ROOT.TFile('ModelA1.root')

acc = ['SR21_acceptance','SR22_acceptance','SR23_acceptance','SR24_acceptance','SR25_acceptance','SR26_acceptance','SR27_acceptance','SR28_acceptance']
lab = ['SR21_HighPt.txt','SR22_HighPt.txt','SR23_HighPt.txt','SR24_HighPt.txt','SR25_HighPt.txt','SR26_HighPt.txt','SR27_HighPt.txt','SR28_HighPt.txt']

Glu = [ x*25 for x in range (0,50)]
Neu = [ x*25 for x in range (0,35)]


for ACC,LAB in zip(acc,lab):
    histo = FILE.Get(ACC)
    OUT = open(LAB,'w')
    OUT.write('# SR '+ LAB[:-4] + ' Gluino , Neutralino , Acc \n')
    for glu in Glu:
        for neu in Neu:
         if(glu > neu):
            eff = float(histo.GetBinContent(histo.FindBin(glu,neu)) )/100
            OUT.write(str(glu) + ' ' + str(neu) + ' ' + str(eff) +'\n')


