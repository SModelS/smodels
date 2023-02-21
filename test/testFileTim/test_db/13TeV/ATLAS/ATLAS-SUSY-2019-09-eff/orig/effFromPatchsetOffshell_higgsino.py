#!/usr/bin/env python3

import json
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.append('/home/pascal/SModelS/smodels-utils')
sys.path.append('/home/pascal/SModelS/smodels')
import pandas as pd
from smodels.tools.physicsUnits import fb, pb
from scipy.interpolate import griddata

with open('bkg_offshell.json', 'r') as f:
    bkg = json.load(f)
# Getting the paths of each SR from the background JSON
chPath = {}
iCh = 0
for ch in bkg['channels']:
    srName = ch['name']
    if srName[:2] == 'SR':
        chPath[srName] = '/channels/%d/samples/' % iCh
    iCh += 1
# Preparing the cross sections
lumi = 139/fb
n2c1m = pd.read_csv('/home/pascal/SModelS/smodels-utils/smodels_utils/morexsecs/tables/xsecN2C1mnondegenp13.txt', sep='\s+', comment='#', header=None)
n2c1m.columns = ['mN2', 'mC1', 'mN1', 'xsec', 'unc']
n2c1p = pd.read_csv('/home/pascal/SModelS/smodels-utils/smodels_utils/morexsecs/tables/xsecN2C1pnondegenp13.txt', sep='\s+', comment='#', header=None)
n2c1p.columns = ['mN2', 'mC1', 'mN1', 'xsec', 'unc']
plt.scatter(n2c1p['mN2'], n2c1p['mN2']-n2c1p['mN1'], color='C0', marker='+')
# Extracting the signals from the patchset and writing the exlusive eff. files
selectSR = 'SR' # allows to select a set of SRs based on a common substring
selectMass = [320., 300.]
selectEff = 0
with open('offshell_higgsino_patchset.json', 'r') as f:
    patchset = json.load(f)
for sr in chPath:
    name = "%s_higgsino_efficiency.csv" % sr
    # print(name)
    goodpath = chPath[sr]
    with open(name, 'w') as out:
        out.write('M(n2),M(n2)-M(n1),Efficiency\n')
        for pa in patchset['patches']:
            massv = pa['metadata']['values']
            minus = griddata(n2c1m[['mC1', 'mN1']], n2c1m['xsec'], massv)[0]
            if np.isnan(minus): minus = None
            plus = griddata(n2c1p[['mC1', 'mN1']], n2c1p['xsec'], massv)[0]
            if np.isnan(plus): plus = None
            xsec = None
            if plus and minus: # could be None if we are out of the interpolation hull
                xsec = plus + minus
            signal = None
            for op in pa['patch']:
                if goodpath in op['path']:
                    signal = op['value']['data'][0]
                    break
            if signal != None and xsec == None:
                plt.scatter(2*massv[0]-massv[1], 2*massv[0]-2*massv[1], color='tab:red', alpha=0.2)
            if signal != None and xsec != None:
                plt.scatter(2*massv[0]-massv[1], 2*massv[0]-2*massv[1], color='tab:orange', alpha=0.2)
                xsec *= pb
                eff = signal/xsec/lumi
                if eff > 1. :
                    print("ERROR. sr : ", sr," | signal : ", signal, " | xs : ", xsec, " | massv : ", massv)
                    continue
                if selectSR in sr and selectMass == massv:
                    print(sr)
                    selectEff += eff
                out.write("{},{},{}\n".format(2*massv[0]-massv[1], 2*(massv[0]-massv[1]), eff))
print('"""%s"""' % selectSR)
print("{},{},{}".format(2*selectMass[0]-selectMass[1], 2*(selectMass[0]-selectMass[1]), selectEff))
customlabels = \
[Line2D([0], [0], marker='+', color='tab:blue', label='tabulated xsecs'),
Line2D([0], [0], marker='o', color='tab:orange', label='signal and xsec'),
Line2D([0], [0], marker='o', color='tab:red', label='signal and no xsec')]
plt.legend(handles=customlabels)
plt.show()
