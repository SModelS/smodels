#!/usr/bin/env python3

import json
import sys
sys.path.append('/home/pascal/SModelS/smodels-utils')
sys.path.append('/home/pascal/SModelS/smodels')
from smodels_utils.morexsecs.refxsecComputer import RefXSecComputer
from smodels.tools.physicsUnits import fb, pb

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
xs = RefXSecComputer()
n2c1m = xs.getXSecsFrom('/home/pascal/SModelS/smodels-utils/smodels_utils/morexsecs/tables/xsecN2C1m13.txt')
n2c1p = xs.getXSecsFrom('/home/pascal/SModelS/smodels-utils/smodels_utils/morexsecs/tables/xsecN2C1p13.txt')
# Extracting the signals from the patchset and writing the exlusive eff. files
# selectSR = 'high_0J'
# selectMass = [150., 60.]
# selectEff = 0
with open('offshell_winobino_plus_patchset.json', 'r') as f:
    patchset = json.load(f)
for sr in chPath:
    name = "%s_winobino(+)_efficiency.csv" % sr
    goodpath = chPath[sr]
    with open(name, 'w') as out:
        out.write('M(c1,n2),M(c1,n2)-M(n1),Efficiency\n')
        for pa in patchset['patches']:
            massvector = pa['metadata']['values']
            xsec = None
            plus = xs.interpolate(massvector[0], n2c1m)
            minus = xs.interpolate(massvector[0], n2c1p)
            if plus and minus: # could be None if we are out of the interpolation hull
                xsec = plus + minus
            signal = None
            for op in pa['patch']:
                if goodpath in op['path']:
                    signal = op['value']['data'][0]
            if signal != None and xsec != None:
                xsec *= fb
                eff = signal/xsec/lumi
                if eff > 1. :
                    print("ERROR. sr : ", sr," | signal : ", signal, " | xs : ", xsec, " | massv : ", massvector)
                    continue
                # if selectSR in sr and selectMass == massvector:
                #     print(sr)
                #     selectEff += eff
                # print("{},{},{}".format(massvector[0], massvector[1], eff))
                out.write("{},{},{}\n".format(massvector[0], massvector[0]-massvector[1], eff))
# print('"""%s"""' % selectSR)
# print("{},{},{}".format(selectMass[0], selectMass[0]-selectMass[1], selectEff))
