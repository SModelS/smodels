#!/usr/bin/env python3

import json
import sys
sys.path.append('/home/pascal/SModelS/smodels-utils')
sys.path.append('/home/pascal/SModelS/smodels')
from smodels_utils.morexsecs.refxsecComputer import RefXSecComputer
from smodels.tools.physicsUnits import fb, pb

with open('C1C1WW_patchset.json', 'r') as f:
    js = json.load(f)

lumi = 139/fb
xs = RefXSecComputer()
c1c1 = xs.getXSecsFrom('/home/pascal/SModelS/smodels-utils/smodels_utils/morexsecs/tables/xsecC1C113.txt')


srb = ['SR-DF-0J','SR-DF-1J','SR-SF-0J','SR-SF-1J']
mt2 = ['100,105','105,110','110,120','120,140','140,160','160,180','180,220','220,260','260,inf']

iSR = 2
for i in srb:
    for j in mt2:
        name = 'EffFromPatch_' + i + '_mt2=[' + j + ').txt'
        iSR += 1
        goodpath = "/channels/%d/samples/" % iSR
        with open(name, 'w') as out:
            out.write('# M(c1), M(n1), Efficiency\n')
            for pa in js['patches']:
                massvector = pa['metadata']['values']
                xsec = xs.interpolate(massvector[0], c1c1)
                xsec *= fb
                for op in pa['patch']:
                    if goodpath in op['path']:
                        signal = op['value']['data'][0]
                eff = signal/xsec/lumi
                if eff > 1. :
                    print("ERROR. sr : ", name," | signal : ", signal, " | xs : ", xsec, " | massv : ", massvector)
                    continue
                out.write("{} {} {}\n".format(massvector[0], massvector[1], eff))
