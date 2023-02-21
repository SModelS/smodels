#!/usr/bin/env python3

import json
import sys
sys.path.append('/home/pascal/SModelS/smodels-utils')
sys.path.append('/home/pascal/SModelS/smodels')
from smodels_utils.morexsecs.refxsecComputer import RefXSecComputer
from smodels.tools.physicsUnits import fb, pb

with open('onshell_winobino_plus_patchset.json', 'r') as f:
    js = json.load(f)

lumi = 139/fb
xs = RefXSecComputer()
n2c1m = xs.getXSecsFrom('/home/pascal/SModelS/smodels-utils/smodels_utils/morexsecs/tables/xsecN2C1m13.txt')
n2c1p = xs.getXSecsFrom('/home/pascal/SModelS/smodels-utils/smodels_utils/morexsecs/tables/xsecN2C1p13.txt')

for iSR in range(1, 21):
    name = "SRWZ_%d_efficiency.csv" % iSR
    goodpath = "/channels/%d/samples/" % (iSR+2)
    with open(name, 'w') as out:
        out.write('M(c1,n2),M(n1),Efficiency\n')
        for pa in js['patches']:
            massvector = pa['metadata']['values'] #[C1(=N2),N1]
            xsec = xs.interpolate(massvector[0], n2c1m) + xs.interpolate(massvector[0], n2c1p)
            xsec *= fb
            for op in pa['patch']:
                if goodpath in op['path']:
                    signal = op['value']['data'][0]
            eff = signal/xsec/lumi
            if eff > 1. :
                print("ERROR. sr : ", name," | signal : ", signal, " | xs : ", xsec, " | massv : ", massvector)
                continue
            out.write("{},{},{}\n".format(massvector[0], massvector[1], eff))
