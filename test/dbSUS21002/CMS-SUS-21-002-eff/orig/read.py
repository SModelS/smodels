#!/usr/bin/env python3

import uproot, IPython

f=uproot.open("AllBinAccXEff_bTag_TChiWZ.root")
histo = f.get("AccXEff_WHSR")
for xi,x in enumerate ( histo.axis(0).centers() ):
    for yi,y in enumerate ( histo.axis(1).centers() ):
        for zi,z in enumerate ( histo.axis(2).centers() ):
            v = histo.values()[xi][yi][zi]
            if y < 40:
            # if abs(v)<1e-10:
                print ( f"x,y,z={x},{y},{z} value is: {v}" )

IPython.embed(colors="neutral")
