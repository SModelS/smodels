#!/usr/bin/env python3

import uproot

f=uproot.open("CMS-SUS-20-001_Figure_012-a.root")
t=f["theory_xsec"]
x,y = t.values()
f.close()
f=open("myC1C113.txt","wt" )
f.write ( "# C1C1 production xsecs @ 13 TeV\n " )
f.write ( "# extracted from CMS-SUS-20-001_Figure_012-a.root\n" )
f.write ( "# mC1 [GeV]  xsec [fb]\n" )
for xi,yi in zip ( x, y ):
    f.write ( f"{xi} {1000.*yi:.2f}\n" )
f.close()
