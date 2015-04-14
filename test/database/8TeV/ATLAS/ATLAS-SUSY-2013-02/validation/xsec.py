#!/usr/bin/python

from smodels.theory import crossSection
from smodels.tools.physicsUnits import pb, fb

## File = "TGQ_1435_100_1495_100.slha"
File = "TGQ_1372_100_1430_100.slha"

print "file",File

totalNLL=0*pb
totalLO=0*pb

List=crossSection.getXsecFromSLHAFile ( File )
Dict=List.getDictionary()
for key,value in Dict.items():
    if 1000021 in key and not key in [(1000021,1000021),(1000022,1000021),(1000021,1000022)]:
        print key,value
        if '8 TeV (NLL)' in value:
            totalNLL+=value['8 TeV (NLL)']
        elif '8 TeV (LO)' in value:
            totalNLL+=value['8 TeV (LO)']
        if '8 TeV (LO)' in value:
            totalLO+=value['8 TeV (LO)']

print totalNLL.asUnit(fb)
print totalLO.asUnit(fb)
