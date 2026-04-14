#!/usr/bin/env python3

import json
import sys
sys.path.append('/Users/altakach/SModelS3/smodels-utils')
sys.path.append('/Users/altakach/SModelS3/smodels')
from smodels_utils.morexsecs.refxsecComputer import RefXSecComputer
from smodels_utils import SModelSUtils
from smodels.base.physicsUnits import fb, pb

lumi = 136/fb
xs = RefXSecComputer()

#++++++++load all the relevant cross-sections+++++++++++
ss = xs.getXSecsFrom(f'{SModelSUtils.installDirectory()}/smodels_utils/morexsecs/tables/xsecstop13.txt',columns={"mass":0, "xsec":1})

#list of all patchset json files to be read
Bkg_files = ['SRMET_patchset.json','SRMU_patchset.json']

#list of all file names to write the output efficiencies into
file_name = ['SRMET','SRMU']

#list of all the signal region names occuring in the relevant patchset files
#sr_name = ['SRC_mll', 'SRHigh_mll', 'SRLow_mll', 'SRMed_mll', 'SRZHigh_cuts','SRZLow_cuts','SRZMed_cuts']

#name of path in the patchset under 'patch'
goodpath = "/channels/0/samples/"

#name of the topologies in the patchset
topo_name = 'StopRHadron'

signal_events = []
#name of the topologies in the output file
#topos = ['T5ZZ','T6ZZ']
for j in range(len(Bkg_files)):
    signal_sr = 0
    #load the json file to be read
    file_json = open(f'orig/likelihoods/{Bkg_files[j]}', 'r')
    js = json.load(file_json)
    #no of values in 'data' under 'value' in 'patch' = no of bins , here only 1
    #output file of efficiencies for each bin
    name = 'orig/EffFromPatch_'+ file_name[j] + '.csv'
    out = open(name, 'w')
    out.write('# M_stop (GeV), Width (GeV), Efficiency\n')
    print("writing file")
    for pa in js['patches']:
        if topo_name in pa['metadata']['name']:
            #extract efficiencies
            massvector = pa['metadata']['values']
            mass = massvector[0]
            lifetime = massvector[1]*10**(-3)
            width = 6.582*10**(-16)/lifetime
            xsec = xs.interpolate(mass, ss)
            xsec *= pb
            for op in pa['patch']:
                if goodpath in op['path']:
                    signal = op['value']['data'][0]
                    signal_sr += signal
                    eff = signal/(xsec*lumi)
                    if eff > 1. :
                        print("ERROR. sr : ", name," | signal : ", signal, " | xs : ", xsec, " | massv : ", massvector)
                        continue
                    out.write("{} , {} , {}\n".format(mass, width, eff))
    signal_events.append(signal_sr)


print(signal_events)
