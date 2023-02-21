#!/usr/bin/env python3

import json
import sys
sys.path.append('/Users/sahanan/smodels-utils')
sys.path.append('/Users/sahanan/smodels')
from smodels_utils.morexsecs.refxsecComputer import RefXSecComputer
from smodels_utils import SModelSUtils
from smodels.tools.physicsUnits import fb, pb

#++++++load the json file to be read++++++++
with open('llh_jsons/SRC_patchset.json', 'r') as f:
    js = json.load(f)

lumi = 139/fb
xs = RefXSecComputer()

#++++++++load all the relevant cross-sections+++++++++++
#gg = xs.getXSecsFrom(f'{SModelSUtils.installDirectory()}/smodels_utils/morexsecs/tables/xsecgluino13.txt',columns={"mass":0, "xsec":2})
#ss = xs.getXSecsFrom(f'{SModelSUtils.installDirectory()}/smodels_utils/morexsecs/tables/xsecsquark13.txt',columns={"mass":0, "xsec":2})
c1n2 = xs.getXSecsFrom(f'{SModelSUtils.installDirectory()}/smodels_utils/morexsecs/tables/xsecN2C1pm13.txt',columns={"mass":0, "xsec":1})
cross_sections = [c1n2]

#list of all patchset json files to be read
Bkg_files = ['ewk_signal_patchset.json']

#list of all file names to write the output efficiencies into
file_name_1 = ["SRHigh4_cuts", "SRHigh8_1_cuts", "SRHigh8_2_cuts", "SRllbb_cuts", "SRInt_1_cuts", "SRInt_2_cuts", "SRLow_2_cuts"]
file_name_2 = ["SRHigh16_1_cuts", "SRHigh16_2_cuts", "SRLow_1_cuts", "SRLow2_cuts", "SROffShell_1_cuts", "SROffShell_2_cuts"]
#list of all the signal region names occuring in the relevant patchset files
#sr_name = ['SRC_mll', 'SRHigh_mll', 'SRLow_mll', 'SRMed_mll', 'SRZHigh_cuts','SRZLow_cuts','SRZMed_cuts']

#name of path in the patchset under 'patch'
goodpath_1 = ["/channels/4/samples/7/data/0", "/channels/5/samples/7/data/0", "/channels/6/samples/7/data/0",  "/channels/9/samples/7/data/0", "/channels/10/samples/7/data/0", "/channels/11/samples/7/data/0",  "/channels/13/samples/7/data/0" ]
goodpath_2 = ["/channels/7/samples/8", "/channels/8/samples/8", "/channels/12/samples/8", "/channels/14/samples/8", "/channels/15/samples/8", "/channels/16/samples/8"]
#name of the topologies in the patchset
topo_name = ['C1N2']

#name of the topologies in the output file
topos = ['TChiWZ']


for j in range(len(file_name_1)):
    with open(f'llh_jsons/{Bkg_files[0]}', 'r') as file:
        js = json.load(file)
        name = 'orig/EffFromPatch_'+ file_name_1[j] + '_' + topos[0] + '.csv'
        with open(name, 'w') as out:
            out.write('# M(C1), M(LSP), Efficiency\n')
            for pa in js['patches']:
                if topo_name[0] in pa['metadata']['values'][0]:
                    #extract efficiencies
                    massvector = pa['metadata']['values']
                    xsec = xs.interpolate(massvector[1], cross_sections[0])
                    xsec *= fb
                    for op in pa['patch']:
                        if goodpath_1[j] in op['path']:
                            signal = op['value']
                            eff = signal/xsec/lumi
                            if eff > 1. :
                                print("ERROR. sr : ", name," | signal : ", signal, " | xs : ", xsec, " | massv : ", massvector)
                                continue
                            out.write("{} , {} , {}\n".format(massvector[1], massvector[2], eff))


for k in range(len(topos)):
    for j in range(len(file_name_2)):
        #load the json file to be read
        with open(f'llh_jsons/{Bkg_files[0]}', 'r') as file:
            js = json.load(file)
            #no of values in 'data' under 'value' in 'patch' = no of bins
            bin_size = len(js['patches'][0]['patch'][0]['value']['data'])
            for i in range(bin_size):
                #output file of efficiencies for each bin
                name = 'orig/EffFromPatch_'+ file_name_2[j] + '_' + topos[k]+'_Bin_'+ str(i+1) + '.csv'
                with open(name, 'w') as out:
                    out.write('# M(C1), M(LSP), Efficiency\n')
                    for pa in js['patches']:
                        if topo_name[k] in pa['metadata']['values'][0]:
                            #extract efficiencies
                            massvector = pa['metadata']['values']
                            xsec = xs.interpolate(massvector[1], cross_sections[k])
                            xsec *= fb
                            for op in pa['patch']:
                                if goodpath_2[j] in op['path'] and "add" in op['op']:
                                    signal = op['value']['data'][i]
                                    eff = signal/xsec/lumi
                                    if eff > 1. :
                                        print("ERROR. sr : ", name," | signal : ", signal, " | xs : ", xsec, " | massv : ", massvector)
                                        continue
                                    out.write("{} , {} , {}\n".format(massvector[1], massvector[2], eff))

