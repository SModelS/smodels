#!/usr/bin/env python3

import json
import sys
import math
sys.path.append('/Users/sahanan/smodels-utils')
sys.path.append('/Users/sahanan/smodels')
from smodels_utils.morexsecs.refxsecComputer import RefXSecComputer
from smodels_utils import SModelSUtils
from smodels.tools.physicsUnits import fb, pb



goodpath = "/channels/0/samples/"
sr_name = ["SRC_mll", "SRHigh_mll", "SRLow_mll", "SRMed_mll", "SRZHigh_cuts","SRZLow_cuts","SRZMed_cuts"]

#list of all patchset json files to be read
Bkg_files = ['STR-SRC_bkg.json','STR-SRHigh_bkg.json','STR-SRLow_bkg.json','STR-SRMed_bkg.json','STR-SRZHigh_bkg.json','STR-SRZLow_bkg.json','STR-SRZMed_bkg.json']

#list of all file names to write the output efficiencies into
file_name = ['SRC','SRHigh','SRLow','SRMed','SRZHigh','SRZLow','SRZMed']

for j in range(len(Bkg_files)):
    sum_expbg = 0.0
    sum_bgerr = 0.0
    sum_obs = 0
    #load json file
    with open(Bkg_files[j], 'r') as f:
        js = json.load(f)
        name = "orig/" + file_name[j] + '_Yield_Bin.csv'
        with open(name, 'w') as out:
            out.write('#Obs, ExpBg, Bgerr\n')
            #no of values in 'data' under 'value' in 'patch' = no of bins
            bin_size = len(js['channels'][0]['samples'][0]['data'])
            for i in range(bin_size):
                expbg = 0.0
                bgerr = 0.0
                for ch in js['channels']:
                    for pa in ch['samples']:
                        if sr_name[j] == ch['name']:
                            expbg += pa['data'][i]
                            sum_expbg += pa['data'][i]
                            for mo in pa['modifiers']:
                                if "staterror" in mo['name'] :
                                    bgerr += mo['data'][i]**2
                                    sum_bgerr += mo['data'][i]**2
                if sr_name[j] == js['observations'][0]['name'] :
                    obs = js['observations'][0]['data'][i]
                    sum_obs += obs
                out.write("{} , {} , {} \n".format(obs,expbg, math.sqrt(bgerr)))
            out.write("#{} , {} , {} \t Total Yield \n".format(sum_obs, sum_expbg, (math.sqrt(sum_bgerr))))


