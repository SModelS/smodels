#!/usr/bin/env python3

import json
from re import search

SRlist = ['SRDF_0a_cuts','SRDF_0b_cuts','SRDF_0c_cuts','SRDF_0d_cuts','SRDF_0e_cuts','SRDF_0f_cuts','SRDF_0g_cuts','SRDF_0h_cuts','SRDF_0i_cuts','SRDF_1a_cuts','SRDF_1b_cuts','SRDF_1c_cuts','SRDF_1d_cuts','SRDF_1e_cuts','SRDF_1f_cuts','SRDF_1g_cuts','SRDF_1h_cuts','SRDF_1i_cuts','SRSF_0a_cuts','SRSF_0b_cuts','SRSF_0c_cuts','SRSF_0d_cuts','SRSF_0e_cuts','SRSF_0f_cuts','SRSF_0g_cuts','SRSF_0h_cuts','SRSF_0i_cuts','SRSF_1a_cuts','SRSF_1b_cuts','SRSF_1c_cuts','SRSF_1d_cuts','SRSF_1e_cuts','SRSF_1f_cuts','SRSF_1g_cuts','SRSF_1h_cuts','SRSF_1i_cuts']

with open('C1C1WW_patchset.json', 'r') as f:
    js = json.load(f)

dict={}
nbCR = 3
somme = 0
for pa in js['patches']:
    if pa["metadata"]["name"] == "C1C1_WW_500_100":
        for patch in pa["patch"]:
            srNumb = int(search('/channels/(.*)/samples/', patch["path"]).group(1)) - nbCR
            if srNumb >= 0:
                dict[SRlist[srNumb]] = patch["value"]["data"][0]
                somme += patch["value"]["data"][0]

f.close()

print(dict)
