#!/usr/bin/env python3

import json
from re import search

SRlist = ["SR_HM_Low_MCT", "SR_HM_Med_MCT", "SR_HM_High_MCT", "SR_MM_Low_MCT", "SR_MM_Med_MCT", "SR_MM_High_MCT", "SR_LM_Low_MCT", "SR_LM_Med_MCT", "SR_LM_High_MCT"]

with open('/home/pascal/SModelS/1Lbb-likelihoods/patchset.json', 'r') as f:
    js = json.load(f)

dict={}
nbCR = 5
somme = 0
for pa in js['patches']:
    if pa["metadata"]["name"] == "C1N2_Wh_hbb_350_100":
        for patch in pa["patch"]:
            srNumb = int(search('/channels/(.*)/samples/', patch["path"]).group(1)) - nbCR
            if srNumb >= 0:
                dict[SRlist[srNumb]] = patch["value"]["data"]
                somme += patch["value"]["data"][0]

f.close()

print(dict)
