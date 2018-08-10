#!/usr/bin/env python

import os
from subprocess import call

files=os.listdir('./')
files=[f for f in files if 'TChiChipmHW' in f]
print files

for f in files:
	call('cp ./%s ./%s' %(f,f.replace('TChiChipmHW', 'TChiWH')), shell=True)
	print 'cp', f, '->', f.replace('TChiChipmHW', 'TChiWH')

