#!/usr/bin/env python

import os
from subprocess import call

files=os.listdir('./')
files=[f for f in files if 'T6ttZZ' in f]
print files

for f in files:
	call('cp ./%s ./%s' %(f,f.replace('T6ttZZ', 'T6ZZtt')), shell=True)
	print 'cp', f, '->', f.replace('T6ttZZ', 'T6ZZtt')

