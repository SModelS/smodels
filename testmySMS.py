#!/usr/bin/python

import mySMSDecomposition, os

topolist=os.listdir('regression')

topolist.remove('slha')

print "myTopo, Topo according to filename"

for topo in topolist:
	lhe=topo.split('_')[0]
    	print mySMSDecomposition.getSMS('%s' %lhe),' , ',lhe

