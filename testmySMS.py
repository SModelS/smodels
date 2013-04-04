#!/usr/bin/python

import mySMSDecomposition, os

green='\033[0;32m'
red='\033[0;31m'
reset='\033[;0m'

#topolist=os.listdir('regression')

#topolist.remove('slha')

topolist = ['T1','T2','T1tttt', 'T2tt', 'T5ww', 'TChiwz']

#print "myTopo, Topo according to filename"

def ok ( A,B ):
  if A==B: return "%sok.%s" % ( green, reset )
  return "%sfailed. [%s]%s" % ( red, B, reset )

for topo in topolist:
#	lhe=topo.split('_')[0]
#    	print mySMSDecomposition.getSMS('%s' %topo),' , ',topo
	res = mySMSDecomposition.getSMS('%s' %topo)
	print "Checking %7s: %s" % ("["+res+"]",ok(res,topo))
