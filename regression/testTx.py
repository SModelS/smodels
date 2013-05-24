#!/usr/bin/python

""" starting from a Tx.lhe file, we SmodelS-decompose the event, and 
obtain its Tx name from the SModelS description of the event: a closure test. """

import set_path
from Experiment import TxNames
from Theory import SMSmethods, LHEReader, TopologyBuilder
from TestTools import ok

topolist = ['T1','T2','T1tttt', 'T2tt','T3W', 'T5WW', 'TChiWZ', 'T1bbbb', 'T2bb', 'T5WZ', 'T3Wb', 'T3Z', 'T5ZZ', 'T6bbZZ', 'TChiWW', 'TSlepSlep']

for topo in topolist:
        File=open("%s_1.lhe" % topo)
        reader = LHEReader.LHEReader("%s_1.lhe" % topo)
        Event = reader.next()
        SMSTop = TopologyBuilder.fromEvent(Event, {})
	res = TxNames.getTx(SMSTop[0].leadingElement())
	print "Checking %7s: %s" % ("["+res+"]",ok(res,topo))
