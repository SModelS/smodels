#!/usr/bin/python

""" Closure: Tx -> SModelS description -> Tx """

import set_path
from Experiment import TxNames
from Theory import LHEReader, TopologyBuilder
from TestTools import ok

topolist = ['T1','T2','T1tttt', 'T2tt','T3W', 'T5WW', 'TChiWZ', 'T1bbbb', 'T2bb', 'T5WZ', 'T3Wb', 'T3Z', 'T5ZZ', 'T6bbZZ', 'TChiWW', 'TSlepSlep']

for topo in topolist:
  reader = LHEReader.LHEReader("../lhe/%s_1.lhe" % topo)
  Event = reader.next()
  SMSTop = TopologyBuilder.fromEvent(Event, {})
  res = TxNames.getTx(SMSTop[0].leadingElement())
  print "Checking %7s: %s" % ("["+res+"]",ok(res,topo))
