#!/usr/bin/python

import set_path, argparse, types

argparser = argparse.ArgumentParser(description='simple tool that is meant to draw lessagraphs, as a pdf feynman plot') 
argparser.add_argument ( '-T', nargs='?', help='Tx name, will look up lhe file in ../regression/Tx_1.lhe. Will be overriden by the "--lhe" argument', type=types.StringType, default='T1' )
argparser.add_argument ( '-l', '--lhe', nargs='?', help='lhe file name, supplied directly. Takes precedence over "-T" argument.', type=types.StringType, default='' )
argparser.add_argument ( '-o', '--output', nargs='?', help='output file, can be pdf or eps', type=types.StringType, default='out.pdf' )
args=argparser.parse_args()

from Theory import LHEReader, SMSmethods, TopologyBuilder
from Tools import SMSFeynmanGraphs

filename="../lhe/%s_1.lhe" % args.T
if args.lhe!="": filename=args.lhe

reader = LHEReader.LHEReader( filename )
Event = reader.next()
SMSTop = TopologyBuilder.fromEvent(Event, {} )
SMSFeynmanGraphs.draw ( SMSTop[0].leadingElement(), args.output )
