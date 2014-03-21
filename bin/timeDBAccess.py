#!/usr/bin/python

import set_path, argparse, types
from experiment import smsResults, smsHelpers
from tools import rcFile
from timeit import timeit

argparser = argparse.ArgumentParser(description='simple tool that times the database access') 
argparser.add_argument ( '-d', '--database', nargs='?', help='database directory', type=types.StringType, default='' )
argparser.add_argument ( '-n', '--nqueries', nargs='?', help='number of queries', type=types.IntType, default=10 )
args=argparser.parse_args()

if args.database!="":
  smsHelpers.Base=args.database
n=args.nqueries

results={ "n": 0 }

def query():
  res=smsResults.getAllResults()

  for  (ana,topos) in res.items():
    for topo in topos:
      constraint=smsResults.getConstraints ( ana, topo )
      condition=smsResults.getConditions ( ana, topo )
      results["n"]+=1
      # print ana,topo,constraint,condition

T=timeit(query, number=n )
print "n=%d base=%s" % ( n, smsHelpers.Base )
print "%d results in %f seconds" % ( results["n"],T )
