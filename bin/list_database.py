#!/usr/bin/python

""" This is a simple tool that I (WW) need to work out with the CMS susy group
what digitized results are not yet published """

import set_path
from Experiment import SMSResults, SMSHelpers

def header ( File, latex=True ):
  File.write ( "\\begin{table}\n" )
  File.write ( "\\begin{tabular}{|c|c|c|}\n" )

def footer ( File, latex=True ):
  File.write ( "\\end{tabular}\n" )
  File.write ( "\\end{table}\n" )

def printResult ( File, analysis, topo, latex=True ):
  sqrts=SMSResults.getSqrts ( analysis )
  pas=SMSResults.getPAS ( analysis )
  constr=SMSResults.getConstraints ( analysis, topo )
  File.write ( "%s & %s & %s \\\\ \n" % ( analysis, topo, pas ) )

def run( latex=True ):
  File=open("database.tex","w")
  header ( File, latex )
  results=[]
  All=SMSResults.getAllResults()
  for analysis in All:
    for topo in All[analysis]:
      if not SMSResults.isPublic ( analysis ): continue
      constr=SMSResults.getConstraints ( analysis, topo )
      if not constr: continue
      run=SMSResults.getRun ( analysis )
      results.append ( [ analysis, topo ] )
  for result in results:
    printResult ( File, result[0], result[1], latex )
  footer ( File, latex )
  File.close()

run()
