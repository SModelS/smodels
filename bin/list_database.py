#!/usr/bin/python

""" This is a simple tool that I (WW) need to work out with the CMS susy group
what digitized results are not yet published """

import set_path
from Experiment import SMSResults, SMSHelpers
from Tools.PhysicsUnits import rmvunit
  
fields=3

def header ( File, latex=True ):
  File.write ( "\\documentclass{article}\n" )
  File.write ( "\\begin{document}\n" )
  File.write ( "\\thispagestyle{empty}\n" )
  File.write ( "\\begin{table}\n" )
  flds="{"
  for i in range(fields):
    flds+="c|"
  flds=flds[:-1]
  flds+="}"
  File.write ( "\\begin{tabular}%s\n" % flds )
  if fields==2:
    File.write ( "{\\bf analysis} & {\\bf topologies} \\\\ \\hline \n" )
  else:
    File.write ( "{\\bf analysis} & $\\sqrt{s}$ & {\\bf topologies} \\\\ \\hline \n" )

def footer ( File, latex=True ):
  File.write ( "\\end{tabular}\n" )
  File.write ( "\\end{table}\n" )
  File.write ( "\\end{document}" )

def compileTex ( Filename="database.tex" ):
  import os
  os.system ( "pdflatex database.tex" )
  os.system ( "convert database.pdf database.png" )


def printResult ( File, analysis, topos, latex=True ):
  sqrts=int(rmvunit(SMSResults.getSqrts ( analysis ),"TeV"))
  pas=SMSResults.getPAS ( analysis )
  Topos=""

  ctr=0
  for topo in topos:
    constr=SMSResults.getConstraints ( analysis, topo )
    sqrts=SMSResults.getSqrts ( analysis, topo )
    if ctr>2:
      Topos+="\\\\ & "
      ctr=0
    if constr: 
      Topos+=topo+", "
      ctr+=1
    
  if len(Topos)>2: Topos=Topos[:-2]
  if fields==2:
    File.write ( "%s & %s \\\\ \n" % ( pas, Topos ) )
  else:
    File.write ( "%s & %s & %s \\\\ \n" % ( pas, sqrts, Topos ) )
  # File.write ( "%s & %s & %s \\\\ \n" % ( analysis, topo, pas ) )

def run( latex=True ):
  File=open("database.tex","w")
  header ( File, latex )
  results={}
  All=SMSResults.getAllResults()
  for analysis in All:
    for topo in All[analysis]:
      if not SMSResults.isPublic ( analysis ): continue
      constr=SMSResults.getConstraints ( analysis, topo )
      if not constr: continue
      sqrts=SMSResults.getSqrts ( analysis )
      run=SMSResults.getRun ( analysis )
      if not results.has_key ( analysis ): results[analysis]=[]
      results[analysis].append ( topo )

  keys=results.keys()
  keys.sort()
  print keys
  for analysis in keys:
    topos=results[analysis]
    printResult ( File, analysis, topos, latex )
  footer ( File, latex )
  File.close()
  compileTex()

run()
