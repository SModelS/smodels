#!/usr/bin/python

""" This is a simple tool that I (WW) need to work out with the CMS susy group
what digitized results are not yet published """

import set_path
from experiment import SMSResults, SMSHelpers
from tools.PhysicsUnits import rmvunit
  
Fields= [ "analysis", "sqrts", "lumi" ]
NiceName= { "nick name": "nick name", "analysis": "analysis", "sqrts":"$\sqrt{s}$", "topologies": "topologies", "constraints":"constraints", "lumi": "lumi" }


def isATLAS ( analysis ):
  if analysis.find("ATLAS")==0: return True
  return False

def begintable ( File ):
  File.write ( "\\begin{table}\n" )
  flds="{"
  for i in range(len(Fields)):
    flds+="l|"
  flds=flds[:-1]
  flds+="}"
  File.write ( "\\begin{center}\n" )
  File.write ( "\\begin{tabular}%s\n" % flds )
  for (ct,F) in enumerate(Fields):
    File.write ( "{\\bf %s} " % NiceName[F] )
    if ct<len(Fields)-1: File.write ( "&" )

  File.write ( "\\\\ \\hline \n" )


def header ( File ):
  #File.write ( "\\documentclass[8pt]{article}\n" )
  #File.write ( "\\usepackage{multirow}\n" )
  #File.write ( "\\begin{document}\n" )
  #File.write ( "\\thispagestyle{empty}\n" )
  begintable ( File )

def endtable ( File, caption, experiment ):
  File.write ( "\\end{tabular}\n" )
  File.write ( "\\caption{%s}\n" % caption )
  File.write ( "\\label{list_%s}\n" % experiment )
  File.write ( "\\end{center}\n" )
  File.write ( "\\end{table}\n" )

def footer ( File, caption, experiment ):
  endtable ( File, caption, experiment )
  #File.write ( "\\end{document}" )

def line ( File, caption, experiment ):
  # going from ATLAS to CMS
  endtable ( File, caption, experiment )
  #File.write ( "\\newpage\n \\thispagestyle{empty}\n" )
  begintable ( File )

def compileTex ( Filename="database.tex" ):
  import os
  os.system ( "pdflatex database.tex" )
  os.system ( "convert database.pdf database.png" )

def pprint ( constr ):
  constr=constr.replace("'","")
  if len(constr)>40:
    constr=constr[:37]+"..."
  return constr

def printFirstTopo ( File, analysis ):
  sqrts=int(rmvunit(SMSResults.getSqrts ( analysis ),"TeV"))
  lumi=SMSResults.getLumi ( analysis ).asNumber()
  pas=SMSResults.getPAS ( analysis )
  File.write ( "%s & %s & %s \\\\ \n" % ( pas, sqrts, lumi ) ) 

def printAnalysis ( File, analysis, topos ):
  for topo in topos:
    constr=SMSResults.getConstraints ( analysis, topo )
    if not constr or constr=="Not yet assigned": continue
    printFirstTopo ( File, analysis )
    break

def run( ):
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
  wasATLAS=True
  File=open("database.tex","w")
  header ( File )
  for analysis in keys:
    topos=results[analysis]
    if wasATLAS and not isATLAS ( analysis ):
      line( File, "List of ATLAS analyses used by SModelS", "ATLAS" )
      wasATLAS=False
    printAnalysis ( File, analysis, topos )
  footer ( File, "List of CMS analyses used by SModelS", "CMS" )
  File.close()
#  compileTex()

run()
