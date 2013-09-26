#!/usr/bin/env python

"""
.. module:: SMSAnalysisFactory
    :synopsis:  unit that creates a list of Analysis objects from a results database

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import SMSResults
from Tools.PhysicsUnits import rmvunit
import types

def getRealTopo ( Tx ):
  """ T3w025 -> T3w, etc """
  ret=Tx
  ret.replace("050","").replace("x1C180","").replace("025","")
  if ret.find("x")>-1: ret=ret[:ret.find("x")]
  return ret

def getArray ( constraint ):
  """ a helper method that is intended to be used to make it possible
      to extract the number of vertices, branches, and insertions from
      the constraint string. It maps e.g.
      2*([[['L'],['L']],[['L'],['nu']]] + [[['L'],['L']],[['nu'],['L']]])
      to
      [[['L'],['L']],[['L'],['nu']]] """
  # print "get array for",constraint
  c=constraint.replace(" ","")
  c=c.strip("1234567890.*%()")
  if c.find("+[")>5:
    c=c[:c.find("+[")]
  if c.find("-[")>5:
    c=c[:c.find("-[")]
  #print "constraint=",constraint,"c=",c
  ret=eval(c)
  return ret

def load( anas = None, topos=None, sqrts= [ 7, 8 ] ):
  """ This method creates the ana objects from the info given in the SMS
    results database.

    :param anas: if given as a list, then we create only objects for these analyses
       (the database naming convention is used).
    :param topos: if given as a list, then only these topos are considered.
    :param sqrts: array of center-of-mass energies of the analyses that are to be considered.
    :returns: list of Analyses

  """
  
  from analysis import Analysis, AnalysisPlot
  
  ## lets make sure the user can also supply a single topo/ana
  ## without having to code an array
  if type(topos)==types.StringType: topos = [ topos ]
  if type(anas)==types.StringType: anas = [ anas ]
  if type(sqrts)==types.IntType: sqrts = [ sqrts ]
  if type(sqrts)==types.FloatType: sqrts = [ int(sqrts) ]

  listOfAnalyses = []

  debug=False
  if debug:
    SMSResults.verbosity ( "info" )
  else:
    SMSResults.verbosity ( "error" )
  #anas=[ "alphaT", "RA48TeV", "RA2b8TeV", "Weakinos8TeV", "alphaT8TeV", "ATLAS_CONF_2013-024" ]
  if anas==None: anas=SMSResults.getAllResults().keys()
  for ana in anas:
    if debug:
      print
      print "Building ana",ana
    Ss=rmvunit(SMSResults.getSqrts( ana ),"TeV")
    if Ss==None:
      print "SS=",Ss,ana
      continue
    Ss=int(Ss)
    if not Ss in sqrts:
      continue
    for Tx in SMSResults.getTopologies(ana):
      if topos!=None and Tx not in topos: continue
      if debug: print Tx,
      # if Tx[:2]!="T2": continue
            
      newAnalysis = Analysis()
      newAnalysis.sqrts=SMSResults.getSqrts( ana )
      stopo=getRealTopo ( Tx )
      newAnalysis.label = ana
      newAnalysis.run = SMSResults.getRun ( ana ) ##  "2012"
      constraint=SMSResults.getConstraints ( ana, topo=stopo )
      if not constraint or constraint=="Not yet assigned":
        if debug:
          print "dont have a constraint for",ana,Tx,"(",stopo,")"
        continue     

      cond=SMSResults.getConditions ( ana, topo=stopo )
      if cond==None: cond=""
      andreconstraint=constraint.replace("'","").replace(" ","")
      andrecond=cond.replace("'","").replace(" ","")
      plot = AnalysisPlot()
      plot.result = andreconstraint
      plot.condition = andrecond
      plot.label = stopo
      newAnalysis.listOfPlots.append(plot)
      
      # ul=SMSResults.getUpperLimit ( ana, Tx, mx=400., my=100., interpolate=True )
#      print "[SMSAnalysisFactory] test for %s:%s ul(mx=400,my=100)=%s" % ( ana, Tx, ul )

  #Add ana to list of analyses:
      listOfAnalyses.append(newAnalysis)
      #print "done."

  return listOfAnalyses


if __name__ == "__main__":
  load()
  print "List of analyses/results: "
  listOfAnalyses=load()
  for (ct,ana) in enumerate(listOfAnalyses):
    print ct,ana.label,ana.Top.vertnumb,ana.Top.vertparts # .label,ana.sqrts
