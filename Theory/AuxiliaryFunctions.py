#!/usr/bin/env python

"""
.. module:: AuxiliaryFunctions
    :synopsis: A collection of functions used in the conditions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""



def Csim(*els):
  """ Defines the auxiliar similar function
    returns the relative difference between any elements 
    of the list normalized to [0,1]. FIXME relative difference of what? """
  res = 0.
  for i in range(len(els)-1):
    a = els[i]
    for j in range(i,len(els)):
      b = els[j]
      if a == b: continue
      res = max(res,abs(a-b)/abs(a+b))
  return res

#Defines the auxiliary greater function
#returns a number between 0 and 1 depending on how much it is violated (0 = A > B, 1 = A << B)
def Cgtr(a,b):
  from Tools.PhysicsUnits import addunit
  if type(a) == type(addunit(1.,'GeV')) and a.asNumber() + b.asNumber()  == 0.: return 'N/A'
  if type(a) != type(addunit(1.,'GeV')) and a + b == 0.: return 'N/A'
  res = (abs(a-b) - (a-b))/(2.*(a+b))
  return res


def similar(els):
  """ Defines similar function when comparing two list of numbers.
      If any of the elements differ by < 10%, returns True
      FIXME should this go to ClusterTools? """
  for i in range(len(els)):
    for j in range(i+1,len(els)):
      if els[i] != els[j]:
        if 2.*abs(els[i]-els[j])/abs(els[i]+els[j]) > 0.1: return False
  return True

def Ceval(instring,nEl):
  """ Routine to evaluate the analyses conditions and constraints.
      Flexible version of eval to allow for additional operators,
      such as ~ (= similar) """
   
  run = instring.replace(" ","")  #Remove blanks
  if "~" in run:
    simels = run.split("~")
    run = 'similar(' + str(simels) + ')'
    run = run.replace("'","")
  if not run or run=="None": return None
  return eval(run)

def getelements(instring):
  """ Parses instring (or a list of strings) and return a list of elements (in string format) appearing
  in instring"""
  import copy, SMSDataObjects
  from ParticleNames import Reven, PtcDic
  
  if type(instring) == type('st'):
    outstr = copy.deepcopy(instring)
  elif type(instring) == type([]):
    outstr = ""
    for st in instring:
      if type(st) != type('st'):
        print "getlements: Input must be a string or a list of strings"
        return False    
      outstr += st
  
  elements = []
  

  outstr = outstr.replace(" ","")
  while "[[[" in outstr:  #String has element
      st = outstr[outstr.find("[[["):outstr.find("]]]")+3] #Get duplet
      element = SMSDataObjects.EElement ( st )
      ptclist = element.allParticles()
  #Syntax checks:
      for ib in range(2):
        for ptcL in ptclist[ib]:
          for ptc in ptcL:
            if not ptc in Reven.values() and not PtcDic.has_key(ptc):
              print "EvalRes: Unknown particle ",ptc
              return False
      outstr = outstr.replace(st,"")  # delete element
      elements.append(st)   #Store elements
      
  return elements    

def eltonum(instring,dic):
  """Replace all elements in string by their respective numerical values given in the dictionary dic"""
  
  outstr = instring
  for key in dic.keys():
    outstr = outstr.replace(key,str(dic[key]))
  
  
  if len(getelements(outstr)) > 0:
    print "eltonum: Missing dictionary entries for ",str(getelements(outstr))
    return False
  
  return outstr
    
