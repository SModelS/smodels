#!/usr/bin/python

"""
.. module:: AuxiliaryFunctions
    :synopsis: ...

.. moduleauthor:: someone <email@example.com>

"""



#Defines the auxiliar similar function
#returns the relative difference between any elements of the list normalized to [0,1]
def Csim(*els):
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


#Defines similar function when comparing two list of numbers.
#If any of the elements differ by < 10%, returns True
#FIXME should this go to ClusterTools?
def similar(els):
  for i in range(len(els)):
    for j in range(i+1,len(els)):
      if els[i] != els[j]:
        if 2.*abs(els[i]-els[j])/abs(els[i]+els[j]) > 0.1: return False
  return True

#Routine to evaluate the analyses conditions and constraints.
#Flexible version of eval to allow for additional operators,
#such as ~ (= similar)
def Ceval(instring,El):

  run = instring.replace(" ","")  #Remove blanks
  if "~" in run:
    simels = run.split("~")
    run = 'similar(' + str(simels) + ')'
    run = run.replace("'","")
  return eval(run)
