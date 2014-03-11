#!/usr/bin/python

import time,copy,sys
import element
import ParticleNames

ptcList = [[['nu'],['e-']],[['e+'],['ta-']]]
elementList = []
for ptc in ParticleNames.Reven.values():
  newList = copy.deepcopy(ptcList)
  newList[0][0] = [ptc]
  newElement = element.Element(str(newList))
  newElement.branches[0].masses = [300.+2.*len(elementList),200.,100.]
  newElement.branches[1].masses = [300.,200.,100.+1.*len(elementList)]
  elementList.append(newElement)


el = elementList[1]
elc = el.copy()
elc.branches = [elc.branches[1],elc.branches[0]]
elc.weight = elc.weight*2.
print el.getParticles(),el.getMasses(),el.weight
print elc.getParticles(),elc.getMasses(),elc.weight
sys.exit()

start = time.clock()
for el in elementList:
  for el2 in elementList:
    res =  el.isEqual(el2)
print 'time for isEqual=',time.clock()-start

start = time.clock()
for el in elementList:
  for el2 in elementList:
    res =  el.isEqual2(el2)
print 'time for isEqual2=',time.clock()-start
