#!/usr/bin/python

from pyfeyn.user import *
import sys, random

processOptions()
#Pair production diagram
nvertices = 2
fd = FeynDiagram()
#R-odd vertices:
in_vtx = Circle(-2,0., radius=0.4).setFillStyle(HATCHED135)
vtxs = []
for i in range(nvertices/2):
  vtx = Vertex(-0.5+float(i)*1.5, 1., mark=CIRCLE)
  vtxs.append(vtx)
  vtx = Vertex(-0.5+float(i)*1.5, -1., mark=CIRCLE)
  vtxs.append(vtx)


#R-odd particles:
odds = []
odd1 = Fermion(in_vtx, vtxs[0])
odd2 = Fermion(in_vtx, vtxs[1])
odds.append(odd1)
odds.append(odd2)

for i in range(0,nvertices-2,2):
  odd = Fermion(vtxs[i],vtxs[i+2])
  odds.append(odd)
  odd = Fermion(vtxs[i+1],vtxs[i+3])
  odds.append(odd)

fd.draw("Dia0.pdf")

#General diagrams
random.seed(3.)
for ngraphs in range(5):
  fd = FeynDiagram()
  nvertices1 = int(random.uniform(2,4))
  nvertices2 = int(random.uniform(2,4))
#R-odd vertices:
  in_vtx = Circle(-2,0., radius=0.4).setFillStyle(HATCHED135)
  vtxs1 = []
  vtxs2 = []
  for i in range(nvertices1):
    vtx = Vertex(-0.5+float(i)*1.5, 1., mark=CIRCLE)
    vtxs1.append(vtx)
  for i in range(nvertices2):
    vtx = Vertex(-0.5+float(i)*1.5, -1., mark=CIRCLE)
    vtxs2.append(vtx)

#R-odd particles:
  odds = []
  odd1 = Fermion(in_vtx, vtxs[0])
  odd2 = Fermion(in_vtx, vtxs[1])
  odds.append(odd1)
  odds.append(odd2)

  for i in range(0,nvertices1-1):
    odd = Fermion(vtxs1[i],vtxs1[i+1])
    odds.append(odd)
  for i in range(0,nvertices2-1):
    odd = Fermion(vtxs2[i],vtxs2[i+1])
    odds.append(odd)

#R-even outer points
  outpts1 = []
  outpts2 = []
  for i in range(nvertices1):
    out = Point(-0.25 + float(i)*1.5,2)
    outpts1.append(out)
  for i in range(nvertices2):
    out = Point(-0.25 + float(i)*1.5,-2)
    outpts2.append(out)

#R-even particles
  evens = []
  for i in range(0,nvertices1): 
    disp = 0.2
    ran = random.uniform(0.,1.)
    if abs(ran) < 0.5:
      even = Fermion(vtxs1[i],Point(outpts1[i].getX()-0.25,outpts1[i].getY()))
      evens.append(even)
      even = Fermion(vtxs1[i],outpts1[i])
      evens.append(even)
    else:
      even = Fermion(vtxs1[i],outpts1[i])
      evens.append(even)

  for i in range(0,nvertices2):
    disp = 0.2
    ran = random.uniform(0.,1.)
    if abs(ran) < 0.5:  
      even = Fermion(vtxs2[i],Point(outpts2[i].getX()-0.25,outpts2[i].getY()))
      evens.append(even)
      even = Fermion(vtxs2[i],outpts2[i])
      evens.append(even)
    else:
      even = Fermion(vtxs2[i],outpts2[i])
      evens.append(even)


  fd.draw("Dia"+str(ngraphs+1)+".pdf")
















