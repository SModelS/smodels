#!/usr/bin/python

from pyfeyn.user import *
import sys, random

processOptions()
fd = FeynDiagram()

#R-odd vertices:
in_vtx = Circle(-2,0., radius=0.4).setFillStyle(HATCHED135)
vtx1 = Vertex(-0.5, 1., mark=CIRCLE)
vtx2 = Vertex(0.5, -1., mark=CIRCLE)
vtx3 = Vertex(1., 1., mark=CIRCLE)
vtx4 = Vertex(2.5, -1., mark=CIRCLE)
vtx5 = Vertex(2.5, 1., mark=CIRCLE)
#R-odd particles:
odd1 = Fermion(in_vtx, vtx1).addLabel("$\\tilde{\chi}_1^+$", pos=0.5, displace=-0.2)
odd2 = Fermion(in_vtx, vtx2).addLabel("$\\tilde{\chi}_2^0$", pos=0.5, displace=0.2)
odd3 = Sfermion(vtx1, vtx3).addLabel("$\\tilde{\\nu}$", pos=0.5, displace=-0.2)
odd4 = Fermion(vtx2, vtx4).addLabel("$\\tilde{\chi}_1^0$", pos=0.5, displace=0.2)
odd5 = Fermion(vtx3, vtx5).addLabel("$\\tilde{\chi}_1^0$", pos=0.5, displace=-0.2)

#R-even outer points
out1 = Point(-0.25,  2)
out2 = Point(0.75, -2)
out3 = Point(1.25,  2)
out4 = Point(0.25, -2)
#R-even particles
even1 = Fermion(vtx1, out1).addLabel("$l^+$", pos=0.7, displace=-0.2)
even2 = Fermion(vtx2, out2).addLabel("$l+$", pos=0.7, displace=-0.2)
even3 = Fermion(vtx3, out3).addLabel("$\\nu$", pos=0.7, displace=-0.2)
even4 = Fermion(vtx2, out4).addLabel("$l^-$", pos=0.7, displace=0.2)

fd.draw("C1N2.pdf")


fd = FeynDiagram()
#R-odd vertices:
in_vtx = Circle(-2,0., radius=0.4).setFillStyle(HATCHED135)
vtx1 = Vertex(-0.5, 1., mark=CIRCLE)
vtx2 = Vertex(0.5, -1., mark=CIRCLE)
vtx3 = Vertex(1., 1., mark=CIRCLE)
vtx4 = Vertex(2.5, -1., mark=CIRCLE)
vtx5 = Vertex(2.5, 1., mark=CIRCLE)
#R-odd particles:
odd1 = Sfermion(in_vtx, vtx1).addLabel("$\\tilde{l}^+$", pos=0.5, displace=-0.2)
odd2 = Sfermion(in_vtx, vtx2).addLabel("$\\tilde{\\nu}_2$", pos=0.5, displace=0.2)
odd3 = Fermion(vtx1, vtx3).addLabel("$\\tilde{\chi}_2^0$", pos=0.5, displace=-0.2)
odd4 = Sfermion(vtx2, vtx4).addLabel("$\\tilde{\\nu}_1$", pos=0.5, displace=0.2)
odd5 = Sfermion(vtx3, vtx5).addLabel("$\\tilde{\\nu}_1$", pos=0.5, displace=-0.2)

#R-even outer points
out1 = Point(-0.25,  2)
out2 = Point(0.75, -2)
out3 = Point(1.25,  2)
out4 = Point(0.25, -2)
#R-even particles
even1 = Fermion(vtx1, out1).addLabel("$l^+$", pos=0.7, displace=-0.2)
even2 = Fermion(vtx2, out2).addLabel("$l+$", pos=0.7, displace=-0.2)
even3 = Fermion(vtx3, out3).addLabel("$\\nu$", pos=0.7, displace=-0.2)
even4 = Fermion(vtx2, out4).addLabel("$l^-$", pos=0.7, displace=0.2)

fd.draw("SLSN.pdf")



fd = FeynDiagram()
#R-odd vertices:
in_vtx = Circle(-2,0., radius=0.4).setFillStyle(HATCHED135)
vtx1 = Vertex(-0.5, 1., mark=CIRCLE)
vtx2 = Vertex(0.5, -1., mark=CIRCLE)
vtx3 = Vertex(1., 1., mark=CIRCLE)
vtx4 = Vertex(2.5, -1., mark=CIRCLE)
vtx5 = Vertex(2.5, 1., mark=CIRCLE)
#R-odd particles:
odd1 = Fermion(in_vtx, vtx1)
odd2 = Fermion(in_vtx, vtx2)
odd3 = Fermion(vtx1, vtx3)
odd4 = Fermion(vtx2, vtx4)
odd5 = Fermion(vtx3, vtx5)

#R-even outer points
out1 = Point(-0.25,  2)
out2 = Point(0.75, -2)
out3 = Point(1.25,  2)
out4 = Point(0.25, -2)
#R-even particles
even1 = Fermion(vtx1, out1).addLabel("$l^+$", pos=0.7, displace=-0.2)
even2 = Fermion(vtx2, out2).addLabel("$l+$", pos=0.7, displace=-0.2)
even3 = Fermion(vtx3, out3).addLabel("$\\nu$", pos=0.7, displace=-0.2)
even4 = Fermion(vtx2, out4).addLabel("$l^-$", pos=0.7, displace=0.2)

fd.draw("C1N2naked.pdf")


fd = FeynDiagram()

#l1 = Label("Topology Label: [ [[$l^+$],[$\\nu$]] , [[$l^+$,$l^-$]] ]", x=0, y=2.5)
#R-odd vertices:
in_vtx = Circle(-2,0., radius=0.4).setFillStyle(HATCHED135)
vtx1 = Vertex(-0.5, 1., mark=CIRCLE).addLabel("[$l^+$]", displace=-0.4, angle=90)
vtx2 = Vertex(0.5, -1., mark=CIRCLE).addLabel("[$l^+$,$l^-$]", displace=0.4, angle=90)
vtx3 = Vertex(1., 1., mark=CIRCLE).addLabel("[$\\nu$]", displace=-0.4, angle=90)
vtx4 = Vertex(2.5, -1., mark=CIRCLE).addLabel("= [[$l^+$,$l^-$]]", displace=1.2)
vtx5 = Vertex(2.5, 1., mark=CIRCLE).addLabel("= [[$l^+$],[$\\nu$]]", displace=1.2)
#R-odd particles:
odd1 = Fermion(in_vtx, vtx1)
odd2 = Fermion(in_vtx, vtx2)
odd3 = Fermion(vtx1, vtx3)
odd4 = Fermion(vtx2, vtx4)
odd5 = Fermion(vtx3, vtx5)

#R-even outer points
out1 = Point(-0.25,  2)
out2 = Point(0.75, -2)
out3 = Point(1.25,  2)
out4 = Point(0.25, -2)
#R-even particles
even1 = Fermion(vtx1, out1)
even2 = Fermion(vtx2, out2)
even3 = Fermion(vtx3, out3)
even4 = Fermion(vtx2, out4)

fd.draw("C1N2labels.pdf")


#General Topology
nvertices = 6
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

odd = Higgs(vtxs[nvertices-2],Point(vtxs[nvertices-2].getX()+0.7,1.))
odd = Higgs(vtxs[nvertices-1],Point(vtxs[nvertices-2].getX()+0.7,-1.))
odd = Fermion(Point(vtxs[nvertices-2].getX()+1.2,1.),Vertex(vtxs[nvertices-2].getX()+1.7, 1., mark=CIRCLE))
odd = Fermion(Point(vtxs[nvertices-1].getX()+1.2,-1.),Vertex(vtxs[nvertices-1].getX()+1.7, -1., mark=CIRCLE))


#R-even outer points
outpts = []
for i in range(nvertices/2):
  out = Point(-0.25 + float(i)*1.5,2)
  outpts.append(out)
  out = Point(-0.25 + float(i)*1.5,-2)
  outpts.append(out)

#R-even particles
evens = []
random.seed(3.)
for i in range(0,nvertices,2): 
  disp = 0.2
  ran = random.uniform(-1.,1.)
  if abs(ran) < 0.5 and ran < 0.:
    even = Fermion(vtxs[i],Point(outpts[i].getX()-0.25,outpts[i].getY())).addLabel("$P_{"+str(len(evens)+1)+"}$", pos=0.7, displace=-disp)
    evens.append(even)
    even = Fermion(vtxs[i],outpts[i]).addLabel("$P_{"+str(len(evens)+1)+"}$", pos=0.7, displace=disp)
    evens.append(even)
  else:
    even = Fermion(vtxs[i],outpts[i]).addLabel("$P_{"+str(len(evens)+1)+"}$", pos=0.7, displace=disp)
    evens.append(even)

  if abs(ran) < 0.5 and ran > 0.:  
    even = Fermion(vtxs[i+1],Point(outpts[i+1].getX()-0.25,outpts[i+1].getY())).addLabel("$P_{"+str(len(evens)+1)+"}$", pos=0.7, displace=disp)
    evens.append(even)
    even = Fermion(vtxs[i+1],outpts[i+1]).addLabel("$P_{"+str(len(evens)+1)+"}$", pos=0.7, displace=-disp)
    evens.append(even)
  else:
    even = Fermion(vtxs[i+1],outpts[i+1]).addLabel("$P_{"+str(len(evens)+1)+"}$", pos=0.7, displace=-disp)
    evens.append(even)


fd.draw("GeneralTop.pdf")


