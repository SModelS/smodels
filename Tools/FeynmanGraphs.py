""" This unit contains two simple routines that draw feynman graphs """
#from pyfeyn.user import FeynDiagram,Vertex,Point,Fermion,Scalar,CIRCLE,SQUARE,\
#  HATCHED135,Circle,pyx
from pyfeyn.user import *
import numpy, math

def printParticle_ ( label ):
  """ very simple method to rename a few particles for the asciidraw
      routine, do not call directly """
  if label=="jet": label="q"
  label=label+"   "
  return label[:3]
  
def segment ( p1, p2, spin, Bend=None ):
  l=NamedLine[spin](p1,p2)# 
  if Bend: l.bend(Bend)
  return l

def zero ():
  """ a super simple convenience thing to mark the (0,0) coordinates """
  c=Vertex(0.,0., mark=CIRCLE, fill=[ RED ] ) ## , radius=.01)
  c1=Vertex(1.0,0., mark=CIRCLE, fill=[ BLUE ] ) ## , radius=.01)
  c1=Vertex(0.,1., mark=CIRCLE, fill=[ BLUE ] ) ## , radius=.01)

def connect ( p1, p2, straight=True, label=None, spin="fermion", bend=True, \
              verbose=False, nspec=None, displace=None ):
  """ simple: draw a line from p1 to p2 
      :param p1: starting point
      :param p2: end point
      :param straight: straight lines or xkcd style?
      :param label: add a label?
      :param nspec: specify the number of segments, None is draw the number randomly
      :param displace: displace at fixed distance?
     
      :returns: array of all line segments
  """
  if straight:
    fl=NamedLine[spin](p1,p2)
    if displace==None: displace=.05
    if label: fl.addLabel ( label, pos=0.1, displace=displace )
    return fl

  fl=Fermion(p1,p2)
  fl.setStyles( [ WHITE ] )
  n=nspec
  if n==None:
    n=int ( math.floor( numpy.random.uniform(1.5,3.75) ) )
  points = [ p1 ]
  f=.0001
  for i in range (1,n):
    pt=fl.fracpoint( float(i) / float(n) )
    pt.setX ( pt.x() + numpy.random.normal( 0, f) )
    pt.setY ( pt.y() + numpy.random.normal( 0, f) )
    points.append ( pt )
  points.append ( p2 )
  b=.015
  a=numpy.random.normal ( 0, 1 )
  if a<0.: b=-b
  segs=[]
  if verbose: print "[FeynmanGraphs.py] ----- "
  for i in range(n):
    br=b * (-1)**i
    if not bend: br=None
    segs.append ( segment(points[i],points[i+1],spin, Bend=br ) )
    if verbose:
      print "[FeynmanGraphs.py] draw line from (%f,%f) to (%f,%f)" % ( points[i].x(), points[i].y(), points[i+1].x(), points[i+1].y() )
  if displace==None: displace=-.08
  if label: segs[-1].addLabel ( label, pos=0.7, displace=displace )
  return segs

def draw ( element, filename="bla.pdf", straight=False ):
  """ plot a lessagraph, write into pdf/eps/png file called <filename> """
  from pyx import text
  #import os
  text.set(mode="latex")
  #text.set(fontmaps="ttfonts.map" ) 
  # text.preamble(r"\usepackage[T1]{fontenc}")
  #text.preamble(r"\font\ttfgeorgia georgia at10pt")
  text.preamble(r"\usepackage{times}")
  fd = FeynDiagram()
  f=1.0

  in1  = Point(-1*f, -.75*f)
  in2  = Point(-1*f, 1.75*f)
  vtx1 = Circle(0,.5*f, radius=0.3*f).setFillStyle(HATCHED135)
  #P1a = Fermion(in1, vtx1 ).addLabel("\\ttfgeorgia P$_1$")
  # P1a = Fermion(in1, vtx1 ).addLabel("P$_1$")
  P1a = connect ( vtx1, in1, straight=straight, label="P$_1$", displace=.3 )
  for i in P1a:
    a1=i.addParallelArrow( pos=.44,displace=.0003,length=pyx.unit.length(1.60*f / float(len(P1a))), size=.0001)
    a2=i.addParallelArrow( pos=.44,displace=-.0003,length=pyx.unit.length(1.60*f / float(len(P1a))), size=.0001)
  #  i.setStyles ( [ BLUE ] )
  # P1a.addParallelArrow( pos=.44,displace=.0003,length=pyx.unit.length(1.75*f), size=.0001)
  # P1a.addParallelArrow( pos=.44,displace=-.0003,length=pyx.unit.length(1.75*f), size=.0001)
  # P2a = Fermion(in2, vtx1 ).addLabel("P$_2$",displace=.3)
  P2a = connect ( vtx1, in2, straight=straight, label="P$_2$", displace=.3 )
  for i in P2a:
    a1=i.addParallelArrow( pos=.44,displace=.0003,length=pyx.unit.length(1.60*f / float(len(P2a))), size=.0001)
    a2=i.addParallelArrow( pos=.44,displace=-.0003,length=pyx.unit.length(1.60*f / float(len(P2a))), size=.0001)
  #P2a.addParallelArrow( pos=.44,displace=.0003,length=pyx.unit.length(1.75*f), size=.0001)
  #P2a.addParallelArrow( pos=.44,displace=-.0003,length=pyx.unit.length(1.75*f), size=.0001) 

  # nbranches=len(element.B)

  for (ct,branch) in enumerate(element.B):
    # print "branch",ct,branch,"with",branch.vertnumb,"vertices"
    # p1 = Point(0, ct)
    lastVertex=vtx1
    for ( nvtx,insertions) in enumerate(branch.particles):
      mark=None
      if len(insertions)>0: mark=CIRCLE
      v1=Vertex ( f*(nvtx+1),f*ct,mark=mark)
      # f1 = Scalar  ( lastVertex,v1) ## .addLabel ( "x")
      f1 = connect ( lastVertex,v1, straight=straight, spin="scalar", bend=True, verbose=False, nspec=3 )
      if nvtx==0:
        b=.10
        if ct==1: b=-.1
        #for xf in f1: xf.bend(b)
      lastVertex=v1
      # print "particles",particles,"ct=",ct
      y=-1.0*f ## y of end point of SM particle
      if ct==1: y=2.*f
      dx=(len(insertions)-1)*(-.5)*f ## delta_x 
      #dx=(particles-1)*(-.5)*f ## delta_x 
      for (i,insertion) in enumerate(insertions):
        p=Point ( f*(nvtx + 1 +  dx + i), f*y )
        ## print "branch=",branch
        label=printParticle_ ( insertion )
        ## ff=Fermion(v1,p).addLabel ( label )
        connect ( v1, p, straight=straight, label=label )
         
    pl = Point ( nvtx+2,ct )
    # fl = Scalar ( lastVertex,pl ) ## .addLabel( "lsp" )
    connect ( lastVertex,pl, straight=straight, spin="scalar" )
    

  # zero()
  pdffile=filename.replace("png","pdf")
  fd.draw( pdffile )
  if pdffile!=filename:
    import os
    os.system ( "convert %s %s" % ( pdffile, filename ) )

def drawBranch_ ( branch, upwards, labels ):
  """ draws a single branch, should only be used via asciidraw, 
      not directly """
  lines=["  ","---"]
  labels="  "
  ## for ( nvtx,particles) in enumerate(branch.particles):
  for insertions in branch.particles:
    if len(insertions)==0: continue
    lines[1]+="*---"
    if len(insertions)==1: 
      labels+=" "+printParticle_(insertions[0])
      lines[0]+=" |  "
    if len(insertions)==2: 
      labels+=printParticle_(insertions[0])+" "+printParticle_(insertions[1])
      if upwards:
        lines[0]+="\ /"
      else:
        lines[0]+="/ \\"
    if len(insertions)>2:
      print "[SMSFeynmanGraphs.py] case for n-body decay, n>3 not yet. implemented. Please implement."
      sys.exit(0)

  order=[0,1]
  if not upwards: order=[1,0]
  if upwards and labels: print labels
  for i in order: print lines[i]
  if not upwards and labels: print labels

def asciidraw ( element, labels=True ):
  """ draw a simple ascii graph on the screen """
  for (ct,branch) in enumerate(element.B):
    drawBranch_ ( branch, upwards=(ct==0), labels=labels )
    print 
