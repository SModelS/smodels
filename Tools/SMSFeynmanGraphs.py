""" This unit contains two simple routines that draw feynman graphs """

def draw ( top, filename="bla.pdf", elementnr=0 ):
  """ plot a lessagraph, write into pdf file called <filename> """
  from pyfeyn.user import FeynDiagram,Vertex,Point,Fermion,Scalar,CIRCLE,SQUARE,\
    HATCHED135,Circle,pyx
  fd = FeynDiagram()
  f=1.0

  in1  = Point(-1*f, -.75*f)
  in2  = Point(-1*f, 1.75*f)
  vtx1 = Circle(0,.5*f, radius=0.3*f).setFillStyle(HATCHED135)
  P1a = Fermion(in1, vtx1 ).addLabel("P$_1$")
  P1a.addParallelArrow( pos=.44,displace=.0003,length=pyx.unit.length(1.75*f), size=.0001)
  P1a.addParallelArrow( pos=.44,displace=-.0003,length=pyx.unit.length(1.75*f), size=.0001)
  P2a = Fermion(in2, vtx1 ).addLabel("P$_2$",displace=.3)
  P2a.addParallelArrow( pos=.44,displace=.0003,length=pyx.unit.length(1.75*f), size=.0001)
  P2a.addParallelArrow( pos=.44,displace=-.0003,length=pyx.unit.length(1.75*f), size=.0001) 

  nbranches=len(top.leadingElement().B)
  # print nbranches,"branches"

  for (ct,branch) in enumerate(top.leadingElement().B):
    # print "branch",ct,branch,"with",branch.vertnumb,"vertices"
    # p1 = Point(0, ct)
    lastVertex=vtx1
    for ( nvtx,particles) in enumerate(branch.particles):
      mark=None
      if particles>0: mark=CIRCLE
      v1=Vertex ( f*(nvtx+1),f*ct,mark=mark)
      f1 = Scalar  ( lastVertex,v1) ## .addLabel ( "x")
      if nvtx==0:
        b=.25
        if ct==1: b=-.25
        f1.bend(b)
      lastVertex=v1
      # print "particles",particles,"ct=",ct
      y=-1.0*f ## y of end point of SM particle
      if ct==1: y=2.*f
      dx=(len(particles)-1)*(-.5)*f ## delta_x 
      #dx=(particles-1)*(-.5)*f ## delta_x 
      for i in range(len(particles)):
        p=Point ( f*(nvtx + 1 +  dx + i), f*y )
        print "branch=",branch
        label=particles[i].replace("jet","q")
        #label=branch.ElList[elementnr].particles[i].replace("jet","q")
        ff=Fermion(v1,p).addLabel ( label )
         
    #pl = Point ( nvtx+2,ct )
    #fl = Scalar ( lastVertex,pl ) ## .addLabel( "lsp" )

  fd.draw( filename )

def asciidraw ( top, elementnr=0, labels=True ):
  """ draw a simple ascii graph on the screen """

  nbranches=len(top.leadingElement().B)
  # print nbranches,"branches"

  def drawBranch ( branch, up, labels ):
    """ draws a single branch """
    print "branch",branch,"up=",up
    lines=["  ","---"]
    labels="  "
    idx=0
    for ( nvtx,particles) in enumerate(branch.particles):
      if particles==0: continue
      lines[1]+="*---"
      if particles==1: 
        labels+=" "+branch.ElList[elementnr].particles[idx][0]
	idx+=1
        lines[0]+=" | "
      if particles==2: 
        labels+=branch.ElList[elementnr].particles[idx][0]+" "+branch.ElList[elementnr].particles[idx+1][0]
	idx+=2
        if up:
          lines[0]+="\ /"
        else:
          lines[0]+="/ \\"

    order=[0,1]
    if not up: order=[1,0]
    if up and labels: print labels
    for i in order: print lines[i]
    if not up and labels: print labels

  for (ct,branch) in enumerate(top.leadingElement().B):
    up=False
    if ct==0: up=True 
    drawBranch ( branch, up=up, labels=labels )
    print 

def simpledraw ( top, filename="bla.pdf", elementnr=0 ):
  """ plot a simple lessagraph, write into pdf file called <filename> """
  from pyfeyn.user import FeynDiagram,Vertex,Point,Fermion,Scalar,CIRCLE
  fd = FeynDiagram()

  nbranches=len(top.B)
  # print nbranches,"branches"

  for (ct,branch) in enumerate(top.B):
    # print "branch",ct,branch,"with",branch.vertnumb,"vertices"
    p1 = Point(0, ct)
    lastVertex=p1
    for ( nvtx,particles) in enumerate(branch.vertparts):
      v1=Vertex ( nvtx+1,ct,mark=None )
      f1 = Fermion  ( lastVertex,v1) ## .addLabel ( "x")
      lastVertex=v1
      # print "particles",particles,"ct=",ct
      y=-1 ## y of end point of SM particle
      if ct==1: y=2
      dx=(particles-1)*(-.5) ## delta_x 
      for i in range(particles):
        p=Point ( nvtx + 1 +  dx + i, y )
        f=Scalar(v1,p).addLabel ( branch.ElList[elementnr].particles[i] )
         
    pl = Point ( nvtx+2,ct )
    fl = Fermion ( lastVertex,pl ) ## .addLabel( "lsp" )

  fd.draw( filename )
