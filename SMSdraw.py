def draw ( top, filename="bla.pdf", elementnr=0 ):
  """ plot a lessagraph, write into pdf file called <filename> """
  from pyfeyn.user import FeynDiagram,Vertex,Point,Fermion,Scalar,CIRCLE
  fd = FeynDiagram()

  nbranches=len(top.B)
  # print nbranches,"branches"

  for (ct,branch) in enumerate(top.B):
    # print "branch",ct,branch,"with",branch.vertnumb,"vertices"
    p1 = Point(0, ct)
    lastVertex=p1
    for ( nvtx,particles) in enumerate(branch.vertparts):
      v1=Vertex ( nvtx+1,ct,mark=CIRCLE)
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

def asciidraw ( top, elementnr=0, labels=True ):
  """ draw a simple ascii graph on the screen """

  nbranches=len(top.B)
  # print nbranches,"branches"

  def drawBranch ( branch, up, labels ):
    """ draws a single branch """
    # print "branch",branch,"up=",up
    lines=["  ","---"]
    labels="  "
    idx=0
    for ( nvtx,particles) in enumerate(branch.vertparts):
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

  for (ct,branch) in enumerate(top.B):
    up=False
    if ct==0: up=True 
    drawBranch ( branch, up=up, labels=labels )
    print 
