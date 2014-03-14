#!/usr/bin/env python

"""
.. module:: asciiGraphs
    :synopsis: This unit contains a simple routine to draw ASCII-art Feynman graphs.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def printParticle_ ( label ):
  """ very simple method to rename a few particles for the asciidraw
      routine, do not call directly """
  if label=="jet": label="q"
  label=label+"   "
  return label[:2]

def drawBranch_ ( branch, upwards, labels, html, border, L ):
  """ draws a single branch, should only be used via .asciidraw,
      not directly """
  length=0
  lines=["   ","----"]
  labels="   "
  if border and upwards:
    lines=[" |    "," | ----"]
    labels=" |    "
  if border and not upwards:
    lines=[" |    "," | ----"]
    labels=" |    "

  for insertions in branch.particles:
    if len(insertions)==0:
      lines[0]+=" "
      lines[1]+="*"
      continue
    lines[1]+="*----"
    if len(insertions)==1:
      labels+=" "+printParticle_(insertions[0])+"  "
      lines[0]+=" |   "
    if len(insertions)==2:
      labels+=printParticle_(insertions[0])+" "+printParticle_(insertions[1])
      if upwards:
        lines[0]+="\\ /  "
      else:
        lines[0]+="/ \\  "
    if len(insertions)>2:
      print "[SMSFeynmanGraphs.py] case for n-body decay, n>3 not yet. implemented. Please implement."
      sys.exit(0)

  order=[0,1]
  if not upwards: order=[1,0]
  HTML="<br>"
  lengthdiff=L-len(lines[0])/5
  if border:
    if L==2:
      lines[0]+=" "
      lines[1]+=" "
      labels+=" "
    labels+=" "+" "*(5*lengthdiff)+" |"
    lines[0]+=" "*(5*lengthdiff+0)+"  |"
    lines[1]+=" "*(5*lengthdiff+0)+" |"
  if border and upwards: print " /"+"-"*(4*L+4)+"\\"
  if html: print HTML
  if upwards and labels: print labels
  if html: print HTML
  for i in order: print lines[i]
  if html: print HTML
  if not upwards and labels: print labels
  if html: print HTML
  if border and not upwards: print " \\"+"-"*(4*L+4)+"/"

def asciidraw ( element, labels=True, html=False, border=False ):
  """ draw a simple ascii graph on the screen """
  L=[]
  for (ct,branch) in enumerate(element.B):
    L.append ( int( str(branch).count("[") ) )
  for (ct,branch) in enumerate(element.B):
    drawBranch_ ( branch, upwards=(ct==0), labels=labels, html=html, border=border, L=max(L) )

if __name__ == "__main__":
    import set_path, argparse, types
    import SModelS

    argparser = argparse.ArgumentParser(description='simple tool that is meant to draw lessagraphs, as an ascii plot') 
    argparser.add_argument ( '-T', nargs='?', help='Tx name, will look up lhe file in ../regression/Tx_1.lhe. Will be overriden by the "--lhe" argument', type=types.StringType, default='T1' )
    argparser.add_argument ( '-l', '--lhe', nargs='?', help='lhe file name, supplied directly. Takes precedence over "-T" argument.', type=types.StringType, default='' )
    argparser.add_argument ( '-b', '--border', help='draw a border around the graph', action='store_true' )
    args=argparser.parse_args()

    from theory import LHEReader, lheDecomposer, crossSection

    filename="%s/lhe/%s_1.lhe" % (SModelS.installDirectory(), args.T )
    if args.lhe!="": filename=args.lhe

    reader = LHEReader.LHEReader( filename )
    Event = reader.next()
    SMSTop = lheDecomposer.elementFromEvent( Event, crossSection.XSectionList() )
    asciidraw ( SMSTop[0].leadingElement(), border=args.border )
            
