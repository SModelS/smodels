#!/usr/bin/env python

"""
.. module:: FeynmanGraphs
    :synopsis: This unit contains two simple routines that draw feynman graphs.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""


def printParticle_ ( label ):
  """ very simple method to rename a few particles for the asciidraw
      routine, do not call directly """
  if label=="jet": label="q"
  label=label+"   "
  return label[:2]
  
def drawBranch_ ( branch, upwards, labels, html ):
  """ draws a single branch, should only be used via .asciidraw, 
      not directly """
  lines=["   ","----"]
  labels="   "
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
  if html: print HTML
  if upwards and labels: print labels
  if html: print HTML
  for i in order: print lines[i]
  if html: print HTML
  if not upwards and labels: print labels
  if html: print HTML

def asciidraw ( element, labels=True, html=False ):
  """ draw a simple ascii graph on the screen """
  for (ct,branch) in enumerate(element.B):
    drawBranch_ ( branch, upwards=(ct==0), labels=labels, html=html )
