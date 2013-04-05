#!/usr/bin/env python


""" this little script SMS-decomposes an LHE file """

import SMSDecomposition, os, sys, CountCascades, types, math
import lheSimplifier

def getLSP ( particles, name=False ):
    """ obtain the pdgid of the lsp """
    mylsp=1 # look for the lsp
    m=99999999
    for p in particles:
      if p.pdgid > 1000000 and abs(p.M) < abs(m):
        mylsp=p.pdgid
        m=abs(p.M)
    if name: return SMSDecomposition.particle ( mylsp )
    return mylsp

def readOneEvent ( reader, debug, countcascades=True, simplify=False ):
  """ read and discuss a single event """
  if debug:
    print "[lheDecomposer.py] next event!"
  import numpy, ROOT
  vec=ROOT.TLorentzVector()
  nn=reader.getN()
  particles=[ reader.get(i) for i in range(nn) ]
  if simplify: particles=lheSimplifier.simplify(particles, debug)
  lsp=getLSP ( particles, True )

  n={}
  for i in SMSDecomposition.countingvariables: n[i]=0
  for p in particles:
    pname=SMSDecomposition.particle (p.pdgid)
    if n.has_key ( pname ): n[pname]+=1
    if int(p.mother)==0:
      n["p"+pname]+=1
  ret={}
  ret["sms"]=SMSDecomposition.getSMS ( n, lsp )
  ret["prodmode"]=SMSDecomposition.prodmode ( n )
  if debug:
    s=SMSDecomposition.describe ( n )
    print "[lheDecomposer.py] sms=",ret["sms"],"s=",s
  c0,c1=0,0
  if countcascades:
    mylsp=getLSP ( particles )
    if mylsp > 1000000: 
      if debug: 
        print "[lheDecomposer.py] now count cascades, lsp=",mylsp
      chain=CountCascades.reconstructChains ( particles, lsp=mylsp )
    if len(chain)>1:
      c0=len(chain[0])
      c1=len(chain[1])
  ret["c0"]=c0
  ret["c1"]=c1
  if debug:
    print "[lheDecomposer.py] returning:",ret
  return ret


def inc ( dic, key ):
  if not dic.has_key ( key ): dic[key]=0
  dic[key]+=1

def createReader ( filename, reader="LHEReader_cc.so" ):
  ## get the C++ LHE reader
  import ROOT
  ROOT.gROOT.SetBatch()
  ROOT.gSystem.Load(reader)
  return ROOT.LHEReader( filename )

def read ( filename, nmax=-1, debug=False, readerso="LHEReader_cc.so", simplify=False, countcascades=True ):
  """ read an lhe file name <filename>, but max nmax events. Returns dictionary with SMS counts."""
  if debug:
    print "[lheDecomposer.py] start reading",filename
  Reader=createReader ( filename, readerso )
  count={"base":0,"total":0 }
  ctr=0
  while Reader.next():
    ctr+=1
    if nmax>0 and ctr>nmax:
      break
    tmp=readOneEvent ( Reader, debug, countcascades=countcascades, simplify=simplify )
    if countcascades:
      inc(count,"c_%d%d"  % (tmp["c0"],tmp["c1"]) )
    inc(count,tmp["sms"])
    inc(count,tmp["prodmode"])
    inc(count,"base")
    inc(count,"total")
  Reader.close()
  return count

def readSequential ( filename, nmax=-1, readerso="LHEReader_cc.so", debug=False ):
  """ reads an LHE file, just like the read method. But returns an array of the SMSes. """
  ret=[]
  ctr=0
  Reader=createReader ( filename, readerso )
  while Reader.next():
    ctr+=1
    if nmax>0 and ctr>nmax:
      break
    ret.append ( readOneEvent ( Reader, debug=debug, countcascades=False, simplify=False )["sms"] )
  return ret


def translateCountsIntoXSecs ( totalxsec, counts ):
  """ given the total xsec and counts, 
      we translate into xsec * brs """
  xsecs={}
  if totalxsec==None: 
    for i in counts: xsecs[i]=None
    return xsecs
  if not counts.has_key("base"):
    print "[lheDecomposer.py:translateIntoXSecs] error: cannot find ``base'' in counts=",counts
    sys.exit(0)
  base=float(counts["base"])
  if base and base!=0.:
    for (key,value) in counts.items():
      if not value: value=0.
      xsecs[key]=value
      if value.__class__==types.FloatType:
        xsecs[key]=float(value)/base * totalxsec
  return xsecs

if __name__ == "__main__":
  import argparse, types
  argparser = argparse.ArgumentParser(description='SMS-decompose an LHE source')
  argparser.add_argument ( '-f', '-l', '--lhe', nargs='?', help='LHE input file', type=types.StringType, default="1420_60_10_0_1.lhe" )
  argparser.add_argument ( '-r', '--reader', nargs='?', help='so file of LHEReader', type=types.StringType, default="./LHEReader_cc.so" )
  argparser.add_argument ( '-d', '--debug', help='turn on debugging', action='store_true' )
  argparser.add_argument ( '-s', '--simplify', help='perform simplification step', action='store_true' )
  argparser.add_argument ( '-c', '--countcascades', help='count cascades', action='store_true' )
  argparser.add_argument ( '-n', '--nmax', nargs='?', help='maximum number of events', type=types.IntType, default=-1 )
  args=argparser.parse_args()

  count=read ( args.lhe, args.nmax, args.debug, args.reader, countcascades=args.countcascades, simplify=args.simplify )
  print count
