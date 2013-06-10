#!/usr/bin/python

#with -t give input topology
#run runRecreate.py for all analyses


import os, argparse, set_path
from Experiment import SMSResults, SMSInterpolation

argparser=argparse.ArgumentParser()
argparser.add_argument('-t','--topo',help='input topology [T2bb]',default='T2bb')
argparser.add_argument('-mz','--mz',help='intermediate mass information')
argparser.add_argument('-axes','--axes',help='axes information')
argparser.add_argument('-n','--nevts',help='number of events per point in refXSec [10000]', type=int, default=10000)
argparser.add_argument('-b','--binsize',help='binsize in GeV', type=int)
argparser.add_argument('-text','--text',help='write numbers in 100 GeV distances',action='store_true')
args=argparser.parse_args()

anas = SMSResults.getAnalyses(args.topo)

ourAnas = ["alphaT", "alphaT8TeV", "Weakinos8TeV", "RA48TeV", "RA2b8TeV", "LeptonicStop8TeV", "MultiLepton8TeV", "SUS13008", "ATLAS_CONF_2013_024", "ATLAS_CONF_2013_037", "ATLAS_CONF_2013_007", "ATLAS_CONF_2013_001", "ATLAS_CONF_2012_105", "ATLAS_CONF_2012_166"]

toponame = args.topo
if args.mz: toponame = SMSInterpolation.gethistname(args.topo, args.mz)

for ana in anas:
  if ana not in ourAnas: continue
  options="-t %s -a %s" %(args.topo, ana)
  if args.mz: options=options+" -mz %s" % args.mz
  if args.axes: options=options+" -axes %s" % args.axes
  if args.nevts: options=options+" -n %d" % args.nevts
  if args.binsize: options=options+" -b %d" % args.binsize
  if args.text: options=options+" -text"
  os.system("./runRecreate.py %s" % options)
  if not args.binsize:
    plotname= "%s_%s_%devts.png" %(ana,toponame,args.nevts)
  else:
    plotname= "%s_%s_%devts_%dGeVbin.png" %(ana,toponame,args.nevts,args.binsize)
  os.system("cp ../plots/%s /var/www/gallery/%s" %(plotname, plotname))
