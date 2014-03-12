#!/usr/bin/python

""" another sandbox to try things out """

import set_path, tempfile
from theory import SLHATools
from tools import xsecComputer
from tools.VariousHelpers import logging
log = logging.getLogger(__name__)

slhafile="../slha/xsec.slha"


nevts=10000
Tmp=tempfile.mkdtemp()
log.info ( "now run pythia in "+Tmp )
Wv=xsecComputer.compute(nevts,slhafile,rpythia = True, datadir=Tmp)
log.info ( "done running pythia" )
#print "Wv=",Wv
print "1000025 1000035 8 tev:",Wv.getCrossSection ( 1000025, 1000035, order="NLL", sqrts=8 )

xsec=SLHATools.xSecFromSLHAFile ( slhafile )
print "1000025 1000035 8 tev:",xsec.getCrossSection ( 1000025, 1000035, order="NLL", sqrts=8 )

