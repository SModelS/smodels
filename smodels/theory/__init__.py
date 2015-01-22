"""
.. module:: theory.__init__
   :synopsis: This Package is intended to contain everything related to theory:   
      * cross section calculation code
      * sms decomposition code (LHE-based, SLHA-based)
      * some more tools, e.g. for reading/writing slha files, or particle names
   
"""

import os
import logging.config

basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
#if basepath.find(".egg")>0:
#    basepath=basepath[:-7] ## remove /smodels, if we're in an egg
## print "basepath=",basepath
logging.config.fileConfig('%s/etc/logging.conf' % basepath,
                          disable_existing_loggers=False)
logger = logging.getLogger(__name__)
