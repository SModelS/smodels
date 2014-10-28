"""
.. module:: experiment.__init__
   :synopsis: This package is intended to contain everything that has to do
              with the experimental results.


"""
import os
import logging.config

basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
#if basepath.find(".egg")>0:
#    basepath=basepath[:-7] ## remove /smodels, if we're in an egg
logging.config.fileConfig('%s/etc/logging.conf' % basepath,
                          disable_existing_loggers=False)
logger = logging.getLogger(__name__)
