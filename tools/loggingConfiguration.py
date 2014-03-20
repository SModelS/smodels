#!/usr/bin/python

"""
.. module:: loggingConfiguration
   :synopsis: Creates a logger object that is configured according to
   <current-working-directory>/etc/logging.conf. If that file is not found,
   revert to <smodels-installation>/etc/logging.conf.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import logging.config

def configure():
    basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    logfile= '%s/etc/logging.conf' % basepath
    if not os.path.exists(logfile):
        # print "[tools/__init__.py] %s does not exist." % logfile
        import inspect
        base=os.path.realpath ( inspect.getabsfile(configure) )
        base=base.replace("tools/loggingConfiguration.py","")
        logfile=base+"etc/logging.conf"
        # print "base=",base
    logging.config.fileConfig(logfile, disable_existing_loggers=False)
    return logfile


def createLogger(name):
    configure()
    return logging.getLogger(name)


if __name__ == "__main__": 
    """
    If this module is called as a script, it shows the path to the log file.
    
    """
    print "Using log file", configure()
