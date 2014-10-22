"""
.. module:: tools.__init__
    :synopsis: This package contains tools for handling results obtained with the
    main SModelS code.
"""

import os
import logging.config

# print "tools.__init__: basepath=",os.path.abspath(os.path.join(os.path.dirname(__file__)))
basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
# print "tools.__init__: basepath=",basepath
logging.config.fileConfig('%s/etc/logging.conf' % basepath,
                          disable_existing_loggers=False)
logger = logging.getLogger(__name__)
