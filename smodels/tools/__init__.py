"""
.. module:: tools.__init__
    :synopsis: This package contains tools for handling results obtained with the
       main SModelS code.
"""

import os
import logging.config

basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
logging.config.fileConfig('%s/etc/logging.conf' % basepath,
                          disable_existing_loggers=False)
logger = logging.getLogger(__name__)
