"""
This package contains all code that cannot be classified as being part of
*experiment* or *theory*.

"""
import os
import sys
import logging.config

basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

if sys.version_info < (2, 5):
    logging.config.fileConfig('%s/etc/logging.conf' % basepath)
else:
    logging.config.fileConfig('%s/etc/logging.conf' % basepath, disable_existing_loggers=False)

logger = logging.getLogger(__name__)
