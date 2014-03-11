"""
This package contains all code that cannot be classified as being part of
*experiment* or *theory*.
"""
import os
import logging.config

basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
logging.config.fileConfig('%s/logging.conf' % basepath, disable_existing_loggers=False)
logger = logging.getLogger(__name__)
