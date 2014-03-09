"""
This package is intended to contain everything that has to do with the
experimental results.
"""

import SMSResults
import SMSAnalysisFactory
import LimitGetter
import TxNames
import os
import logging.config

basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
logging.config.fileConfig('%s/logging.conf' % basepath, disable_existing_loggers=False)
logger = logging.getLogger(__name__)


