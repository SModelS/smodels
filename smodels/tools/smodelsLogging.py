"""
.. module:: smodelsLogging
   :synopsis: Simple code that creates and configures a central logger

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import logging

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger("smodels")

def setLogLevel ( level ):
    """ set the log level of the central logger. 
        can either be directly an integer ( e.g. logging.DEBUG ),
        or "debug", "info", "warning", or "error".
    """
    if level == None: return
    if type ( level ) == int:
        logger.setLevel ( level=level )
        return
    level = level.lower()
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warning": logging.WARNING, "error": logging.ERROR }
    if not level in levels:
        logger.error ( "Unknown log level ``%s'' supplied!" % level )
        return
    logger.setLevel ( level= levels[level] )
