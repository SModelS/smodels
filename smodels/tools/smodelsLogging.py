"""
.. module:: smodelsLogging
   :synopsis: Simple code that creates and configures a central logger

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
    
from smodels.tools.colors import colors
import logging

class ColoredLogger:

    def __init__ ( self ):
        FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in' \
                 ' %(lineno)s: %(message)s'
        logging.basicConfig(format=FORMAT)
        self._logger = logging.getLogger("smodels")
    def setLevel ( self, level ):
        return self._logger.setLevel ( level )
    def set_level ( level ):
        return self._logger.setLevel ( level )
    def error ( self, msg, *args, **kwargs):
        msg = "%s%s%s" % ( colors.error, msg, colors.reset )
        return self._logger.error ( msg, *args, **kwargs )
    def warning ( self, msg, *args, **kwargs):
        msg = "%s%s%s" % ( colors.warn, msg, colors.reset )
        return self._logger.warning ( msg, *args, **kwargs )
    def info ( self, msg, *args, **kwargs):
        msg = "%s%s%s" % ( colors.info, msg, colors.reset )
        return self._logger.info ( msg, *args, **kwargs )
    def debug ( self, msg, *args, **kwargs):
        msg = "%s%s%s" % ( colors.debug, msg, colors.reset )
        return self._logger.debug ( msg, *args, **kwargs )

logger = ColoredLogger()

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
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if not level in levels:
        logger.error ( "Unknown log level ``%s'' supplied!" % level )
        return
    logger.setLevel ( level = levels[level] )
