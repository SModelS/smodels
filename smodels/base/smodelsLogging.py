"""
.. module:: smodelsLogging
   :synopsis: Simple code that creates and configures a central logger

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import logging


class Colors:
    def __init__ ( self ):
        self.on = False

    @property
    def debug ( self ):
        return ""

    @property
    def error ( self ):
        return self.red

    @property
    def warn ( self ):
        return self.yellow

    @property
    def info ( self ):
        return self.green

    @property
    def magenta( self ):
        if self.on == False: return ""
        return '\033[0;35m'

    @property
    def green( self ):
        if self.on == False: return ""
        return '\033[0;32m'

    @property
    def red ( self ):
        if self.on == False: return ""
        return '\033[0;31m'

    @property
    def yellow ( self ):
        if self.on == False: return ""
        return '\033[0;33m'

    @property
    def cyan ( self ):
        if self.on == False: return ""
        return '\033[0;36m'

    @property
    def blue ( self ):
        if self.on == False: return ""
        return '\033[0;34m'

    @property
    def reset ( self ):
        if self.on == False: return ""
        return '\033[;0m'

colors = Colors()


class ColorizedStreamHandler(logging.StreamHandler):
    def _color_wrap(self, *c):
        def wrapped(inp):
            return "".join(list(c) + [inp, colors.reset])
        return wrapped

    def __init__(self, stream=None):
        logging.StreamHandler.__init__(self, stream)

    def should_color(self):

        # If the stream is a tty we should color it
        if hasattr( self.stream, "isatty") and self.stream.isatty():
            return True

        # If we have an ASNI term we should color it
        if os.environ.get("TERM") == "ANSI":
            return True

        # If anything else we should not color it
        return False

    def format(self, record):
        msg = logging.StreamHandler.format(self, record)

        if self.should_color():
            COLORS = [
                # This needs to be in order from highest logging level to lowest.
                (logging.ERROR, self._color_wrap(colors.error)),
                (logging.WARNING, self._color_wrap(colors.warn)),
                (logging.INFO, self._color_wrap(colors.info)),
            ]
            for level, color in COLORS:
                if record.levelno >= level:
                    msg = color(msg)
                    break

        return msg

def getLogger ():
    FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in' \
       ' %(lineno)s: %(message)s'
    logging.basicConfig(format=FORMAT)
    formatter = logging.Formatter( FORMAT )
    ch = ColorizedStreamHandler()
    ch.setFormatter ( formatter )
    logger = logging.getLogger("smodels")
    logger.addHandler(ch)
    logger.propagate = False
    return logger

logger = getLogger()

def getLogLevel( asString=False ):
    """ obtain the current log level.
    :params asString: return string, not number.
    """
    # return logger.level
    ret = logger.getEffectiveLevel()
    if not asString:
        return ret
    lvlNames= { logging.DEBUG: "debug", logging.INFO: "info",
                logging.WARNING: "warning", logging.ERROR: "error",
                logging.CRITICAL: "critical", logging.FATAL: "fatal" }
    lvl = list ( lvlNames.keys() )
    lvl.sort( reverse = True)
    for l in lvl:
        if ret >= l:
            return lvlNames[l]
    return "unknown log level"


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
               "warn": logging.WARNING, "warning": logging.WARNING,
               "error": logging.ERROR, "critical": logging.CRITICAL,
               "fatal": logging.FATAL }
    if not level in levels:
        logger.error ( f"Unknown log level ``{level}'' supplied!" )
        return
    logger.setLevel ( level = levels[level] )
