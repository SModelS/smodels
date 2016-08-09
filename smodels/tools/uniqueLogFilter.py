"""
.. module:: uniqueLogFilter
   :synopsis: Contains a stolen code snippet for a logging filter to
       have identical log messages appear only once.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import logging

class UniqueFilter(logging.Filter):
    """Messages are allowed through just once.
    The 'message' includes substitutions, but is not formatted by the 
    handler. If it were, then practically all messages would be unique!
    stolen from: http://code.activestate.com/recipes/412552-using-the-logging-module/
    """
    def __init__(self, name=""):
        logging.Filter.__init__(self, name)
        self.reset()
    def reset(self):
        """Act as if nothing has happened."""
        self.__logged = {}
    def filter(self, rec):
        """logging.Filter.filter performs an extra filter on the name."""
        return logging.Filter.filter(self, rec) and self.__is_first_time(rec)
    def __is_first_time(self, rec):
        """Emit a message only once."""
        msg = rec.msg %(rec.args)
        if msg in self.__logged:
            self.__logged[msg] += 1
            return False
        else:
            self.__logged[msg] = 1
            return True

