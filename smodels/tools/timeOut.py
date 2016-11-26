"""
.. module:: timeOut
   :synopsis: Facility to implement a time out option when running smodels

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import signal
from smodels.tools.smodelsLogging import logger

class NoTime(Exception):
    """
    The time out exception. Raised when the running time exceeds timeout
    """
    def __init__(self, value=None):
        self.value = value
        Exception.__init__(self, value)
        
    def __str__(self):
        return '%.2es time out exceeded' %float(self.value)


class Timeout():
    """Timeout class using ALARM signal."""
    
    def __init__(self, sec):
        self.sec = sec
        if type ( sec ) != int:
            logger.warning ( "timeout set to a non-integral number of seconds."
                             " Will try to cast to integer." )
            self.sec = int ( sec )
 
    def __enter__(self):
        if self.sec:            
            signal.signal(signal.SIGALRM, self.raise_timeout)
            signal.alarm(self.sec)            
 
    def __exit__(self, *args):
        signal.alarm(0)    # disable alarm
 
    def raise_timeout(self, *args):
        raise NoTime(self.sec)
