"""
.. module:: basicPrinter
   :synopsis: Base class for defining printer classes

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import sys
import os
from smodels.base.smodelsLogging import logger
from smodels.statistics.basicStats import observed, apriori, aposteriori
import numpy as np
import time


class BasicPrinter(object):
    """
    Super class to handle the basic printing methods
    """

    def __init__(self, output, filename, outputFormat = 'current'):
        """
        :ivar str typeofexpectedvalues: what type of evaluationType values to print,
              apriori or posteriori
        """

        self.name = "basic"
        self.time = time.time()  # time stamps

        self.outputList = []
        self.filename = filename
        self.output = output
        self.printingOrder = []
        self.typeofexpectedvalues = apriori
        self.toPrint = []
        self.outputFormat = outputFormat

        if filename and os.path.isfile(filename):
            logger.warning( f"Removing file {filename}" )
            os.remove(filename)

    def getTypeOfExpected(self):
        """ tiny convenience function for what evaluationType values to print,
            apriori (True) or posteriori """
        evaluationType = apriori
        if self.typeofexpectedvalues in [ "aposteriori", "posteriori", aposteriori ]:
            evaluationType = aposteriori
        return evaluationType

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, fn):
        self._filename = fn
        self.mkdir()

    def mkdir(self):
        """ create directory to file, if necessary """
        if not self.filename:
            return
        dirname = os.path.dirname(self.filename)
        if dirname != "" and not os.path.exists(dirname):
            os.makedirs(dirname)

    def setOptions(self, options):
        """
        Store the printer specific options to control the output of each printer.
        Each option is stored as a printer attribute.

        :param options: a list of (option,value) for the printer.
        """

        for opt, value in options:
            setattr(self, opt, eval(value))

    def addObj(self, obj):
        """
        Adds object to the Printer.

        :param obj: A object to be printed. Must match one of the types defined in formatObj

        :return: True if the object has been added to the output. If the object does not belong
                to the pre-defined printing list toPrint, returns False.
        """

        for iobj, objType in enumerate(self.printingOrder):
            if isinstance(obj, objType):
                self.toPrint[iobj] = obj
                return True
        return False

    def openOutFile(self, filename, mode):
        """ creates and opens a data sink,
            creates path if needed """
        d = os.path.dirname(filename)
        if not os.path.exists(d):
            os.makedirs(d)
            logger.info(f"creating directory {d}")
        return open(filename, mode)

    def flush(self):
        """
        Format the objects added to the output, print them to the screen
        or file and remove them from the printer.
        """
        ret = ""

        for obj in self.toPrint:
            if obj is None:
                continue
            output = self._formatObj(obj)
            if not output:
                continue  # Skip empty output
            ret += output
            if self.output == 'stdout':
                sys.stdout.write(output)
            elif self.output == 'file':
                if not self.filename:
                    logger.error('Filename not defined for printer')
                    return False
                with self.openOutFile(self.filename, "a") as outfile:
                    outfile.write(output)
                    outfile.close()

        self.toPrint = [None]*len(self.printingOrder)  # Reset printing objects
        self.time = time.time()  # prepare next timestamp
        return ret

    def _formatObj(self, obj):
        """
        Method for formatting the output depending on the type of object
        and output.

        :param obj: A object to be printed. Must match one of the types defined in formatObj

        """

        typeStr = type(obj).__name__
        try:
            formatFunction = getattr(self, '_format'+typeStr)
            ret = formatFunction(obj)
            # print ( " `-", len(ret))
            return ret
        except AttributeError as e:
            logger.warning(f'Error formating object {typeStr}: \n {e}')
            return False

    def _round(self, number, n=6):
        """ round a number to n significant digits, if it *is* a number """
        if type(number) not in [float, np.float64]:
            return number
        if not np.isfinite(number):
            return f'float("{number}")'
        if np.isnan(number) or not np.isfinite(number):
            return number
        try:
            if abs(number) < 1e-40:
                return number
            return round(number, -int(np.floor(np.sign(number) * np.log10(abs(number)))) + n)
        except Exception:
            pass
        return number
        # return round ( number, n )
