
"""
.. module:: printerRegistry
   :synopsis: A facility where all printers should register,
   for the masterPrinter to access them

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from smodels.tools.printers.basicPrinter import BasicPrinter
from typing import Optional

class PrinterRegistry:
    """
    Facility for all printers to register
    """
    printers = {}

    @classmethod
    def register(cls, printer : BasicPrinter, extension : str ) -> bool:
        """
        register this printer
        :param printer: Printer object to register
        :param extension: Extension this printer will be in charge of
        :returns: False if printer already exists (or other error),
        else True
        """
        if extension in cls.printers:
            return False
        cls.printers[extension] = printer
        return True

    @classmethod
    def has(cls, extension : str ) -> bool:
        """ check if a printer is registered for given extensions
        """
        return extension in cls.printers

    @classmethod
    def get(cls, extension : str ) -> Optional[BasicPrinter]:
        """ get the printer that serves a specific extension
        :returns: None, if no printer registered for that extension
        """
        if not extension in cls.printers:
            return None
        return cls.printers[extension]
