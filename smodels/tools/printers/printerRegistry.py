
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
    def check_api( cls, printer ):
        """ this class method checks that the printer the user registers
        implements the required methods

        :raises NotImplementedError: if method is missing
        """
        required_methods = [ "setOutPutFile", "addObj", "flush" ]

        missing = [m for m in required_methods if not hasattr(printer, m) or \
            not callable(getattr(printer, m))]
        if missing:
            comment = f"{printer.__name__} is missing methods: {' '.join(missing)}"
            raise NotImplementedError( comment )

    @classmethod
    def register( cls, printer : type[BasicPrinter], printer_label : str,
                  allow_overwrite : bool = False ) -> bool:
        """
        register this printer
        :param printer: Printer object to register
        :param printer_label: Label for the printer, e.g. "python", "xml", "summary", "example"
        :param allow_overwrite: if true, then allow overwriting existing
        entries
        :returns: False if printer already existed, else True
        """
        cls.check_api ( printer )
        if printer_label in cls.printers: # we allow overwrites though
            if allow_overwrite:
                cls.printers[printer_label] = printer
            return False
        cls.printers[printer_label] = printer
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
