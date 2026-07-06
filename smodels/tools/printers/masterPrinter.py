
"""
.. module:: masterPrinter
   :synopsis: Class to handle the distinct printer formats

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import os
from smodels.base.smodelsLogging import logger
from smodels.tools.printers.pythonPrinter import PyPrinter
from smodels.tools.printers.xmlPrinter import XmlPrinter
from smodels.tools.printers.summaryPrinter import SummaryPrinter
from smodels.tools.printers.txtPrinter import TxTPrinter
from smodels.tools.printers.slhaPrinter import SLHAPrinter

from smodels.tools.printers.printerRegistry import PrinterRegistry
PrinterRegistry.register ( PyPrinter, "python" )
PrinterRegistry.register ( XmlPrinter, "xml" )
PrinterRegistry.register ( SummaryPrinter, "summary" )
PrinterRegistry.register ( TxTPrinter, "log" )
PrinterRegistry.register ( TxTPrinter, "stdout" )
PrinterRegistry.register ( SLHAPrinter, "slha" )

class MPrinter(object):
    """
    Master Printer class to handle the Printers (one printer/output type)
    """

    def __init__(self):

        self.name = "master"
        self.Printers = {}
        self.outputFormat = "version3"

    def setPrinterOptions(self, parser):
        """
        Define the printer types and their options.

        :param parser: ConfigParser storing information from the parameters file
        """
        # Define the printer types and the printer-specific options:
        printerTypes = [prt.strip() for prt in parser.get(
            "printer", "outputType").split(",")]

        # Copy stdout options to log options:
        if 'log' in printerTypes:
            if parser.has_section('stdout-printer') and not \
                    parser.has_section('log-printer'):
                parser.add_section('log-printer')
                for option, val in parser.items('stdout-printer'):
                    parser.set('log-printer', option, val)

        if parser.has_option("printer","outputFormat"):
            self.outputFormat = parser.get("printer","outputFormat")
        for prt in printerTypes:
            if not PrinterRegistry.has ( prt ):
                logger.warning(f"Unknown printer format: {str(prt)}")
                continue
            PrinterClass = PrinterRegistry.get(prt)

            if prt == 'python':
                newPrinter = PrinterClass(output='file', 
                                          outputFormat=self.outputFormat)
            elif prt == 'summary':
                newPrinter = PrinterClass(output='file', 
                                          outputFormat=self.outputFormat)
            elif prt == 'stdout':
                newPrinter = PrinterClass(output='stdout', 
                                          outputFormat=self.outputFormat)
            elif prt == 'log':
                newPrinter = PrinterClass(output='file', 
                                          outputFormat=self.outputFormat)
            elif prt == 'xml':
                newPrinter = PrinterClass(output='file', 
                                          outputFormat=self.outputFormat)
            elif prt == 'slha':
                newPrinter = PrinterClass(output='file', 
                                          outputFormat=self.outputFormat)
                if parser.getboolean("options", "doCompress") or \
                        parser.getboolean("options", "doInvisible"):
                    newPrinter.docompress = 1
                if parser.has_option("options", "combineSRs") and \
                        parser.getboolean("options", "combineSRs"):
                    newPrinter.combinesr = 1
                if parser.has_option("options", "combineAnas") and \
                        parser.get("options", "combineAnas"):
                    newPrinter.combineanas = 1
            else:
                newPrinter = PrinterClass()

            # Set printer-specific options:
            if parser.has_section(prt+'-printer'):
                newPrinter.setOptions(parser.items(prt+'-printer'))
            self.Printers[prt] = newPrinter

    def addObj(self, obj):
        """
        Adds the object to all its Printers:

        :param obj: An object which can be handled by the Printers.
        """
        for prt in self.Printers.values():
            prt.addObj(obj)

    def setOutPutFiles(self, filename : os.PathLike, silent : bool = False ):
        """
        Set the basename for the output files. Each printer will
        use this file name appended of the respective extension
        (i.e. .py for a python printer, .smodels for a summary printer,...)

        :param filename: Input file name
        :param silent: dont comment removing old files
        """

        for printer in self.Printers.values():
            printer.setOutPutFile(filename, silent=silent)

    def flush(self) -> dict:
        """
        Ask all printers to write the output and clear their cache.
        If the printers return anything other than None,
        we pass it on.
        """
        ret = {}
        for printerType, printer in self.Printers.items():
            ret[printerType] = printer.flush()
        return ret
