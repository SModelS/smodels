"""
.. module:: crashReport
   :synopsis: Facility used in runSModelS.py to create and read SModelS crash report files.

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

import os
from datetime import datetime
import platform
import traceback
from smodels.installation import installDirectory
from smodels.tools.smodelsLogging import logger

class CrashReport(object):
    """
    Class that handles all crash report information.
    
    """    
    def __init__(self):        
        timestamp = datetime.now().strftime('%Y%m%d%H%M%S%f')
        self.timestampHuman = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
        self.crashReportFileName = 'smodels-' + timestamp + '.crash'


    def createCrashReportFile(self, inputFileName, parameterFileName):
        """
        Create a new SModelS crash report file.
        
        A SModelS crash report file contains:
        
        - a timestamp
        - SModelS version
        - platform information (CPU architecture, operating system, ...)
        - Python version
        - stack trace
        - input file name
        - input file content
        - parameter file name
        - parameter file content
        
        :param inputFileName: relative location of the input file
        :param parameterFileName: relative location of the parameter file
        
        """
    
        with open ( installDirectory()+'/smodels/version', 'r') as versionFile:
            version = versionFile.readline()
    
        with open(inputFileName, 'r') as inputFile:
            inputFileContent = inputFile.read()
    
        with open(parameterFileName, 'r') as parameterFile:
            parameterFileContent = parameterFile.read()
    
        with open(self.crashReportFileName, 'w') as crashReportFile:
            crashReportFile.write("="*80+"\n")
            crashReportFile.write("SModelS Crash Report File\n")
            crashReportFile.write("="*80+"\n")
            crashReportFile.write("Timestamp: " + self.timestampHuman + "\n\n")
            crashReportFile.write("SModelS Version: " + version + "\n")
            crashReportFile.write("Platform: " + platform.platform() + "\n")
            crashReportFile.write("Python Version: " + platform.python_version() + "\n\n")
            crashReportFile.write("="*80+"\n\n")
            crashReportFile.write("-"*80+"\n")
            crashReportFile.write("* Output\n")
            crashReportFile.write("-"*80+"\n\n")
            crashReportFile.write(traceback.format_exc() + "\n\n")
            crashReportFile.write("-"*80+"\n")
            crashReportFile.write("* Input File\n")
            crashReportFile.write("  " + os.path.basename(inputFileName) + "\n")
            crashReportFile.write("-"*80+"\n\n")
            crashReportFile.write(inputFileContent + "\n")
            crashReportFile.write("-"*80+"\n")
            crashReportFile.write("* Parameter File\n")
            crashReportFile.write("  " + os.path.basename(parameterFileName) + "\n")
            crashReportFile.write("-"*80+"\n\n")
            crashReportFile.write(parameterFileContent)
    
    def createUnknownErrorMessage(self):
        """
        Create a message for an unknown error.
        
        """
        message = ("\n\n\n" +"="*80+ "\n\n"
          "SModelS quit unexpectedly due to an unforeseen error.\n"
          "The error has been written to\n"
          + self.crashReportFileName + ".\n\n"
          "If you want to help make SModelS better, then please send this file to\n"
          "smodels-users@lists.oeaw.ac.at and shortly describe what you did!\n\n"
          "Alternatively, use the '--development' option when running runSModelS.py\n"
          "to prevent this message from showing up again.\n\n"
          + 80*"=" )
        return message
    
    
def readCrashReportFile(crashReportFileName):
    """
    Read a crash report file to use its input and parameter file sections for a
    SModelS run.
    
    :param crashReportFileName: relative location of the crash report file
    
    """
    with open(crashReportFileName, 'r') as crashReportFile:
        crashReportFileContent = crashReportFile.readlines()
        
    lineNumber = 0
    inputStartLine = 0
    inputEndLine = 0
    parameterStartLine = 0
    
    for line in crashReportFileContent:
        if lineNumber == 1:
            if not line.rstrip() == "SModelS Crash Report File":
                logger.error("ERROR: Not a SModelS crash report file!")
                break
            
        if line.rstrip() == "* Input File":
            inputStartLine = lineNumber + 4
            
        if line.rstrip() == "* Parameter File":
            inputEndLine = lineNumber - 2
            parameterStartLine = lineNumber + 4
        
        lineNumber += 1
        
    parameterEndLine = lineNumber
    
    crashReportInputFileName = 'crash_report_input'
    crashReportParameterFileName = 'crash_report_parameter'
    
    crashReportInputFile = open(crashReportInputFileName, 'w')
    crashReportParameterFile = open(crashReportParameterFileName, 'w')
    
    for i in range(inputStartLine, inputEndLine):
        crashReportInputFile.write(crashReportFileContent[i])
    
    for i in range(parameterStartLine, parameterEndLine):
        crashReportParameterFile.write(crashReportFileContent[i])
        
    crashReportInputFile.close()
    crashReportParameterFile.close()
    
    return crashReportInputFileName, crashReportParameterFileName
    
    
def createStackTrace():
    """
    Return the stack trace.
    
    """
    return traceback.format_exc()

