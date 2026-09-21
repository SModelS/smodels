from smodels.tools.printers.basicPrinter import BasicPrinter
from smodels.matching.theoryPrediction import TheoryPredictionList

class ExamplePrinter(BasicPrinter):
    """ A simple example of a custom printer, which only prints the analysis ID and the r-value """
    
    def __init__ ( self, output : str = "file", 
                   filename : str = "my.file" ):
        """ constructor, do what you will.
        In this example, we cachine in a self.toPrint
        object cache, write when .flush is called
        """
        BasicPrinter.__init__ ( self, output, filename )
        self.toPrint = []

    def setOutPutFile( self, filename : os.PathLike, overwrite : bool = True,       
        silent : bool = False ): 
        """ need to implement. Can implement your own logic here 
        :param filename: slha filename
        :param overwrite: does the user want to overwrite?
        :param silent: usually used to comment on removing old files
        """
        
        self.filename = filename + '_example.txt'
        if os.path.isfile(self.filename):
            os.remove(self.filename)

    def addObj ( self, obj ):
        """ add an object, either do something immediately with it,
        or write to an object cache 
        """
        self.toPrint.append ( obj )
        
    def _formatTheoryPredictionList(self, obj: object) -> dict:
        """
        Format data of the TheoryPredictionList object.

        :param obj: A TheoryPredictionList object to be printed.
        """
        obj.sortTheoryPredictions()
        
        outputDict = {}
        for theoryPrediction in obj._theoryPredictions:
            expID = theoryPrediction.analysisId()
            r = theoryPrediction.getRValue()
            outputDict[expID] = r
            
        return outputDict
            

    def flush ( self ) -> dict:
        """ this method is called at the end of a model point """
        if hasattr(self,"rmin") and self.rmin > 0:
            rmin = self.rmin
        else:
            rmin = 0.0

        printerOutput = f"Example printer (r > {rmin:1.1e}):\n"
        for obj in self.toPrint:
            if not isinstance(obj,TheoryPredictionList):
                continue
            output = self._formatObj(obj)
            if not output:
                continue
            if not isinstance(output,dict):
                continue
            for expID,r in output.items():
                if r < rmin:
                    continue
                printerOutput += f"{expID} : r-value = {r:1.3f}\n"
        with open(self.filename, "a") as outfile:
            outfile.write(printerOutput)
            outfile.close()
        return {}

# these lines register the printer with smodels, to handle the "example" extension
# The handle name has to be defined in [printer]:outputType in order for the printer to be called.
from smodels.tools.printers.printerRegistry import PrinterRegistry                  
PrinterRegistry.register ( ExamplePrinter, "example" )

