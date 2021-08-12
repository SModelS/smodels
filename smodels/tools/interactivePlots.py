"""
.. module:: interactive_plots
   :synopsis: Main module of the interactive plots.

.. moduleauthor:: Humberto Reyes <humberto.reyes-gonzalez@lpsc.in2p3.fr>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Sabine Kraml <sabine.kraml@gmail.com>

"""

from __future__ import print_function

from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import os,glob,pathlib
import imp
from smodels.tools import interactivePlotsHelpers as helpers
import smodels

class Plotter(object):
    """
    A class to store the required data and produce the interactive plots
    """

    def __init__(self,smodelsFolder,slhaFolder,parameterFile,modelFile=None ):
        """
        Initializes the class.

        :parameter smodelsFolder: path to the folder or tarball containing the smodels 
                                  (python) output files
        :parameter slhaFolder: path to the folder or tarball containing the SLHA input files
        :parameter parameterFile: path to the file containing the plotting definitions
        :parameter modelFile: path to the model file, e.g smodels/share/models/mssm.py
        """
        
        self.data_dict = []
        self.smodelsFolder = smodelsFolder
        self.slhaFolder = slhaFolder
        self.parameterFile = parameterFile
        self.modelFile = modelFile

        self.slha_hover_information = None
        self.ctau_hover_information = None
        self.BR_hover_information = None
        self.SModelS_hover_information = None
        self.plot_data = None
        self.variable_x = None
        self.variable_y = None
        self.plot_list = None
        self.min_BR = None

        if not os.path.isfile(parameterFile):
            raise SModelSError('Parameters file %s not found' %parameterFile)
          
        if modelFile != None:
            if not os.path.isfile(self.modelFile):
                raise SModelSError('model.py file %s not found' % modelFile )
            
        if not os.path.exists(smodelsFolder):
            raise SModelSError("%s not found" %smodelsFolder)
        if not os.path.exists(slhaFolder):
            raise SModelSError("%s not found" %slhaFolder)

        self.loadParameters()
        self.loadModelFile()
        self.initializeDataDict()

    def loadParameters(self):
        """
        Reads the parameters from the plotting parameter file.
        """

        logger.info("Reading parameters from %s ..." %(self.parameterFile))

        parFile = self.parameterFile
        import imp

        try:
            with open(self.parameterFile, 'rb') as fParameters: ## imports parameter file
                parameters = imp.load_module("parameters",fParameters,self.parameterFile,('.py', 'rb', imp.PY_SOURCE))
        # except Exception as e:
        except (IOError,ValueError,ImportError,SyntaxError) as e:
            logger.error("Error loading parameters file %s: %s" % (self.parameterFile,e) )
            raise SModelSError()

        if not hasattr(parameters, 'slha_hover_information'):
            logger.debug("slha_hover_information dictionary was not found in %s. SLHA data will not be included in info box." %parFile)
            self.slha_hover_information = {}
        else:
            self.slha_hover_information = parameters.slha_hover_information

        if not hasattr(parameters, 'ctau_hover_information'):
            logger.debug("ctau_hover_information dictionary was not found in %s. Lifetime data will not be included in info box." %parFile)
            self.ctau_hover_information = {}
        else:
            self.ctau_hover_information = parameters.ctau_hover_information

        if not hasattr(parameters, 'BR_hover_information'):
            logger.debug("BR_hover_information dictionary was not found in %s. Branching ratio data will not be included in info box." %parFile)
            self.BR_hover_information = {}
        else:
            self.BR_hover_information = parameters.BR_hover_information

        if not hasattr(parameters, 'SModelS_hover_information'):
            logger.debug("SModelS_hover_information dictionary was not found in %s. SModelS data will not be included in info box." %parFile)
            self.SModelS_hover_information = {}
        else:
            self.SModelS_hover_information = list(set(parameters.SModelS_hover_information))

        if not hasattr(parameters, 'plot_data'):
            logger.debug("plot_data list was not found in %s. All points will be plotted" %parFile)
            self.plot_data = ['all']
        else:
            self.plot_data = list(set(parameters.plot_data))

        if not hasattr(parameters, 'variable_x'):
            raise SModelSError("variable_x was not found in %s. Please define the variable to be plotted in the x-axis." %parFile)
        else:
            self.variable_x = parameters.variable_x
        if not hasattr(parameters, 'variable_y'):
            raise SModelSError("variable_y was not found in %s. Please define the variable to be plotted in the y-axis." %parFile)
        else:
            self.variable_y = parameters.variable_y
        if not hasattr(parameters, 'plot_list'):
            raise SModelSError("plot_list was not found in %s. Please define the list of plots to be plotted." %parFile)
        else:
            self.plot_list = list(set(parameters.plot_list))

        if not hasattr(parameters,'min_BR'):
            logger.debug("min_BR not found in %s. Will include all decay channels")
            self.min_BR = 'all'
        else:
            self.min_BR = parameters.min_BR

        if not hasattr(parameters,'plot_title'):
            logger.warning("plot_title not defined in %s. Using default title" %parFile)
            self.plot_title = 'interactive-plots'
        else:
            self.plot_title = parameters.plot_title
        

    def loadModelFile(self):
        """
        Reads the parameters from the plotting parameter file.
        """
        
        logger.info("Reading model.py file from %s ..." %(self.modelFile))

        if self.modelFile==None:
            self.particle_names=None
        else:

            try:
                with open(self.modelFile, 'rb') as fparticles:
                    ## imports parameter file
                    self.particle_names = imp.load_module("BSMparticles",fparticles,self.modelFile,('.py', 'rb', imp.PY_SOURCE))
                    # except Exception as e:
            except:
                logger.warning("Error loading model.py file %s , will use pdgs instead.",  self.modelFile)
                self.particle_names=None
           

    def getParticleName(self,pdg):
        """ looks for the particle label in the model.py file """
        found=False
        
        full_list=self.particle_names.BSMList
        for particle in full_list:
            #print(particle.pdg)
            if isinstance(particle,smodels.theory.particle.MultiParticle):
                
                for sub_pdg in particle.pdg:
                    if sub_pdg==pdg:
                        particle_name=particle.label
                        found=True
               
            else:
                if particle.pdg==pdg:
                    particle_name=particle.label
                    found=True
            if found:
                break
        if not found:
            particle_name=pdg
            
        return particle_name

    def editSlhaInformation(self):
        """Edits slha_hover_information,ctau_hover_information,BR_hover_information,variable_x,variable_y if they are defined as a list. The function transforms it in a dict whose keys are the object names """
        
        #variable_x
        if  isinstance(self.variable_x,list):
            variable_x_dict={}
            
            if self.variable_x[0]=='MASS' and self.particle_names!=None:
                particle_name=Plotter.getParticleName(self,self.variable_x[1])
                var_name='m('+particle_name+')'
            else:
                var_name=str(self.variable_x[0])+str(self.variable_x[1])
                    
            variable_x_dict[var_name]=self.variable_x
                
            self.variable_x=variable_x_dict
        
        #variable_y
        if  isinstance(self.variable_y,list):
            variable_y_dict={}
            
            if self.variable_y[0]=='MASS' and self.particle_names!=None:
                particle_name=Plotter.getParticleName(self,self.variable_y[1])
                var_name='m('+particle_name+')'
            else:
                var_name=str(self.variable_y[0])+str(self.variable_y[1])
                    
            variable_y_dict[var_name]=self.variable_y
                
            self.variable_y=variable_y_dict
        
        
        #slha_hover_information
        if  isinstance(self.slha_hover_information,list):
            slha_hover_information_dict={}
            for slha_info in self.slha_hover_information:
                if slha_info[0]=='MASS' and self.particle_names!=None:
                    particle_name=Plotter.getParticleName(self,slha_info[1])
                    var_name='m('+particle_name+')'
                else:
                    var_name=str(slha_info[0])+str(slha_info[1])
                    
                slha_hover_information_dict[var_name]=slha_info
                
            self.slha_hover_information=slha_hover_information_dict
            
        #ctau hover information
        if  isinstance(self.ctau_hover_information,list):
            ctau_hover_information_dict={}
            for slha_info in self.ctau_hover_information:
                if self.particle_names!=None:
                    particle_name=Plotter.getParticleName(self,slha_info)
                    var_name='ctau('+particle_name+')'
                else:
                    var_name='ctau('+str(slha_info)+')'
                    
                ctau_hover_information_dict[var_name]=slha_info
                
            self.ctau_hover_information=ctau_hover_information_dict
        
         #BR hover information
        if  isinstance(self.BR_hover_information,list):
            BR_hover_information_dict={}
            for slha_info in self.BR_hover_information:
                if self.particle_names!=None:
                    particle_name=Plotter.getParticleName(self,slha_info)
                    var_name='BR('+particle_name+')'
                else:
                    var_name='BR('+str(slha_info)+')'
                    
                BR_hover_information_dict[var_name]=slha_info
                
            self.BR_hover_information=BR_hover_information_dict

        return

    def initializeDataDict(self):
        """
        Initializes an empty dictionary with the plotting options.
        """
        Plotter.editSlhaInformation(self)
        self.data_dict = {}
        self.data_dict['SModelS_status']=[]
        for smodels_names in sorted(self.SModelS_hover_information):
            if smodels_names=='SModelS_status':
                continue
            self.data_dict[smodels_names]=[]
        for plot_name in self.plot_list:
            self.data_dict[plot_name]=[]
        for slha_names in self.slha_hover_information:
            self.data_dict[slha_names]=[]
        if list(self.variable_x.keys())[0] not in self.slha_hover_information.keys():
            for variable in self.variable_x:
                self.data_dict[variable]=[]
        if list(self.variable_y.keys())[0] not in self.slha_hover_information.keys():
            for variable in self.variable_y:
                self.data_dict[variable]=[]
        for ctau in self.ctau_hover_information:
            self.data_dict[ctau]=[]
        for BR in self.BR_hover_information:
            self.data_dict[BR]=[]

        self.data_dict['file'] = []
        

    def fillWith(self,smodelsOutput,slhaData):
        """
        Fill the dictionary (data_dict) with the desired data from
        the smodels output dictionary (smodelsDict) and the pyslha.Doc object
        slhaData
        """
        
        filler=helpers.Filler(self,smodelsOutput,slhaData)
        #Fill with smodels data if defined
        if smodelsOutput is None:
            for key in self.SModelS_hover_information:
                if key != 'file':
                    self.data_dict[key].append(False)
        else:
            self.data_dict=filler.getSmodelSData()

        self.data_dict=filler.getSlhaData(self.variable_x,self.variable_y)
        

    def rmFiles ( self, flist ):
        """ remove files in flist """
        for f in flist:
            if os.path.exists ( f ):
                os.remove ( f )

    def loadData(self,npoints=-1):
        """
        Reads the data from the smodels and SLHA folders.
        If npoints > 0, it will limit the number of points in the plot to npoints.

        :parameter npoints: Number of points to be plotted (int).
                            If < 0, all points will be used.
        """
        logger.info( f"Reading data folders {self.smodelsFolder} and {self.slhaFolder} ..." )

        n = 0

        rmfiles = []

        if self.smodelsFolder.endswith(".tar.gz"):
            import tarfile
            with tarfile.open ( self.smodelsFolder, "r:gz" ) as tar:
                tar.extractall()
                files = [ x.name for x in tar.getmembers() ]
                rmfiles += files
                tar.close()
        else:
            files = glob.glob(self.smodelsFolder+'/*')

        slhaFolderIsTarball=False

        if self.slhaFolder.endswith ( ".tar.gz" ):
            slhaFolderIsTarball=True
            import tarfile
            with tarfile.open ( self.slhaFolder, "r:gz" ) as tar:
                tar.extractall()
                slhafiles = [ x.name for x in tar.getmembers() ]
                rmfiles += slhafiles
                tar.close()
        
        for f in files:

            if npoints > 0 and n >= npoints:
                break
            
            smodelsOutput = helpers.importPythonOutput(f)
           
            if not smodelsOutput:
                continue
           
            #Get SLHA file name:
            slhaFile = helpers.getSlhaFile(smodelsOutput)
            files = []
            if not slhaFolderIsTarball:
                slhaFile = os.path.join(self.slhaFolder,os.path.basename(slhaFile))
            slhaData = helpers.getSlhaData(slhaFile)
            if not slhaData:
                continue

            #Data read successfully
            self.data_dict['file'].append(f.split('/')[-1])
            outputStatus = helpers.outputStatus(smodelsOutput)
            
            if outputStatus == -1:
                self.fillWith(None,slhaData)
            else:
                self.fillWith(smodelsOutput,slhaData)
            n += 1

        self.rmFiles ( rmfiles )

        return True

    def display ( self ):
        """ display the pages, works in jupyter notebooks only """
        from IPython.core.display import display,HTML
        mainFile = open( self.outFolder + self.indexfile,'r')
        display(HTML(mainFile.read()))
        mainFile.close()

    def plot(self, outFolder, indexfile = "plots.html" ):
        """
        Uses the data in self.data_dict to produce the plots.

        :parameter outFolder: Path to the output folder.
        :parameter indexfile: name of entry webpage
        """

        if not os.path.isdir(outFolder):
            os.makedirs(outFolder)

        self.outFolder = outFolder
        self.indexfile = indexfile
        logger.info('Making plots...')

        plotter=helpers.PlotlyBackend(self, outFolder )
        plotter.makePlots( indexfile )
        logger.info('Generation of interactive plots finished. Go to: \n %s/%s \n to see the plots.' % ( outFolder, indexfile ) )


def main(args,indexfile= "index.html" ):
    """
    Create the interactive plots using the input from argparse

    :parameter args: argparser.Namespace object containing the options for the plotter

    Main interface for the interactive-plots.

    :parameter smodelsFolder: Path to the folder or tarball containing the 
                              SModelS python output
    :parameter slhaFolder: Path to the folder or tarball containing the SLHA files 
                           corresponding to the SModelS output
    :parameter parameters: Path to the parameter file setting the options for the 
                           interactive plots
    :parameter npoints: Number of points used to produce the plot. If -1, all points 
                        will be used.
    :parameter verbosity: Verbosity of the output (debug,info,warning,error)
    :parameter indexfile: name of the starting web page (index.html)

    :return: True if the plot creation was successfull
    """

    #First check if the needed directories are there
    #inputdirSlha = os.path(args.slhaFolder)
    #args.modelFile='/Users/humberto/Documents/work/smodels-iplots/github/smodels/smodels/share/models/mssm.py'

    if not os.path.exists(args.slhaFolder):
        raise SModelSError("slha directory: "+str(args.slhaFolder)+"' does not exist")

    if not os.path.exists(args.smodelsFolder):
        raise SModelSError("directory of SModelS python output files: "+str(args.smodelsFolder)+"' does not exist")
    if not os.path.exists ( args.outputFolder ):
        os.mkdir ( args.outputFolder )
    if os.path.isdir(args.outputFolder)==False:
        raise SModelSError(f"output directory '{args.outputFolder}' does not exist or is a file")

    if os.path.isfile(args.parameters)==False:
        raise SModelSError("parameter file '"+str(args.parameters)+"' does not exist")
    
    if args.modelFile != None:
        if os.path.isfile(args.modelFile)==False:
            raise SModelSError("model file '"+str(args.modelFile)+"' does not exist")

    #Basic checks:
    smodelsFolder = args.smodelsFolder
    slhaFolder = args.slhaFolder
    parFile = args.parameters
    modelFile=args.modelFile
    verbosity=args.verbosity
    outputFolder=args.outputFolder
    npoints=args.npoints

    try:
        import plotly
    except ImportError:
        raise SModelSError("Plotly is not installed. To use this tool, please install plotly")

    try:
        import pandas
    except ImportError:
        raise SModelSError("Pandas is not installed. To use this tool, please install pandas")

    setLogLevel(verbosity)

    plotter = Plotter(smodelsFolder,slhaFolder,parFile, modelFile )
    loadData = plotter.loadData(npoints)
    if not loadData:
        raise SModelSError("Error loading data from folders:\n %s\n %s" %\
                            (smodelsFolder,slhaFolder))
    plotter.plot(outputFolder, indexfile )
    return outputFolder
