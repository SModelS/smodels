#!/usr/bin/env python3

"""
.. module:: interactive_plots
   :synopsis: Main module of the interactive plots.
   
   .. moduleauthor:: Humberto Reyes <humberto.reyes-gonzalez@lpsc.in2p3.fr>
   .. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
   .. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
"""

from __future__ import print_function


from smodels.tools import interactivePlotsHelpers 
from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import os, imp, glob

class DataHolder(object):
    """
    A simple object to store the required data for producing the interactive plots
    """
    
    def __init__(self,smodelsFolder,slhaFolder,parameterFile):
        self.data_dict = []
        self.smodelsFolder = smodelsFolder
        self.slhaFolder = slhaFolder
        self.parameterFile = parameterFile
        
        self.slha_hover_information = None 
        self.ctau_hover_information = None
        self.BR_hover_information = None
        self.smodels_hover_information = None
        self.plot_data = None
        self.variable_x = None 
        self.variable_y = None
        self.plot_list = None
        self.BR_get_top = None
        
        
        if not os.path.isfile(parameterFile):
            raise SModelSError('Parameters file %s not found' %parameterFile)
        if not os.path.isdir(smodelsFolder):
            raise SModelSError("Folder %s not found" %smodelsFolder)
        if not os.path.isdir(slhaFolder):
            raise SModelSError("Folder %s not found" %slhaFolder)

        self.loadParameters()
        self.initializeDataDict()

    def loadParameters(self):
        
        logger.info("Reading parameters from %s ..." %(self.parameterFile))        
        
        parFile = self.parameterFile
        
        try:
            with open(self.parameterFile, 'rb') as fParameters: ## imports parameter file
                parameters = imp.load_module("parameters",fParameters,self.parameterFile,('.py', 'rb', imp.PY_SOURCE))
        except:
            logger.error("Error loading parameters file %s" %self.parameterFile)
            return False
         
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
    
        if not hasattr(parameters, 'smodels_hover_information'):
            logger.debug("smodels_hover_information dictionary was not found in %s. SModelS data will not be included in info box." %parFile)
            self.smodels_hover_information = {}
        else:
            self.smodels_hover_information = list(set(parameters.smodels_hover_information))
    
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
            
        if not hasattr(parameters,'BR_get_top'):
            logger.debug("BR_get_top not found in %s. Will include all decay channels")
            self.BR_get_top = 'all'
        else:
            self.BR_get_top = parameters.BR_get_top

        if not hasattr(parameters,'plot_title'):
            logger.warning("plot_ttiel not defined in %s. Using default title" %parFile)
            self.plot_title = 'interactive-plots'
        else:
            self.plot_title = parameters.plot_title


    def initializeDataDict(self):
        self.data_dict = {}
        for smodels_names in self.smodels_hover_information:
            self.data_dict[smodels_names]=[]
        for plot_name in self.plot_list:
            self.data_dict[plot_name]=[]
        for slha_names in self.slha_hover_information:
            self.data_dict[slha_names]=[] 
        for variable in self.variable_x:
            self.data_dict[variable]=[]
        for variable in self.variable_y:
            self.data_dict[variable]=[] 
        for ctau in self.ctau_hover_information:
            self.data_dict[ctau]=[] 
        for BR in self.BR_hover_information:
            self.data_dict[BR]=[]
            
        self.data_dict['file'] = []
        
    def fillWith(self,smodelsDict,slhaData):
        """
        Fill the dictionary (data_dict) with the desired data from
        the smodels output dictionary (smodelsDict) and the pyslha.Doc object
        slhaData
        """
    
        #Fill with smodels data if defined
        if smodelsDict is None:
            for key in self.data_dict:
                if key != 'file':
                    self.data_dict[key].append(None)                
        else:    
            self.data_dict = interactivePlotsHelpers.get_expres(self.data_dict,smodelsDict)
            self.data_dict = interactivePlotsHelpers.get_missed_topologies(self.data_dict,smodelsDict)
            self.data_dict = interactivePlotsHelpers.get_asymmetric_branches(self.data_dict,smodelsDict)
            self.data_dict = interactivePlotsHelpers.get_outside_grid(self.data_dict,smodelsDict)


            #Fill with SLHA data:
            self.data_dict =  interactivePlotsHelpers.get_slha_hover_info(self.data_dict,slhaData,
                                                                          self.slha_hover_information)
            self.data_dict = interactivePlotsHelpers.get_ctau(self.data_dict,slhaData,
                                                              self.ctau_hover_information)             
            self.data_dict = interactivePlotsHelpers.get_BR(self.data_dict,slhaData,self.BR_hover_information,
                                                            self.BR_get_top)     
            #Fill with the x and y data:
            self.data_dict = interactivePlotsHelpers.get_variable(self.data_dict,slhaData,self.BR_hover_information,
                                                            self.variable_x) 
            self.data_dict = interactivePlotsHelpers.get_variable(self.data_dict,slhaData,self.BR_hover_information,
                                                            self.variable_y) 
     
    def loadData(self,npoints=-1):
        """
        Reads the data from the smodels and SLHA folders.
        If npoints > 0, it will limit the number of points in the plot to npoints.
        
        :parameter npoints: Number of points to be plotted (int). If < 0, all points will be used.
        """
        
        logger.info("Reading data folders %s and %s ..." %(self.smodelsFolder,self.slhaFolder))
        
        n = 0
        for f in glob.glob(self.smodelsFolder+'/*'):
            
            if npoints > 0 and n >= npoints:
                break
            
            smodelsOutput = interactivePlotsHelpers.import_python_output(f)
            if not smodelsOutput:
                continue
            
            #Get SLHA file name:
            slhaFile = interactivePlotsHelpers.get_slha_file(smodelsOutput)
            slhaFile = os.path.join(self.slhaFolder,os.path.basename(slhaFile))
            #Get SLHA data:
            slhaData = interactivePlotsHelpers.get_slha_data(slhaFile)
            if not slhaData:
                continue
            
            #Data read successfully
            self.data_dict['file'].append(f)       
            outputStatus = interactivePlotsHelpers.output_status(smodelsOutput)
            if outputStatus == -1:
                self.fillWith(None,slhaData)
            else:
                self.fillWith(smodelsOutput,slhaData)
            n += 1

        return True
    
    def makePlots(self,outFolder):
        """
        Uses the data in self.data_dict to produce the plots.
        
        :parameter outFolder: Path to the output folder.
        """
        
        
        if not os.path.isdir(outFolder):
            os.makedirs(outFolder)
        
        
        logger.info('Making plots...')
          
        data_frame_all = interactivePlotsHelpers.make_data_frame(self.data_dict)
     
        data_frame_all = interactivePlotsHelpers.fill_hover(data_frame_all,
                                                            self.smodels_hover_information,
                                                            self.slha_hover_information,
                                                            self.ctau_hover_information,
                                                            self.BR_hover_information) 
     
        data_frame_excluded,data_frame_nonexcluded = interactivePlotsHelpers.data_frame_excluded_nonexcluded(data_frame_all) 
        x_axis,y_axis = interactivePlotsHelpers.get_xy_axis(self.variable_x,self.variable_y) 
        cont_plots,disc_plots = interactivePlotsHelpers.separate_cont_disc_plots(self.plot_list,self.data_dict) 
     
        interactivePlotsHelpers.make_continuous_plots_all(cont_plots,x_axis,
                                                          y_axis,outFolder,data_frame_all,self.plot_data,
                                                          self.plot_title)
         
        interactivePlotsHelpers.make_continuous_plots_excluded(cont_plots,x_axis,
                                                               y_axis,outFolder,data_frame_excluded,self.plot_data,
                                                               self.plot_title)
         
        interactivePlotsHelpers.make_continuous_plots_nonexcluded(cont_plots,x_axis,
                                                                  y_axis,outFolder,data_frame_nonexcluded,self.plot_data,
                                                                  self.plot_title)
         
        interactivePlotsHelpers.make_discrete_plots_all(disc_plots,x_axis,y_axis,
                                                        outFolder,data_frame_all,self.plot_data,
                                                        self.plot_title)
         
        interactivePlotsHelpers.make_discrete_plots_excluded(disc_plots,x_axis,y_axis,
                                                             outFolder,data_frame_excluded,self.plot_data,
                                                             self.plot_title)
         
        interactivePlotsHelpers.make_discrete_plots_nonexcluded(disc_plots,x_axis,y_axis,
                                                                outFolder,data_frame_nonexcluded,self.plot_data,
                                                                self.plot_title)
        
        interactivePlotsHelpers.create_main_html(outFolder,self.plot_data,self.plot_list)
        
        logger.info('Generation of interactive plots finished. Go to: \n %s/main.html \n to see the plots.' %outFolder)
        
        return True




def main(args):
    """
    Main interface for the interactive-plots. 
    """

    try:
        import plotly
    except ImportError as e:
        raise SModelSError("Plotly is not installed. To use this tool, please install plotly")

    try:
        import pandas
    except ImportError as e:
        raise SModelSError("Pandas is not installed. To use this tool, please install pandas")

    setLogLevel(args.verbosity)

    #Basic checks:
    smodelsFolder = args.smodelsFolder
    slhaFolder = args.slhaFolder
    parFile = args.parameters
    
    dataHolder = DataHolder(smodelsFolder,slhaFolder,parFile)
    loadData = dataHolder.loadData(args.npoints)
    if not loadData:
        raise SModelSError("Error loading data from folders:\n %s\n %s" %(smodelsFolder,slhaFolder))
    
    dataHolder.makePlots(args.outputFolder)
    
    return True
