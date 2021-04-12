"""
.. module:: interactivePlotsHelpers
   :synopsis: Main functions for producing interactive plots
   
.. moduleauthor:: Humberto Reyes <humberto.reyes-gonzalez@lpsc.in2p3.fr>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Sabine Kraml <sabine.kraml@gmail.com>
   
"""
import sys
import math
import os
from smodels.tools.smodelsLogging import logger
try:
    import warnings ## we can ignore these warnings, see
    # https://stackoverflow.com/questions/40845304/runtimewarning-numpy-dtype-size-changed-may-indicate-binary-incompatibility
    warnings.filterwarnings("ignore", message="numpy.dtype size changed")
    import plotly.graph_objs as go
    import plotly,imp
except ImportError as e:
    logger.error ( "could not import plotly. Please try to install it, e.g. via:\npip install --user plotly" )
    sys.exit()
try:
    import pandas as pd
except ImportError as e:
    logger.error ( "could not import pandas. Please try to install it, e.g. via:\npip install --user pandas" )
    sys.exit()
import os
import pyslha
from smodels.theory.exceptions import SModelSTheoryError as SModelSError

def importPythonOutput(smodelsFile):
    """
    Imports the smodels output from each .py file.
    """
       
    try:
        with open(smodelsFile, 'rb') as fsmodels: ## imports smodels file
            smodelsOut = imp.load_module("smodelsOutput",fsmodels,smodelsFile,('.py', 'rb', imp.PY_SOURCE))
            smodelsOutput = smodelsOut.smodelsOutput
    except (ImportError,AttributeError,IOError,ValueError,OSError,SyntaxError):
        logger.debug("Error loading smodels file %s. Does it contain a smodelsOutput dictionary?" %smodelsFile)
        return False
    
    if not isinstance(smodelsOutput,dict):    
        logger.warning("smodelsOutput in file %s is not a dictionary." %smodelsFile)
        return False
    
    return smodelsOutput


def getEntry(inputDict,*keys):
    """
    Get entry key in dictionary inputDict.
    If a list of keys is provided, it will assumed nested
    dictionaries (e.g. key1,key2 will return inputDict[key1][key2]). 
    """
    
    keys = list(keys)
    
    if not keys:
        return inputDict
    
    key = keys.pop(0)
    if not key in inputDict:
        logger.debug('Key %s not found in input dictionary' %key)
        return False
    else:
        return getEntry(inputDict[key],*keys)
 
def outputStatus(smodelsDict):
    """
    Check the smodels output status in the file, if it's -1, 
    it will append 'none' to each list in the dictionary.
    """
    
    outputStatus = getEntry(smodelsDict, 'OutputStatus', 'file status')
    if outputStatus is False:
        raise SModelSError()

    return outputStatus

def getSlhaFile(smodelsDict):
    """
    Returns the file name of the SLHA file corresponding to the output in smodelsDict
    """

    slhaFile = getEntry(smodelsDict,'OutputStatus','input file')
    if not slhaFile:
        raise SModelSError()

    return os.path.basename(slhaFile)

def getSlhaData(slhaFile):
    """
    Uses pyslha to read the SLHA file. Return a pyslha.Doc objec, if successful.
    """

    if not os.path.isfile(slhaFile):
        logger.warning("SLHA file %s not found. This point will be ignored" %slhaFile)
        return False
    
    try:
        slhaData = pyslha.readSLHAFile(slhaFile)
    except (OSError,pyslha.ParseError) as e:
        logger.warning("Error reading SLHA file %s: %s" % (slhaFile,e) )
        return False
    
    return slhaData
    
def truncate(number):
    "truncate float to 3 decimal places"
    
    if number<1 and number>1e-90:
        
        provisional_number=number
        order=0
       
        while provisional_number<1:
           
            provisional_number=provisional_number*(10)
            order=order+1
            
        
        factor=10**3
        truncated=(math.trunc(provisional_number * factor) / (factor*10**(order)))
   
    #elif number>1e-90:
    # truncated=0.0
    
    else:
        factor=10**4
        truncated=math.trunc(number * factor) / factor
    return truncated
     
def getExpres(data_dict,smodelsOutput):
    """
    Extracts the Expres info from the .py output. If requested, the data will be appended on each corresponding list
    """   
    
    rmax=0
    decompStatus = getEntry(smodelsOutput,'OutputStatus','decomposition status')
    
    if decompStatus == 1:
        expResList = getEntry(smodelsOutput,'ExptRes')
        if not expResList or not isinstance(expResList,list):
            print(smodelsOutput)
            raise SModelSError("Error reading ExptRes.")    
        for expres in expResList:
            if 'r' in expres:
                
                r = getEntry(expres,'r')
        
            else:
                r = getEntry(expres,'theory prediction (fb)')/getEntry(expres,'upper limit (fb)')
                r = truncate(r)
            if r>rmax:
                rmax = r
                Txname = getEntry(expres,'TxNames')
                Txname=','.join(Txname)
                if getEntry(expres,'chi2')==False:
                    chi_2=False
                else:    
                    chi_2 = getEntry(expres,'chi2')
                    chi_2 = truncate(chi_2)

                analysis_id = getEntry(expres,'AnalysisID')
    else:
        Txname = False
        chi_2 = False
        analysis_id = False
        rmax=False
 
    if 'SModelS_status' in  data_dict:
        if rmax ==False:
            data_dict['SModelS_status'].append(False)
        elif rmax>1:
            data_dict['SModelS_status'].append('Excluded')  
        else:
            data_dict['SModelS_status'].append('Non-excluded')  

    if 'r_max' in data_dict.keys():
        data_dict['r_max'].append(rmax)
        
    if 'Tx' in data_dict.keys(): 
        data_dict['Tx'].append(Txname)
        
    if 'chi2' in data_dict.keys(): 
        data_dict['chi2'].append(chi_2) 
    if 'Analysis' in data_dict.keys(): 
        data_dict['Analysis'].append(analysis_id)   
    return data_dict;   
  
def getMissedTopologies(data_dict,smodelsOuptut):
    """
    Extracts the Missed topologies info from the .py output. If requested, the data will be appended on each corresponding list
    """
    decompStatus = getEntry(smodelsOuptut,'OutputStatus','decomposition status')
    missedtopo_max_xsec=0
    missedtopo_total_xsec=0
    mt_max = 'False'
    if decompStatus >= 0:
        for missed_topo in getEntry(smodelsOuptut, 'Missed Topologies'):
            missedtopo_xsec=missed_topo.get('weight (fb)')            
            missedtopo_total_xsec=missedtopo_total_xsec+missedtopo_xsec
            if missedtopo_xsec>missedtopo_max_xsec:
                missedtopo_max_xsec=missedtopo_xsec
                mt_max=missed_topo.get('element')
        missedtopo_max_xsec=truncate(missedtopo_max_xsec)
        missedtopo_total_xsec=truncate(missedtopo_total_xsec)
        
        
    else:
        missedtopo_max_xsec=False
        missedtopo_total_xsec=False
        missed_topo = False
    if 'MT_max' in data_dict.keys():
        data_dict.get('MT_max').append(str(mt_max))
    if 'MT_max_xsec' in data_dict.keys():
        data_dict.get('MT_max_xsec').append(missedtopo_max_xsec)    
    if 'MT_total_xsec' in data_dict.keys():
        data_dict.get('MT_total_xsec').append(missedtopo_total_xsec)     
        
    return data_dict           
   
def getLongCascades(data_dict,smodelsOutput):
    """
    Extracts the Long cascade info from the .py output. If requested, the data will be appended on each corresponding list
    """
    decompStatus = getEntry(smodelsOutput,'OutputStatus','decomposition status')
    long_cascade_total_xsec=0
    if decompStatus >= 0:
        long_cascade_total_xsec=0
        for long_cascade in smodelsOutput.get('Long Cascades'):
            long_cascade_xsec=long_cascade.get('weight (fb)')
            long_cascade_total_xsec=long_cascade_total_xsec+long_cascade_xsec
        long_cascade_total_xsec=truncate(long_cascade_total_xsec)
    
    else:
        long_cascade_total_xsec=False
    if 'MT_long_xsec' in data_dict.keys():
        data_dict.get('MT_long_xsec').append(long_cascade_total_xsec)
             
    return data_dict
             
def getAsymmetricBranches(data_dict,smodelsOutput):
    """
    Extracts the asymmetric branches info from the .py output. If requested, the data will be appended on each corresponding list
    """
    decompStatus = getEntry(smodelsOutput,'OutputStatus','decomposition status')
    if decompStatus >= 0:
        asymmetric_branch_total_xsec = sum([asym_br['weight (fb)'] for asym_br in smodelsOutput['Asymmetric Branches']])
        asymmetric_branch_total_xsec=truncate(asymmetric_branch_total_xsec)
    else:
        asymmetric_branch_total_xsec=False             
    if 'MT_asym_xsec' in data_dict.keys():
        data_dict.get('MT_asym_xsec').append(asymmetric_branch_total_xsec) 
    return data_dict  
    
def getOutsideGrid(data_dict,smodelsOutput):
    """
    Extracts the outside grid info from the .py output. If requested, the data will be appended on each corresponding list.
    """   
    decompStatus = getEntry(smodelsOutput,'OutputStatus','decomposition status')
    outside_grid_total_xsec=0
    if decompStatus >= 0:
        for outside_grid in smodelsOutput.get('Outside Grid'):
            outside_grid_xsec=outside_grid.get('weight (fb)') 
            outside_grid_total_xsec = outside_grid_total_xsec+outside_grid_xsec
        outside_grid_total_xsec=truncate(outside_grid_total_xsec)
    else:
        outside_grid_total_xsec = False      
    if 'MT_outgrid_xsec' in data_dict.keys():
        data_dict.get('MT_outgrid_xsec').append(outside_grid_total_xsec)  
    return data_dict   

def getSlhaHoverInfo(data_dict,slhaData,slha_hover_information):
        """
        Gets the requested slha info from eachh slha file, to fill the hover.
        """
          
        for key in self.slha_hover_information.keys():
            block = self.slha_hover_information.get(key)[0]
            code_number = self.slha_hover_information.get(key)[1]
            if block=='MASS':
                data_dict.get(key).append(truncate(abs(slhaData.blocks[block][code_number])))
            else:
                data_dict.get(key).append(truncate(slhaData.blocks[block][code_number]))
            
        return data_dict
        
def getCtau(data_dict,slhaData,ctau_hover_information):
    """
    Computes the requested ctaus, that will go into de hover.
    """   
    
    for key in ctau_hover_information.keys():
        value=ctau_hover_information.get(key)
        total_width=float(str(slhaData.decays[value]).split('total width = ')[1].split(' GeV')[0])
        if total_width == 0:
            ctau=float('inf')
        else:
            mean_lifetime=(6.582119e-16)/(total_width*1e9)
            ctau=(mean_lifetime)*(299792458)
            
            ctau=truncate(ctau)
            
        data_dict.get(key).append(ctau)
        
    return data_dict 

              
              
def getBR(data_dict,slhaData,BR_hover_information,BR_get_top):
    """
    Gets the requested branching ratios from the slha file, that will go into de hover.
    """   
    for key in BR_hover_information.keys():
        pdg_number=BR_hover_information.get(key) 
        BR=str(slhaData.decays[pdg_number]).split('\n')[1:]
        if BR_get_top != 'all':
            BR=BR[:int(BR_get_top)]  
        BR=str(BR).split(' ')      
        BR=''.join(BR)
        data_dict.get(key).append(BR)
        
    return data_dict 

def getVariable(data_dict,slhaData,slha_hover_information,variable):
    """
    Gets the variable from the slha file.
    """  
    for key in variable.keys():
        if str(key) not in slha_hover_information.keys():
            block=variable.get(key)[0]
            code_number=variable.get(key)[1]
            if block=='MASS':
                data_dict.get(key).append(truncate(abs(slhaData.blocks[block][code_number])))
            else:
                data_dict.get(key).append(truncate(slhaData.blocks[block][code_number]))
            
    return data_dict


class Plotter(data_dict,SModelS_hover_information,slha_hover_information,ctau_hover_information,BR_hover_information,variable_x,variable_y,plot_list,plot_data,plot_title):
               
    def __init__ ( self, data_dict,SModelS_hover_information,
               slha_hover_information,ctau_hover_information,BR_hover_information,variable_x,variable_y,plot_list,plot_data,plot_title):
               
               
        self.data_dict=data_dict
        self.SModelS_hover_information=SModelS_hover_information
        self.slha_hover_information=slha_hover_information
        self.ctau_hover_information=ctau_hover_information
        self.BR_hover_information=BR_hover_information
        self.variable_x=variable_x
        self.variable_y=variable_y
        self.plot_list=plot_list
        self.plot_data=plot_data
        self.plot_title=plot_title
        
        
        return

    def makeDataFrame(self):
        """
        Transform the main dictionary in a data frame.
        """
        self.data_frame_all = pd.DataFrame(data=self.data_dict)
        self.data_frame_all.to_csv('all_data_frame.txt', sep=' ', index=False,header=True)
      
        return self.data_frame_all
  

    def refiningVariableNames(self):
        ''' Redifining the output variable names to html format  '''
    
        self.html_names={'SModelS_status':'Smodels status',
                      'r_max':'r<sub>max</sub>',
                      'chi2':' &#967;<sup>2</sup>',
                      'Tx':'T<sub>max</sub>',
                      'Analysis':'Analysis',
                      'MT_max':'MT<sub>max</sub>',
                      'MT_max_xsec':'MT<sub>max xsection</sub>.',
                      'MT_total_xsec':'MT<sub>total xsection</sub>',
                      'MT_long_xsec':'MT<sub>long cascade xsection</sub>',
                      'MT_asym_xsec':'MT<sub>asymmetric branch xsection</sub>',
                      'MT_outgrid_xsec':'MT<sub>outside grid xsection</sub>'}
        
        return self.html_names

    def fillHover(self):
        """ Generates the text of the hover, according to users's requests. """
        self.data_frame_all['hover_text']=''
        for column in self.slha_hover_information:
            if  self.slha_hover_information.get(column)[0]=='MASS':
                self.data_frame_all['hover_text'] = self.data_frame_all['hover_text']+column+': '+self.data_frame_all[column].astype('str')+' GeV'+'<br>'
            else:
                self.data_frame_all['hover_text'] = self.data_frame_all['hover_text']+column+': '+self.data_frame_all[column].astype('str')+'<br>'
        for column in self.ctau_hover_information:
            self.data_frame_all['hover_text'] = self.data_frame_all['hover_text']+column+': '+self.data_frame_all[column].astype('str')+' m'+'<br>'
        data_frame_br=pd.DataFrame()
        for column in self.BR_hover_information:
            i=0
            data_frame_br[column] = self.data_frame_all[column]
            for i in range(len(data_frame_all.index)):
                if self.data_frame_all[column][i]:
                    data_frame_br[column][i] = self.data_frame_all[column][i].split("','")
                else:
                    data_frame_br[column][i] = []
                j=0
                brs=''
                for br in data_frame_br[column][i]:
                    if j<=4:
                        brs=brs+' '+br
                    if j>4:
                        brs=brs+'<br>'
                        j=0
                    j=j+1
                data_frame_br[column][i]=brs
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+column+': '+data_frame_br[column].astype('str')+'<br>'
        if 'SModelS_status' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('SModelS_status')+': '+self.data_frame_all['SModelS_status'].astype('str')+'<br>'
        if 'r_max' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('r_max')+': '+self.data_frame_all['r_max'].astype('str')+'<br>'
        if 'Tx' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('Tx')+': '+self.data_frame_all['Tx'].astype('str')+'<br>'
        if 'Analysis' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('Analysis')+': '+self.data_frame_all['Analysis'].astype('str')+'<br>'
        if 'chi2' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('chi2')+': '+self.data_frame_all['chi2'].astype('str')+'<br>'
        if 'MT_max' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_max')+': '+self.data_frame_all['MT_max'].astype('str')+'<br>'
        if 'MT_max_xsec' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_max_xsec')+': '+self.data_frame_all['MT_max_xsec'].astype('str')+' fb'+'<br>'
        if 'MT_total_xsec' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_total_xsec')+': '+self.data_frame_all['MT_total_xsec'].astype('str')+' fb'+'<br>'
        if 'MT_long_xsec' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_long_xsec')+': '+self.data_frame_all['MT_long_xsec'].astype('str')+' fb'+'<br>'
        if 'MT_asym_xsec' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_asym_xsec')+': '+self.data_frame_all['MT_asym_xsec'].astype('str')+' fb'+'<br>'
        if 'MT_outgrid_xsec' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_outgrid_xsec')+': '+self.data_frame_all['MT_outgrid_xsec'].astype('str')+' fb'+'<br>'
        if 'file' in self.SModelS_hover_information:
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+'file'+': '+self.data_frame_all['file'].astype('str')+'<br>'
        

        return self.data_frame_all;
 
 
  
 
    def DataFrameExcludedNonexcluded(self):
        """ Generate sub data frames for excluded and non-excluded points """
        self.data_frame_excluded=self.data_frame_all.loc[self.data_frame_all['SModelS_status']=='Excluded']
        self.data_frame_nonexcluded=self.data_frame_all.loc[self.data_frame_all['SModelS_status']=='Non-excluded']
        return self.data_frame_excluded, self.data_frame_nonexcluded;
 
    def GetXyAxis(self):
        """ Retrieves the names of the x and y axis variables. """
        self.x_axis=list(self.variable_x.keys())[0]
        self.y_axis=list(self.variable_y.keys())[0]
        return self.x_axis,self.y_axis;
     
 
    def SeparateContDiscPlots(self):
        ''' Generate sub lists of plots with discrete and conitnuous z axis variables. '''
        self.cont_plots=[]
        self.disc_plots=[]
        discrete_options=['Tx','Analysis','MT_max','file','SModelS_status']
        for plot in self.plot_list:
            if plot in discrete_options:
                self.disc_plots.append(plot)
            else:
            
                
                self.cont_plots.append(plot)
        return self.cont_plots, self.disc_plots;
    

                
    
 
    def plotDescription(self):
        ''' Generate a description for each plot.'''
        self.plot_descriptions={'SModelS_status':'Excluded or not excluded by SModelS.',
                      'r_max':'highest r-value from SModelS.',
                      'chi2':'&#967;<sup>2</sup> value associated to the highest r-value.',
                      'Tx':'Topology/ies which give the highest r-value.',
                      'Analysis':'Experimental analysis from which the highest r-value comes from.',
                      'MT_max':'Missing topologies with the largest cross section.',
                      'MT_max_xsec':'Cross section of MT_max.',
                      'MT_total_xsec':'Total missing cross section.',
                      'MT_long_xsec':'Missing cross section in long cascade decays.',
                      'MT_asym_xsec':'Missing cross section in decays with asymmetric branches.',
                      'MT_outgrid_xsec':'Missing cross section outside the mass grids of the experimental results.'}
        return self.plot_descriptions;
 
#####continuous plots##############
    def makeContinuousPlotsAll(self):
        """ Generate plots with continuous z axis variables, using all data points """
        if 'all' in self.plot_data:
            for cont_plot in self.cont_plots:
        
                if cont_plot=='chi2':
                    all_false=True
                    for chi2_value in self.data_frame_all['chi2']:
                        if chi2_value!=False:
                            all_false=False
                    if all_false==True:
                        logger.info('No values where found for chi^2. Skipping this plot')
                        continue
                
                plot_desc=self.plot_descriptions.get(cont_plot)
                cont_plot_legend=self.html_names.get(cont_plot)
                if cont_plot=='MT_max_xsec' or cont_plot=='MT_total_xsec' or cont_plot=='MT_long_xsec' or cont_plot=='MT_asym_xsec' or cont_plot=='MT_outgrid_xsec':
                    cont_plot_legend=self.html_names.get(cont_plot)+' (fb)'
                z=self.data_frame_all[cont_plot]
                x=self.data_frame_all[self.x_axis]
                y=self.data_frame_all[self.y_axis]
                hover_text=self.data_frame_all['hover_text']
            
                data = [
                    go.Scatter(
                    x=x,
                    y=y,
                    text=hover_text,
                    hoverinfo='text',
                    mode='markers',
                    marker=dict(
                        size=10,
                        cmax=self.data_frame_all[cont_plot].max(),
                        cmin=self.data_frame_all[cont_plot].max(),
                        color=z,
                        colorbar=dict(
                            title=cont_plot_legend), 
                    colorscale='Jet')
  
                            )
                    ]
               
                if self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis+' (GeV)'),
                                   yaxis = dict(title=self.y_axis+' (GeV)'),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )
            
                elif self.variable_x.get(self.x_axis)[0]!='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis),
                                   yaxis = dict(title=self.y_axis+' (GeV)'),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )            
            
                elif self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]!='MASS':
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis+' (GeV)'),
                                   yaxis = dict(title=self.y_axis),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )             
            
                else:
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis),
                                   yaxis = dict(title=self.y_axis),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )
                
                fig = go.Figure(data=data, layout=layout)
                plotly.offline.plot(fig, filename = path_to_plots+'/'+cont_plot+'_all.html', auto_open=False)
        return;
 
 
    def makeContinuousPlotsExcluded(self):
        """ Generate plots with continuous z axis variables, using excluded data points """
        if 'excluded' in self.plot_data:
            for cont_plot in self.cont_plots:
        
                if cont_plot=='chi2':
                    all_false=True
                    for chi2_value in self.data_frame_excluded['chi2']:
                        if chi2_value!=False:
                            all_false=False
                    if all_false==True:
                        logger.info('No values where found for chi^2 in the excluded region. Skipping this plot')
                        continue
                plot_desc=self.plot_descriptions.get(cont_plot)
                cont_plot_legend=self.html_names.get(cont_plot)
                if cont_plot=='MT_max_xsec' or cont_plot=='MT_total_xsec' or cont_plot=='MT_long_xsec' or cont_plot=='MT_asym_xsec' or cont_plot=='MT_outgrid_xsec':
                    cont_plot_legend=self.html_names.get(cont_plot)+' (fb)'
                z=self.data_frame_excluded[cont_plot]
                x=self.data_frame_excluded[self.x_axis]
                y=self.data_frame_excluded[self.y_axis]
                hover_text=self.data_frame_excluded['hover_text']
             
             
                data = [
                    go.Scatter(
                    x=x,
                    y=y,
                    text=hover_text,
                    hoverinfo='text',
                    marker=dict(
                        size=10,
                        cmax=self.data_frame_excluded[cont_plot].max(),
                        cmin=self.data_frame_excluded[cont_plot].max(),
                        color=z,
                        colorbar=dict(
                            title=cont_plot_legend), 
                    colorscale='Jet'), 
                    mode='markers'  
            )
     
                    ]
     
                if self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis+' (GeV)'),
                                   yaxis = dict(title=self.y_axis+' (GeV)'),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )
            
                elif self.variable_x.get(self.x_axis)[0]!='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis),
                                   yaxis = dict(title=self.y_axis+' (GeV)'),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )            
            
                elif self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]!='MASS':
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis+' (GeV)'),
                                   yaxis = dict(title=self.y_axis),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )             
            
                else:
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis),
                                   yaxis = dict(title=self.y_axis)
                                   )
                fig = go.Figure(data=data, layout=layout)
                plotly.offline.plot(fig, filename=path_to_plots+'/'+cont_plot+'_excluded.html',auto_open=False)
        return;
 
 
    def makeContinuousPlotsNonexcluded(self):
        """ Generate plots with continuous z axis variables, using non-excluded data points """
        if 'non-excluded' in self.plot_data:
            for cont_plot in self.cont_plots:
        
                if cont_plot=='chi2':
                    all_false=True
                    for chi2_value in self.data_frame_nonexcluded['chi2']:
                        if chi2_value!=False:
                            all_false=False
                    if all_false==True:
                        logger.info('No values where found for chi^2 in the non-excluded region. Skipping this plot')
                        continue
                plot_desc=self.plot_descriptions.get(cont_plot)
                cont_plot_legend=self.html_names.get(cont_plot)
                if cont_plot=='MT_max_xsec' or cont_plot=='MT_total_xsec' or cont_plot=='MT_long_xsec' or cont_plot=='MT_asym_xsec' or cont_plot=='MT_outgrid_xsec':
                    cont_plot_legend=self.html_names.get(cont_plot)+' (fb)'
                z=self.data_frame_nonexcluded[cont_plot]
                x=self.data_frame_nonexcluded[self.x_axis]
                y=self.data_frame_nonexcluded[self.y_axis]
                hover_text=self.data_frame_nonexcluded['hover_text']
             
                data = [
                    go.Scatter(
                    x=x,
                    y=y,
                    text=hover_text,
                    hoverinfo='text',
                    marker=dict(
                        size=10,
                        cmax=self.data_frame_nonexcluded[cont_plot].max(),
                        cmin=self.data_frame_nonexcluded[cont_plot].max(),
                        color=z,
                        colorbar=dict(
                            title=cont_plot_legend), 
                    colorscale='Jet'), 
                    mode='markers'  
            )
     
                    ]
     
                if self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis+' (GeV)'),
                                   yaxis = dict(title=self.y_axis+' (GeV)'),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )
            
                elif self.variable_x.get(self.x_axis)[0]!='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis),
                                   yaxis = dict(title=self.y_axis+' (GeV)'),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )            
            
                elif self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]!='MASS':
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis+' (GeV)'),
                                   yaxis = dict(title=self.y_axis),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )             
            
                else:
                    layout = go.Layout(hovermode= 'closest',
                                   title = self.plot_title,
                                   xaxis = dict(title=self.x_axis),
                                   yaxis = dict(title=self.y_axis),
                                   annotations=[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc,
                                                   xref='paper',
                                                   yref='paper'
                                                   )]
                                   )
        
                fig = go.Figure(data=data, layout=layout)
                plotly.offline.plot(fig, filename = path_to_plots+'/'+cont_plot+'_non-excluded.html', auto_open=False)
        return;
 
     
    #########Discrete_plots############
    def makeDiscretePlotsAll(self):
        """ Generate plots with discrete z axis variables, using all data points """
        if 'all' in self.plot_data:
            for disc_plot in self.disc_plots:
                plot_desc=self.plot_descriptions.get(disc_plot)
             
                disc_list=[]
                for value in self.data_frame_all[disc_plot]:
                    if value not in disc_list:
                        disc_list.append(value)
     
                fig = {
                    'data': [
                    {
                        'x': self.data_frame_all.loc[self.data_frame_all[disc_plot]==value][self.x_axis],
                        'y': self.data_frame_all.loc[self.data_frame_all[disc_plot]==value][self.y_axis],
                        'name': value, 'mode': 'markers',
                        'marker':dict(size=10),
                        'text':self.data_frame_all.loc[self.data_frame_all[disc_plot]==value]['hover_text'],
                        'hoverinfo':'text',
                 
                 
                    } for value in disc_list
                ]

            }
                
                if self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis+' (GeV)'},
                    'yaxis': {'title': self.y_axis+' (GeV)'}
                                } 
                
                elif self.variable_x.get(self.x_axis)[0]!='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis},
                    'yaxis': {'title': self.y_axis+' (GeV)'}
                                }  
                
                elif self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]!='MASS':
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis+' (GeV)'},
                    'yaxis': {'title': self.y_axis}
                                } 

                else:
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis+' (GeV)'},
                    'yaxis': {'title': self.y_axis+' (GeV)'}
                                }                 
            
            
     
                plotly.offline.plot(fig, filename = path_to_plots+'/'+disc_plot+'_all.html', auto_open=False)
        return;
 
 
    def makeDiscretePlotsExcluded(self):
        """ Generate plots with discrete z axis variables, using excluded data points """
        if 'excluded' in self.plot_data:
            for disc_plot in self.disc_plots:
                plot_desc=self.plot_descriptions.get(disc_plot)
                disc_list=[]
                for value in self.data_frame_excluded[disc_plot]:
                    if value not in disc_list:
                        disc_list.append(value)
     
                fig = {
                    'data': [
                    {
                        'x': self.data_frame_excluded.loc[self.data_frame_excluded[disc_plot]==value][self.x_axis],
                        'y': self.data_frame_excluded.loc[self.data_frame_excluded[disc_plot]==value][self.y_axis],
                        'name': value, 'mode': 'markers',
                        'marker':dict(size=10),
                        'text':self.data_frame_excluded.loc[self.data_frame_excluded[disc_plot]==value]['hover_text'],
                        'hoverinfo':'text',
                 
                 
                    } for value in disc_list
                ],
                'layout': {'title':self.plot_title,
                    'showlegend':True,
                    'hovermode':'closest',
                    'xaxis': {'title': self.x_axis},
                    'yaxis': {'title': self.y_axis}
                }
            }
                
                if self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis+' (GeV)'},
                    'yaxis': {'title': self.y_axis+' (GeV)'}
                                } 
                
                elif self.variable_x.get(self.x_axis)[0]!='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis},
                    'yaxis': {'title': self.y_axis+' (GeV)'}
                                }  
                
                elif self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]!='MASS':
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis+' (GeV)'},
                    'yaxis': {'title': self.y_axis}
                                } 

                else:
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis+' (GeV)'},
                    'yaxis': {'title': self.y_axis+' (GeV)'}
                                } 
     
                plotly.offline.plot(fig, filename = path_to_plots+'/'+disc_plot+'_excluded.html', auto_open=False)
        return;
 
    def makeDiscretePlotsNonexcluded(self):
        """ Generate plots with discrete z axis variables, using non-excluded data points """
        if 'non-excluded' in self.plot_data:
            for disc_plot in self.disc_plots:
                plot_desc=self.plot_descriptions.get(disc_plot)
                disc_list=[]
                for value in self.data_frame_nonexcluded[disc_plot]:
                    if value not in disc_list:
                        disc_list.append(value)
     
                fig = {
                    'data': [
                    {
                        'x': self.data_frame_nonexcluded.loc[self.data_frame_nonexcluded[disc_plot]==value][self.x_axis],
                        'y': self.data_frame_nonexcluded.loc[self.data_frame_nonexcluded[disc_plot]==value][self.y_axis],
                        'name': value, 'mode': 'markers',
                        'marker':dict(size=10),
                        'text':self.data_frame_nonexcluded.loc[self.data_frame_nonexcluded[disc_plot]==value]['hover_text'],
                        'hoverinfo':'text',
                 
                 
                    } for value in disc_list
                ],
                'layout': {'title':self.plot_title,
                    'showlegend':True,
                    'hovermode':'closest',
                    'xaxis': {'title': self.x_axis},
                    'yaxis': {'title': self.y_axis}
                }
            }
                
                if self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis+' (GeV)'},
                    'yaxis': {'title': self.y_axis+' (GeV)'}
                                } 
                
                elif self.variable_x.get(self.x_axis)[0]!='MASS' and self.variable_y.get(self.y_axis)[0]=='MASS':
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis},
                    'yaxis': {'title': self.y_axis+' (GeV)'}
                                }  
                
                elif self.variable_x.get(self.x_axis)[0]=='MASS' and self.variable_y.get(self.y_axis)[0]!='MASS':
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis+' (GeV)'},
                    'yaxis': {'title': self.y_axis}
                                } 

                else:
                    fig['layout'] =  {'title':self.plot_title,
                    'showlegend':True,       
                    'hovermode':'closest','annotations':[
                                           dict(
                                                   x=0.0,
                                                   y=1.05,
                                                   showarrow=False,
                                                   text=plot_desc+'  ('+self.html_names.get(disc_plot)+')',
                                                   xref='paper',
                                                   yref='paper')],
                    'xaxis': {'title': self.x_axis+' (GeV)'},
                    'yaxis': {'title': self.y_axis+' (GeV)'}
                                }                 
     
                plotly.offline.plot(fig, filename = path_to_plots+'/'+disc_plot+'_non-excluded.html', auto_open=False)
        return;



    def createIndexHtml(self,
                      filename = "index.html"):
        """
        Fills the index.html file with links to the interactive plots.
        """
    
        main_file= open(path_to_plots+'/'+filename, 'w')
        main_file.write('<html><head><font size=6>SModelS interactive plots: '+self.plot_title+'</font></head>')
        hyperlink_format = '<a href={link}>{text}</a>'
        for plot in self.plot_list:
            plot_name=plot.split('.')[0]
            plot_desc=self.plot_descriptions.get(plot_name)
            main_file.write('<p>'+'<strong>'+self.html_names.get(plot_name)+'</strong>'+': '+plot_desc+' <br>')
            for option in self.plot_data:

                
                if plot=='chi2' and os.path.isfile(path_to_plots+'/chi2_'+option+'.html')==False:
                    main_file.write('<p> <i> No &#967;<sup>2</sup> values where found in region '+option+' </i> <br>')
                    continue
                
            
                plot_link=hyperlink_format.format(link=plot_name+'_'+option+'.html', text=option)
                main_file.write(plot_link)
                main_file.write(' ')
            main_file.write('</p>')
        main_file.close()
        return True
