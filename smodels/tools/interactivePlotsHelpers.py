"""
.. module:: interactivePlotsHelpers
   :synopsis: Main functions for producing interactive plots
   
.. moduleauthor:: Humberto Reyes <humberto.reyes-gonzalez@lpsc.in2p3.fr>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

   
"""


import plotly.graph_objs as go
import plotly,imp
import pandas as pd
import os
import pyslha
from smodels.tools.smodelsLogging import logger
from smodels.theory.exceptions import SModelSTheoryError as SModelSError



def import_python_output(smodelsFile): 
    """
    Imports the smodels output from each .py file.
    """
       
    try:
        with open(smodelsFile, 'rb') as fsmodels: ## imports smodels file
            smodelsOut = imp.load_module("smodelsOutput",fsmodels,smodelsFile,('.py', 'rb', imp.PY_SOURCE))
            smodelsOutput = smodelsOut.smodelsOutput
    except:
        logger.debug("Error loading smodels file %s. Does it contain a smodelsOutput dictionary?" %smodelsFile)
        return False
    
    if not isinstance(smodelsOutput,dict):    
        logger.warning("smodelsOutput in file %s is not a dictionary." %smodelsFile)
        return False
    
    return smodelsOutput


def get_entry(inputDict,*keys):
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
        return get_entry(inputDict[key],*keys)
 
def output_status(smodelsDict):
    """
    Check the smodels output status in the file, if it's -1, 
    it will append 'none' to each list in the dictionary.
    """
    
    outputStatus = get_entry(smodelsDict, 'OutputStatus', 'file status')
    if outputStatus is False:
        raise SModelSError()

    return outputStatus

def get_slha_file(smodelsDict):
    """
    Returns the file name of the SLHA file corresponding to the output in smodelsDict
    """

    slhaFile = get_entry(smodelsDict,'OutputStatus','input file')
    if not slhaFile:
        raise SModelSError()

    return os.path.basename(slhaFile)

def get_slha_data(slhaFile):
    """
    Uses pyslha to read the SLHA file. Return a pyslha.Doc objec, if successful.
    """

    if not os.path.isfile(slhaFile):
        logger.warning("SLHA file %s not found. This point will be ignored" %slhaFile)
        return False
    
    try:
        slhaData = pyslha.readSLHAFile(slhaFile)
    except:
        logger.warning("Error reading SLHA file %s." %slhaFile)
        return False
    
    return slhaData
     
def get_expres(data_dict,smodelsOutput):
    """
    Extracts the Expres info from the .py output. If requested, the data will be appended on each corresponding list
    """   
    
    rmax=0
    
    decompStatus = get_entry(smodelsOutput,'OutputStatus','decomposition status')
    if decompStatus == 1:
        expResList = get_entry(smodelsOutput,'ExptRes')
        if not expResList or not isinstance(expResList,list):
            print(smodelsOutput)
            raise SModelSError("Error reading ExptRes.")    
        for expres in expResList:
            if 'r' in expres:
                r = get_entry(expres,'r')
            else:
                r = get_entry(expres,'theory prediction (fb)')/get_entry(expres,'upper limit (fb)')
            if r>rmax:
                rmax = r
                Txname = get_entry(expres,'TxNames')
                if get_entry(expres,'chi2')==False:
                    chi_2=None
                else:    
                    chi_2 = get_entry(expres,'chi2')

                analysis_id = get_entry(expres,'AnalysisID')
    else:
        Txname = 'None'
        chi_2 = None
        analysis_id = 'None'
        rmax=None
 
    if 'SmodelS_excluded' in  data_dict:
        if rmax ==None:
            data_dict['SmodelS_excluded'].append('None')
        elif rmax>1:
            data_dict['SmodelS_excluded'].append('Excluded')  
        else:
            data_dict['SmodelS_excluded'].append('Non-excluded')  

    if 'r_max' in data_dict.keys():
        data_dict['r_max'].append(rmax)
        
    if 'Tx' in data_dict.keys(): 
        data_dict['Tx'].append(str(Txname))
        
    if 'chi2' in data_dict.keys(): 
        data_dict['chi2'].append(chi_2) 
    if 'Analysis' in data_dict.keys(): 
        data_dict['Analysis'].append(analysis_id)   
    return data_dict;   
  
def get_missed_topologies(data_dict,smodelsOuptut):
    """
    Extracts the Missed topologies info from the .py output. If requested, the data will be appended on each corresponding list
    """
    decompStatus = get_entry(smodelsOuptut,'OutputStatus','decomposition status')
    missedtopo_max_xsec=0
    missedtopo_total_xsec=0
    mt_max = 'None'
    if decompStatus >= 0:
        for missed_topo in get_entry(smodelsOuptut, 'Missed Topologies'):
            missedtopo_xsec=missed_topo.get('weight (fb)')
            missedtopo_total_xsec=missedtopo_total_xsec+missedtopo_xsec
            if missedtopo_xsec>missedtopo_max_xsec:
                missedtopo_max_xsec=missedtopo_xsec
                mt_max=missed_topo.get('element')       
    else:
        missedtopo_max_xsec=None
        missedtopo_total_xsec=None
        missed_topo = 'None'
    if 'MT_max' in data_dict.keys():
        data_dict.get('MT_max').append(str(mt_max))
    if 'MT_max_xsec' in data_dict.keys():
        data_dict.get('MT_max_xsec').append(missedtopo_max_xsec)    
    if 'MT_total_xsec' in data_dict.keys():
        data_dict.get('MT_total_xsec').append(missedtopo_total_xsec)     
        
    return data_dict           
   
def get_long_cascades(data_dict,smodelsOutput):
    """
    Extracts the Long cascade info from the .py output. If requested, the data will be appended on each corresponding list
    """
    decompStatus = get_entry(smodelsOutput,'OutputStatus','decomposition status')    
    if decompStatus >= 0:
        long_cascade_total_xsec=0
        for long_cascade in smodelsOutput.get('Long Cascades'):
            long_cascade_xsec=long_cascade.get('weight (fb)')
            long_cascade_total_xsec=long_cascade_total_xsec+long_cascade_xsec
    else:
        long_cascade_total_xsec=None
    if 'MT_long_xsec' in data_dict.keys():
        data_dict.get('MT_long_xsec').append(long_cascade_total_xsec)
             
    return data_dict
             
def get_asymmetric_branches(data_dict,smodelsOutput):
    """
    Extracts the asymmetric branches info from the .py output. If requested, the data will be appended on each corresponding list
    """
    decompStatus = get_entry(smodelsOutput,'OutputStatus','decomposition status')
    if decompStatus >= 0:
        asymmetric_branch_total_xsec = sum([asym_br['weight (fb)'] for asym_br in smodelsOutput['Asymmetric Branches']])
    else:
        asymmetric_branch_total_xsec=None    
    if 'MT_asym_xsec' in data_dict.keys():
        data_dict.get('MT_asym_xsec').append(asymmetric_branch_total_xsec) 
    return data_dict  
    
def get_outside_grid(data_dict,smodelsOutput):
    """
    Extracts the outside grid info from the .py output. If requested, the data will be appended on each corresponding list.
    """    
    decompStatus = get_entry(smodelsOutput,'OutputStatus','decomposition status')
    outside_grid_total_xsec=0
    if decompStatus >= 0:
        for outside_grid in smodelsOutput.get('Outside Grid'):
            outside_grid_xsec=outside_grid.get('weight (fb)')
            outside_grid_total_xsec = outside_grid_total_xsec+outside_grid_xsec
    else:
        outside_grid_total_xsec = None
    if 'MT_outgrid_xsec' in data_dict.keys():
        data_dict.get('MT_outgrid_xsec').append(outside_grid_total_xsec)  
    return data_dict   

def get_slha_hover_info(data_dict,slhaData,slha_hover_information):
        """
        Gets the requested slha info from eachh slha file, to fill the hover.
        """
          
        for key in slha_hover_information.keys():
            block = slha_hover_information.get(key)[0]
            code_number = slha_hover_information.get(key)[1]
            data_dict.get(key).append(slhaData.blocks[block][code_number])
            
        return data_dict
        
def get_ctau(data_dict,slhaData,ctau_hover_information):    
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
            
        data_dict.get(key).append(ctau)
        
    return data_dict 

              
              
def get_BR(data_dict,slhaData,BR_hover_information,BR_get_top):
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

def get_variable(data_dict,slhaData,slha_hover_information,variable):
    """
    Gets the variable from the slha file.
    """  
    
    for key in variable.keys():
        if key not in slha_hover_information:
            block=variable.get(key)[0]
            code_number=variable.get(key)[1]
            data_dict.get(key).append(slhaData.blocks[block][code_number])
            
    return data_dict
    

def make_data_frame(data_dict):
    """
    Transform the main dictionary in a data frame.
    """
    
    data_frame_all = pd.DataFrame(data=data_dict)
    
    return data_frame_all
  
 
def fill_hover(data_frame_all,smodels_hover_information,
               slha_hover_information,ctau_hover_information,BR_hover_information):
    """ Generates the text of the hover, according to users's requests. """
    data_frame_all['hover_text']=''
    for column in slha_hover_information:
        data_frame_all['hover_text'] = data_frame_all['hover_text']+column+': '+data_frame_all[column].astype('str')+'<br>' 
    for column in ctau_hover_information:
        data_frame_all['hover_text'] = data_frame_all['hover_text']+column+': '+data_frame_all[column].astype('str')+'<br>' 
    data_frame_br=pd.DataFrame()
    for column in BR_hover_information:
        i=0
        data_frame_br[column] = data_frame_all[column]
        for i in range(len(data_frame_all.index)):
            if data_frame_all[column][i]:
                data_frame_br[column][i] = data_frame_all[column][i].split("','")
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
        data_frame_all['hover_text']=data_frame_all['hover_text']+column+': '+data_frame_br[column].astype('str')+'<br>'     
    for column in smodels_hover_information:
        data_frame_all['hover_text']=data_frame_all['hover_text']+column+': '+data_frame_all[column].astype('str')+'<br>'
    return data_frame_all;
 
 
  
 
def data_frame_excluded_nonexcluded(data_frame_all):
    """ Generate sub data frames for excluded and non-excluded points """
    data_frame_excluded=data_frame_all.loc[data_frame_all['SmodelS_excluded']=='Excluded']
    data_frame_nonexcluded=data_frame_all.loc[data_frame_all['SmodelS_excluded']=='Non-excluded']
    return data_frame_excluded, data_frame_nonexcluded;
 
def get_xy_axis(variable_x,variable_y):
    """ Retrieves the names of the x and y axis varaibles. """
    x_axis=list(variable_x.keys())[0]
    y_axis=list(variable_y.keys())[0]
    return x_axis,y_axis;
     
 
def separate_cont_disc_plots(plot_list,data_dict):
    ''' Generate sub lists of plots with discrete and conitnuous z axis variables. '''
    cont_plots=[]
    disc_plots=[]
    discrete_options=['Tx','Analysis','MT_max','file','SmodelS_excluded']
    for plot in plot_list:
        if plot in discrete_options:
            disc_plots.append(plot)
        else:
            cont_plots.append(plot)
    return cont_plots, disc_plots;   
 
 
#####continuous plots##############
def make_continuous_plots_all(cont_plots,x_axis,y_axis,path_to_plots,data_frame_all,plot_data,plot_title):
    """ Generate plots with continuous z axis variables, using all data points """
    if 'all' in plot_data: 
        for cont_plot in cont_plots:
            z=data_frame_all[cont_plot]
            y=data_frame_all[x_axis]
            x=data_frame_all[y_axis]
            hover_text=data_frame_all['hover_text']
            
            data = [
                go.Scatter(
                    x=x,
                    y=y,
                    text=hover_text,
                    hoverinfo='text',
                    marker=dict(
                        size=6,
                        cmax=data_frame_all[cont_plot].max(),
                        cmin=data_frame_all[cont_plot].max(),
                        color=z,
                        colorbar=dict(
                            title=cont_plot), 
                    colorscale='Jet'), 
                    mode='markers'  
            )
                    ]
     
            layout = go.Layout(hovermode= 'closest',
                title = plot_title,
                xaxis = dict(title=x_axis),
                yaxis = dict(title=y_axis)
            )
            fig = go.Figure(data=data, layout=layout)
            plotly.offline.plot(fig, filename = path_to_plots+'/'+cont_plot+'_all.html', auto_open=False)
    return;
 
 
def make_continuous_plots_excluded(cont_plots,x_axis,y_axis,path_to_plots,data_frame_excluded,plot_data,plot_title):
    """ Generate plots with continuous z axis variables, using excluded data points """
    if 'excluded' in plot_data: 
        for cont_plot in cont_plots:
         
            z=data_frame_excluded[cont_plot]
            y=data_frame_excluded[x_axis]
            x=data_frame_excluded[y_axis]
            hover_text=data_frame_excluded['hover_text']
             
             
            data = [
                go.Scatter(
                    x=x,
                    y=y,
                    text=hover_text,
                    hoverinfo='text',
                    marker=dict(
                        size=6,
                        cmax=data_frame_excluded[cont_plot].max(),
                        cmin=data_frame_excluded[cont_plot].max(),
                        color=z,
                        colorbar=dict(
                            title=cont_plot), 
                    colorscale='Jet'), 
                    mode='markers'  
            )
     
                    ]
     
            layout = go.Layout(hovermode= 'closest',
                title=plot_title,
                xaxis = dict(title=x_axis),
                yaxis = dict(title=y_axis)
            )
            fig = go.Figure(data=data, layout=layout)
            plotly.offline.plot(fig, filename=path_to_plots+'/'+cont_plot+'_excluded.html',auto_open=False)
    return;
 
 
def make_continuous_plots_nonexcluded(cont_plots,x_axis,y_axis,path_to_plots,data_frame_nonexcluded,plot_data,plot_title):
    """ Generate plots with continuous z axis variables, using non-excluded data points """
    if 'non-excluded' in plot_data: 
        for cont_plot in cont_plots:
         
            z=data_frame_nonexcluded[cont_plot]
            y=data_frame_nonexcluded[x_axis]
            x=data_frame_nonexcluded[y_axis]
            hover_text=data_frame_nonexcluded['hover_text'] 
             
            data = [
                go.Scatter(
                    x=x,
                    y=y,
                    text=hover_text,
                    hoverinfo='text',
                    marker=dict(
                        size=6,
                        cmax=data_frame_nonexcluded[cont_plot].max(),
                        cmin=data_frame_nonexcluded[cont_plot].max(),
                        color=z,
                        colorbar=dict(
                            title=cont_plot), 
                    colorscale='Jet'), 
                    mode='markers'  
            )
     
                    ]
     
            layout = go.Layout(hovermode= 'closest',
                title=plot_title,
                xaxis = dict(title=x_axis),
                yaxis = dict(title=y_axis)
            )
            fig = go.Figure(data=data, layout=layout)
            plotly.offline.plot(fig, filename = path_to_plots+'/'+cont_plot+'_non-excluded.html', auto_open=False)
    return;
 
     
    #########Discrete_plots############
def make_discrete_plots_all(disc_plots,x_axis,y_axis,path_to_plots,data_frame_all,plot_data,plot_title):
    """ Generate plots with discrete z axis variables, using all data points """
    if 'all' in plot_data:
        for disc_plot in disc_plots:
             
            disc_list=[]
            for value in data_frame_all[disc_plot]:   
                if value not in disc_list:
                    disc_list.append(value)  
     
            fig = {
                'data': [
                    {
                        'x': data_frame_all.loc[data_frame_all[disc_plot]==value][x_axis],
                        'y': data_frame_all.loc[data_frame_all[disc_plot]==value][y_axis],
                        'name': value, 'mode': 'markers',
                        'marker':dict(size=6),
                        'text':data_frame_all.loc[data_frame_all[disc_plot]==value]['hover_text'],
                        'hoverinfo':'text',
                 
                 
                    } for value in disc_list
                ],
                'layout': {'title':plot_title,
                    'showlegend':True,       
                    'hovermode':'closest',
                    'xaxis': {'title': x_axis},
                    'yaxis': {'title': y_axis}
                }
            }
     
            plotly.offline.plot(fig, filename = path_to_plots+'/'+disc_plot+'_all.html', auto_open=False)
    return;
 
 
def make_discrete_plots_excluded(disc_plots,x_axis,y_axis,path_to_plots,data_frame_excluded,plot_data,plot_title):
    """ Generate plots with discrete z axis variables, using excluded data points """
    if 'excluded' in plot_data:
        for disc_plot in disc_plots:
             
            disc_list=[]
            for value in data_frame_excluded[disc_plot]:   
                if value not in disc_list:
                    disc_list.append(value)  
     
            fig = {
                'data': [
                    {
                        'x': data_frame_excluded.loc[data_frame_excluded[disc_plot]==value][x_axis],
                        'y': data_frame_excluded.loc[data_frame_excluded[disc_plot]==value][y_axis],
                        'name': value, 'mode': 'markers',
                        'marker':dict(size=6),
                        'text':data_frame_excluded.loc[data_frame_excluded[disc_plot]==value]['hover_text'],
                        'hoverinfo':'text',
                 
                 
                    } for value in disc_list
                ],
                'layout': {'title':plot_title,
                    'showlegend':True,
                    'hovermode':'closest',
                    'xaxis': {'title': x_axis},
                    'yaxis': {'title': y_axis}
                }
            }
     
            plotly.offline.plot(fig, filename = path_to_plots+'/'+disc_plot+'_excluded.html', auto_open=False)
    return;
 
def make_discrete_plots_nonexcluded(disc_plots,x_axis,y_axis,path_to_plots,data_frame_nonexcluded,plot_data,plot_title):
    """ Generate plots with discrete z axis variables, using non-excluded data points """
    if 'non-excluded' in plot_data:
        for disc_plot in disc_plots:
             
            disc_list=[]
            for value in data_frame_nonexcluded[disc_plot]:   
                if value not in disc_list:
                    disc_list.append(value)  
     
            fig = {
                'data': [
                    {
                        'x': data_frame_nonexcluded.loc[data_frame_nonexcluded[disc_plot]==value][x_axis],
                        'y': data_frame_nonexcluded.loc[data_frame_nonexcluded[disc_plot]==value][y_axis],
                        'name': value, 'mode': 'markers',
                        'marker':dict(size=6),
                        'text':data_frame_nonexcluded.loc[data_frame_nonexcluded[disc_plot]==value]['hover_text'],
                        'hoverinfo':'text',
                 
                 
                    } for value in disc_list
                ],
                'layout': {'title':plot_title,
                    'showlegend':True,
                    'hovermode':'closest',
                    'xaxis': {'title': x_axis},
                    'yaxis': {'title': y_axis}
                }
            }
     
            plotly.offline.plot(fig, filename = path_to_plots+'/'+disc_plot+'_non-excluded.html', auto_open=False)
    return;



def create_main_html(path_to_plots,plot_data,plot_list):
    """
    Fills the main.html file with links to the interactive plots.
    """
    
    main_file= open(path_to_plots+'/main.html', 'w')
    main_file.write('<html><head><font size=6>Smodels interactive plots.</font></head>')
    hyperlink_format = '<a href={link}>{text}</a>' 
    for plot in plot_list:
        plot_name=plot.split('.')[0]
        main_file.write('<p>'+plot_name+': ')
        for option in plot_data:   
            plot_link=hyperlink_format.format(link=plot_name+'_'+option+'.html', text=option) 
            main_file.write(plot_link)
            main_file.write(' ')  
        main_file.write('</p>')
    main_file.close()
    return True 
