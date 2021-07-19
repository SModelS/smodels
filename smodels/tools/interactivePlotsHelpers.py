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
import imp
import smodels
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
        with open(smodelsFile, 'rb') as fsmodels:
            ## imports smodels file
            smodelsOut = imp.load_module("smodelsOutput",fsmodels,smodelsFile,('.py', 'rb', imp.PY_SOURCE))
            smodelsOutput = smodelsOut.smodelsOutput

    except (ImportError,AttributeError,IOError,ValueError,OSError,SyntaxError):
        logger.debug("Error loading smodels file %s. Does it contain a smodelsOutput dictionary?" %smodelsFile)

        return False

    if not isinstance(smodelsOutput,dict):
        logger.warning("smodelsOutput in file %s is not a dictionary." %smodelsFile)
        return False

    return smodelsOutput

def outputStatus(smodelsDict):
    """
    Check the smodels output status in the file, if it's -1,
    it will append 'none' to each list in the dictionary.
    """

    outputStatus = getEntry(smodelsDict, 'OutputStatus', 'file status')
    if outputStatus is False:
        raise SModelSError()

    return outputStatus

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

def getSlhaFile(smodelsOutput):
    """
    Returns the file name of the SLHA file corresponding to the output in smodelsDict
    """

    slhaFile = getEntry(smodelsOutput,'OutputStatus','input file')
    if not slhaFile:
        raise SModelSError()

    return os.path.basename(slhaFile)

def getSlhaData(slhaFile):
    """
    Uses pyslha to read the SLHA file. Return a pyslha.Doc objec, if successful.
    """

    if not os.path.exists(slhaFile):
        logger.warning("%s not found. This point will be ignored" % slhaFile )
        return False

    try:
        slhaData = pyslha.readSLHAFile(slhaFile)
    except (OSError,pyslha.ParseError) as e:
        logger.warning("Error reading SLHA file %s: %s" % (slhaFile,e) )
        return False

    return slhaData

class Filler:
    """
     A class with the functions required to fill the data dictionary to produce the plots
    """

    def __init__( self,plotter, smodelsOutput,slhaData ):
        self.data_dict=plotter.data_dict
        self.smodelsOutput=smodelsOutput
        self.slhaData=slhaData
        self.slha_hover_information=plotter.slha_hover_information
        self.ctau_hover_information=plotter.ctau_hover_information
        self.BR_hover_information=plotter.BR_hover_information
        self.min_BR=plotter.min_BR
        self.particle_names=plotter.particle_names

        return

    def truncate(self,number):
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
            factor=10**3
            truncated=math.trunc(number * factor) / factor
        return truncated


    def getExpres(self):
        """
        Extracts the Expres info from the .py output. If requested, the data
        will be appended on each corresponding list
        """
        rmax=0
        decompStatus = getEntry(self.smodelsOutput,'OutputStatus','decomposition status')

        if decompStatus == 1:
            expResList = getEntry(self.smodelsOutput,'ExptRes')
            if not expResList or not isinstance(expResList,list):
                print(self.smodelsOutput)
                raise SModelSError("Error reading ExptRes.")
            for expres in expResList:
                if 'r' in expres:

                    r = getEntry(expres,'r')
                    r = Filler.truncate(self,r)

                else:
                    r = getEntry(expres,'theory prediction (fb)')/getEntry(expres,'upper limit (fb)')
                    r = Filler.truncate(self,r)
                if r>rmax:
                    rmax = r
                    Txname = getEntry(expres,'TxNames')
                    Txname=','.join(Txname)
                    #if getEntry(expres,'chi2')==False:
                    #    chi_2=False
                    #else:
                    #    chi_2 = getEntry(expres,'chi2')
                    #    chi_2 = Filler.truncate(self,chi_2)

                    analysis_id = getEntry(expres,'AnalysisID')
        else:
            Txname = False
            #chi_2 = False
            analysis_id = False
            rmax=False

       # if 'SModelS_status' in  self.data_dict:
       #always add 'SModelS_status' to data_dict
        if rmax ==False:
            self.data_dict['SModelS_status'].append('False')
        elif rmax>1:
            self.data_dict['SModelS_status'].append('Excluded')
        else:
            self.data_dict['SModelS_status'].append('Non-excluded')




        if 'r_max' in self.data_dict.keys():
            self.data_dict['r_max'].append(rmax)

        if 'Tx' in self.data_dict.keys():
            self.data_dict['Tx'].append(Txname)

        #if 'chi2' in self.data_dict.keys():
        #    self.data_dict['chi2'].append(chi_2)
        if 'Analysis' in self.data_dict.keys():
            self.data_dict['Analysis'].append(analysis_id)
        return self.data_dict;


    def getTotalMissingXsec(self):
        """ Extracts the total crossection from missing topologies. """

        decompStatus = getEntry(self.smodelsOutput,'OutputStatus','decomposition status')

        if decompStatus >= 0:
            total_xsec = self.smodelsOutput['Total xsec for missing topologies (fb)']
            total_xsec =  Filler.truncate(self,total_xsec)


        else:
            total_xsec=False
        if 'MT_total_xsec' in self.data_dict.keys():
            self.data_dict.get('MT_total_xsec').append(total_xsec)

        return self.data_dict

    def getMaxMissingTopology(self):
        """ Extracts the missing topology with the largest cross section  """

        decompStatus = getEntry(self.smodelsOutput,'OutputStatus','decomposition status')
        if decompStatus >= 0:
            try:
                max_xsec = self.smodelsOutput['missing topologies'][0].get('element')
            except:
                max_xsec='No-missing-topologies'
        else:
            max_xsec='False'

        if 'MT_max' in self.data_dict.keys():
            self.data_dict.get('MT_max').append(max_xsec)

        return self.data_dict

    def getMaxMissingTopologyXsection(self):
        """ Extracts the cross section of the missing topology with the largest
            cross section  """

        decompStatus = getEntry(self.smodelsOutput,'OutputStatus','decomposition status')
        if decompStatus >= 0:
            try:
                max_xsec = self.smodelsOutput['missing topologies'][0].get('weight (fb)')
                max_xsec =  Filler.truncate(self,max_xsec)
            except:
                max_xsec=0

        else:
            max_xsec=False
        if 'MT_max_xsec' in self.data_dict.keys():
            self.data_dict.get('MT_max_xsec').append(max_xsec)

        return self.data_dict

    def getTotalMissingPrompt(self):
        """
        Extracts the Total cross section from missing prompt topologies
        """
        decompStatus = getEntry(self.smodelsOutput,'OutputStatus','decomposition status')
        if decompStatus >= 0:
            total_xsec = self.smodelsOutput['Total xsec for missing topologies with prompt decays (fb)']
            total_xsec =  Filler.truncate(self,total_xsec)


        else:
            total_xsec=False
        if 'MT_prompt_xsec' in self.data_dict.keys():
            self.data_dict.get('MT_prompt_xsec').append(total_xsec)
        return self.data_dict

    def getTotalMissingDisplaced(self):
        """
        Extracts the Total cross section from missing displaced topologies
        """
        decompStatus = getEntry(self.smodelsOutput,'OutputStatus','decomposition status')
        if decompStatus >= 0:
            total_xsec = self.smodelsOutput['Total xsec for missing topologies with displaced decays (fb)']
            total_xsec =  Filler.truncate(self,total_xsec)
        else:
            total_xsec=False
        if 'MT_displaced_xsec' in self.data_dict.keys():
            self.data_dict.get('MT_displaced_xsec').append(total_xsec)
        return self.data_dict

    def getOutsideGrid(self):
        """
        Extracts the outside grid info from the .py output. If requested, the
        data will be appended on each corresponding list.
        """
        decompStatus = getEntry(self.smodelsOutput,'OutputStatus','decomposition status')
        if decompStatus >= 0:
            total_xsec = self.smodelsOutput['Total xsec for topologies outside the grid (fb)']
            total_xsec =  Filler.truncate(self,total_xsec)
        else:
            total_xsec=False
        if 'MT_outgrid_xsec' in self.data_dict.keys():
            self.data_dict.get('MT_outgrid_xsec').append(total_xsec)
        return self.data_dict

    def getSlhaHoverInfo(self):
        """
        Gets the requested slha info from each slha file, to fill the hover.
        """

        for key in self.slha_hover_information.keys():
            block = self.slha_hover_information.get(key)[0]
            code_number = self.slha_hover_information.get(key)[1]
            if block=='MASS':
                self.data_dict.get(key).append(Filler.truncate(self,abs(self.slhaData.blocks[block][code_number])))
            else:
                self.data_dict.get(key).append(Filler.truncate(self,self.slhaData.blocks[block][code_number]))

        return self.data_dict

    def getCtau(self):
        """
        Computes the requested ctaus, that will go into de hover.
        """

        for key in self.ctau_hover_information.keys():
            value=self.ctau_hover_information.get(key)
            total_width=float(str(self.slhaData.decays[value]).split('total width = ')[1].split(' GeV')[0])
            if total_width == 0:
                ctau=float('inf')
            else:
                mean_lifetime=(6.582119e-16)/(total_width*1e9)
                ctau=(mean_lifetime)*(299792458)

                ctau=Filler.truncate(self,ctau)

            self.data_dict.get(key).append(ctau)

        return self.data_dict

    def openSMParticles(self):
        """Loads  SMparticles.py to parse over SM pdg-labels"""
        from smodels import installation
        sm_particle_file=f'{installation.installDirectory()}/smodels/share/models/SMparticles.py'
        with open(sm_particle_file, 'rb') as fsmparticles:
            ## imports parameter file
            self.sm_particle_names = imp.load_module("sm_particles",\
                    fsmparticles,sm_particle_file,('.py', 'rb', imp.PY_SOURCE))
        return

    def getParticleName(self,pdg):
        """ looks for the particle label in the model.py file """
        found=False

        full_list=self.particle_names.BSMList+self.sm_particle_names.SMList
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
            particle_name=str(pdg)

        return particle_name

    def getBR(self):
        """
        Gets the requested branching ratios from the slha file, that will go into de hover.
        """
        for key in self.BR_hover_information.keys():
            pdg_number=self.BR_hover_information.get(key)
            BRs=str(self.slhaData.decays[pdg_number]).split('\n')[1:]

            if self.min_BR == 'all':
                BR_top=BRs
            else:
                BR_top=[]
                for br in BRs:
                    br_val=float(br.split(' [')[0])
                    if br_val<self.min_BR:
                        break
                    BR_top.append(br)

            ##### translating pdg to particle name using particles.py

            if self.particle_names !=None:
                new_BR_top=[]

                for br in BR_top:
                #  print(br)
                    br_val=br.split(' [')[0]
                    list_daughters=br.split(' [')[1][:-1].split(',')
                    #print(list_daughters)
                    new_list_daughters=[]
                    for daughter in list_daughters:
                        #print(int(daughter))
                        new_list_daughters.append(Filler.getParticleName(self,int(daughter)))
                    new_daughters=(',').join(new_list_daughters)
                    br=(' [').join([br_val,new_daughters])+']'
                    #print(br)
                    new_BR_top.append(br)
                BR_top=new_BR_top
            BR_top=str(BR_top).split(' ')
            BR_top=''.join(BR_top)
            BR_top=BR_top[2:-2]
            self.data_dict.get(key).append(BR_top)

        return self.data_dict

    def getVariable(self,variable):
        """
        Gets the variable from the slha file.
        """
        for key in variable.keys():
            if str(key) not in self.slha_hover_information.keys():
                block=variable.get(key)[0]
                code_number=variable.get(key)[1]
                if block=='MASS':
                    self.data_dict.get(key).append(Filler.truncate(self,abs(self.slhaData.blocks[block][code_number])))
                else:
                    self.data_dict.get(key).append(Filler.truncate(self,self.slhaData.blocks[block][code_number]))

        return self.data_dict

    def getSmodelSData(self):
        """ fills data dict with smodels data """

        self.data_dict = Filler.getExpres(self)
        #self.data_dict = Filler.getMissedTopologies(self)
        #self.data_dict = Filler.getAsymmetricBranches(self)
        self.data_dict = Filler.getTotalMissingXsec(self)
        self.data_dict = Filler.getMaxMissingTopology(self)
        self.data_dict = Filler.getMaxMissingTopologyXsection(self)

        self.data_dict = Filler.getTotalMissingPrompt(self)
        self.data_dict = Filler.getTotalMissingDisplaced(self)
        self.data_dict = Filler.getOutsideGrid(self)
       # self.data_dict = Filler.getLongCascades(self)

        return self.data_dict

    def getSlhaData(self,variable_x,variable_y):
        ''' fills data dict with slha data'''

        Filler.openSMParticles(self)

        self.data_dict =  Filler.getSlhaHoverInfo(self)
        self.data_dict = Filler.getCtau(self)
        self.data_dict = Filler.getBR(self)
        #Fill with the x and y data:
        #print(list(self.variable_x.keys())[0])
        if list(variable_x.keys())[0] not in self.slha_hover_information.keys():
            self.data_dict = Filler.getVariable(self, variable_x)
        if list(variable_y.keys())[0] not in self.slha_hover_information.keys():
            self.data_dict = Filler.getVariable(self, variable_y)

        return self.data_dict


class PlotlyBackend:
    def __init__ ( self, master, path_to_plots ):
        self.data_dict=master.data_dict
        self.SModelS_hover_information=master.SModelS_hover_information
        self.slha_hover_information=master.slha_hover_information
        self.ctau_hover_information=master.ctau_hover_information
        self.BR_hover_information=master.BR_hover_information
        self.variable_x=master.variable_x
        self.variable_y=master.variable_y
        self.plot_list=master.plot_list
        self.plot_data=master.plot_data
        self.plot_title=master.plot_title
        self.path_to_plots=path_to_plots

    def makeDataFrame(self):
        """
        Transform the main dictionary in a data frame.
        """
        self.data_frame_all = pd.DataFrame(data=self.data_dict)
        self.data_frame_all.to_csv(self.path_to_plots+'/data_frame.txt', sep=' ', index=False,header=True)

        return self.data_frame_all


    def refiningVariableNames(self):
        ''' Redifining the output variable names to html format  '''

        self.html_names={'SModelS_status':'SModelS status',
                      'r_max':'r<sub>max</sub>',
                      #'chi2':' &#967;<sup>2</sup>',
                      'Tx':'T<sub>max</sub>',
                      'Analysis':'Analysis',
                      'MT_max':'MT<sub>max</sub>',
                      'MT_max_xsec':'MT<sub>max xsection</sub>.',
                      'MT_total_xsec':'MT<sub>total xsection</sub>',
                      'MT_prompt_xsec':'MT<sub>prompt xsection</sub>',
                      'MT_displaced_xsec':'MT<sub>displaced xsection</sub>',
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
            for i in range(len(self.data_frame_all.index)):
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

        #Provisinal data frame to handle changes False-->not available
        data_frame_provisional=self.data_frame_all
        data_frame_provisional=data_frame_provisional.astype(str)

        #if 'SModelS_status' in self.SModelS_hover_information:
        #Always add smodels_status to the dataframe

        data_frame_provisional.loc[data_frame_provisional['SModelS_status'] =='False', 'SModelS_status'] = 'Not-available'
        self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('SModelS_status')+': '+data_frame_provisional['SModelS_status']+'<br>'


        if 'r_max' in self.SModelS_hover_information:

            data_frame_provisional.loc[data_frame_provisional['r_max'] =='False', 'r_max'] = 'Not-available'

            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('r_max')+': '+data_frame_provisional['r_max']+'<br>'

        if 'Tx' in self.SModelS_hover_information:

            data_frame_provisional.loc[data_frame_provisional['Tx'] =='False', 'Tx'] = 'Not-available'


            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('Tx')+': '+data_frame_provisional['Tx']+'<br>'
        if 'Analysis' in self.SModelS_hover_information:

            data_frame_provisional.loc[data_frame_provisional['Analysis'] =='False', 'Analysis'] = 'Not-available'

            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('Analysis')+': '+data_frame_provisional['Analysis']+'<br>'
            
        #if 'chi2' in self.SModelS_hover_information:
#
#
#            data_frame_provisional.loc[data_frame_provisional['chi2'] =='False', 'chi2'] = 'Not-available'
#
#
#
#            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('chi2')+': '+data_frame_provisional['chi2'].astype('str')+'<br>'
        if 'MT_max' in self.SModelS_hover_information:

            data_frame_provisional.loc[data_frame_provisional['MT_max'] =='False', 'MT_max'] = 'Not-available'

            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_max')+': '+data_frame_provisional['MT_max']+'<br>'
        if 'MT_max_xsec' in self.SModelS_hover_information:

            data_frame_provisional.loc[data_frame_provisional['MT_max_xsec'] =='False', 'MT_max_xsec'] = 'Not-available'

            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_max_xsec')+': '+data_frame_provisional['MT_max_xsec'] +' [fb]'+'<br>'
        if 'MT_total_xsec' in self.SModelS_hover_information:

            data_frame_provisional.loc[data_frame_provisional['MT_total_xsec'] =='False', 'MT_total_xsec'] = 'Not-available'

            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_total_xsec')+': '+data_frame_provisional['MT_total_xsec']+' [fb]'+'<br>'

        if 'MT_prompt_xsec' in self.SModelS_hover_information:

            data_frame_provisional.loc[data_frame_provisional['MT_prompt_xsec'] =='False', 'MT_prompt_xsec'] = 'Not-available'

            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_prompt_xsec')+': '+data_frame_provisional['MT_prompt_xsec']+' [fb]'+'<br>'

        if 'MT_displaced_xsec' in self.SModelS_hover_information:

            data_frame_provisional.loc[data_frame_provisional['MT_displaced_xsec'] =='False', 'MT_displaced_xsec'] = 'Not-available'

            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_displaced_xsec')+': '+data_frame_provisional['MT_displaced_xsec']+' [fb]'+'<br>'

        if 'MT_outgrid_xsec' in self.SModelS_hover_information:

            data_frame_provisional.loc[data_frame_provisional['MT_outgrid_xsec'] =='False', 'MT_outgrid_xsec'] = 'Not-available'

            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+self.html_names.get('MT_outgrid_xsec')+': '+data_frame_provisional['MT_outgrid_xsec']+' [fb]'+'<br>'
        if 'file' in self.SModelS_hover_information:
            data_frame_provisional.loc[data_frame_provisional['file'] =='False', 'file'] = 'Not-available'
            self.data_frame_all['hover_text']=self.data_frame_all['hover_text']+'file'+': '+data_frame_provisional['file']+'<br>'
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
              #'chi2':'&#967;<sup>2</sup> value associated to the highest r-value.',
              'Tx':'Topology/ies which give the highest r-value.',
              'Analysis':'Experimental analysis that the highest r-value comes from.',
              'MT_max':'Missing topologies with the largest cross section.',
              'MT_max_xsec':'Cross section of MT_max.',
              'MT_total_xsec':'Total missing cross section.',
              'MT_prompt_xsec':'Extracts the total cross section from missing prompt topologies',
              'MT_displaced_xsec':'Extracts the total cross section from missing displaced topologies',
              'MT_outgrid_xsec':'Missing cross section outside the mass grids of the experimental results.'}
        return self.plot_descriptions;

#####continuous plots##############
    def makeContinuousPlots(self,data_frame,data_selection):
            """ Generate plots with continuous z axis variables, using all data points """
            #if 'all' in self.plot_data:
            for cont_plot in self.cont_plots:




               # if cont_plot=='chi2':
               #     all_false=True
               #     for chi2_value in data_frame['chi2']:
               #         if chi2_value!=False:
               #             all_false=False
               #     if all_false==True:
               #         logger.info('No values where found for chi^2. Skipping this plot')
               #         continue

                ####select only available values
                data_frame_noFalse=data_frame.loc[data_frame[cont_plot]!=False]


                plot_desc=self.plot_descriptions.get(cont_plot)
                cont_plot_legend=self.html_names.get(cont_plot)
                if cont_plot=='MT_max_xsec' or cont_plot=='MT_total_xsec' or cont_plot=='MT_prompt_xsec' or cont_plot=='MT_displaced_xsec' or cont_plot=='MT_outgrid_xsec':
                    cont_plot_legend=self.html_names.get(cont_plot)+' (fb)'
                z=data_frame_noFalse[cont_plot]
                x=data_frame_noFalse[self.x_axis]
                y=data_frame_noFalse[self.y_axis]
                hover_text=data_frame_noFalse['hover_text']

                data = [
                    go.Scatter(
                    x=x,
                    y=y,
                    text=hover_text,
                    hoverinfo='text',
                    mode='markers',
                    marker=dict(
                        size=10,
                        cmax=data_frame_noFalse[cont_plot].max(),
                        cmin=data_frame_noFalse[cont_plot].max(),
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
                plotly.offline.plot(fig, filename = self.path_to_plots+'/'+cont_plot+'_'+data_selection+'.html', auto_open=False)
            return;

    #########Discrete_plots############
    def makeDiscretePlots(self,data_frame,data_selection):
        """ Generate plots with discrete z axis variables, using all data points """

        for disc_plot in self.disc_plots:
            data_frame_noFalse=data_frame.loc[data_frame[disc_plot]!=False]
            plot_desc=self.plot_descriptions.get(disc_plot)

            disc_list=[]
            for value in data_frame_noFalse[disc_plot]:
                if value not in disc_list:
                    disc_list.append(value)

            fig = {
                'data': [
                {
                    'x': data_frame_noFalse.loc[data_frame_noFalse[disc_plot]==value][self.x_axis],
                    'y': data_frame_noFalse.loc[data_frame_noFalse[disc_plot]==value][self.y_axis],
                    'name': value, 'mode': 'markers',
                    'marker':dict(size=10),
                    'text':data_frame_noFalse.loc[data_frame_noFalse[disc_plot]==value]['hover_text'],
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



            plotly.offline.plot(fig, filename = self.path_to_plots+'/'+disc_plot+'_'+data_selection+'.html', auto_open=False)
        return;

    def createIndexHtml(self):
        """
        Fills the index.html file with links to the interactive plots.
        """
        main_file= open(self.path_to_plots+'/'+self.indexfile, 'w')
        main_file.write('<html><head><font size=6>SModelS interactive plots: '+self.plot_title+'</font></head>')
        hyperlink_format = '<a href={link}>{text}</a>'
        for plot in self.plot_list:
            plot_name=plot.split('.')[0]
            plot_desc=self.plot_descriptions.get(plot_name)
            main_file.write('<p>'+'<strong>'+self.html_names.get(plot_name)+'</strong>'+': '+plot_desc+' <br>')
            for option in self.plot_data:


               # if plot=='chi2' and os.path.isfile(self.path_to_plots+'/chi2_'+option+'.html')==False:
               #     main_file.write('<p> <i> No &#967;<sup>2</sup> values where found in region '+option+' </i> <br>')
               #     continue


                plot_link=hyperlink_format.format( link=plot_name+'_'+option+'.html', 
                                                   text=option )
                main_file.write(plot_link)
                main_file.write(' ')
            main_file.write('</p>')
        main_file.close()
        return True

    def makePlots(self, indexfile ):
        """
        Uses the data in self.data_dict to produce the plots.

        :parameter outFolder: Path to the output folder.
        """

        self.data_frame_all = PlotlyBackend.makeDataFrame(self)
        self.indexfile = indexfile
        self.html_names=PlotlyBackend.refiningVariableNames(self)
        self.data_frame_all = PlotlyBackend.fillHover(self)
        data_frame_excluded,data_frame_nonexcluded = PlotlyBackend.DataFrameExcludedNonexcluded(self)
        self.x_axis,self.y_axis = PlotlyBackend.GetXyAxis(self)
        self.cont_plots,self.disc_plots = PlotlyBackend.SeparateContDiscPlots(self)

        plot_descriptions=PlotlyBackend.plotDescription(self)

        if 'all' in self.plot_data:

            if self.data_frame_all.shape[0]==0:
                logger.warning('Empty data frame for dataselection:all. Skipping this plots')
                new_plot_data=[]
                for option in self.plot_data:
                    if option=='all':
                        continue
                    new_plot_data.append(option)
                self.plot_data=new_plot_data
            else:
                PlotlyBackend.makeContinuousPlots(self,self.data_frame_all,'all')
                PlotlyBackend.makeDiscretePlots(self,self.data_frame_all,'all')

        if 'excluded' in self.plot_data:
            if self.data_frame_excluded.shape[0]==0:
                logger.warning('Empty data frame for dataselection:excluded. Skipping this plots')
                new_plot_data=[]
                for option in self.plot_data:
                    if option=='excluded':
                        continue
                    new_plot_data.append(option)
                self.plot_data=new_plot_data
            else:
                PlotlyBackend.makeContinuousPlots(self,data_frame_excluded,'excluded')
                PlotlyBackend.makeDiscretePlots(self,data_frame_excluded,'excluded')

        if 'non-excluded' in self.plot_data:

            if self.data_frame_nonexcluded.shape[0]==0:
                logger.warning('Empty data frame for dataselection:non-excluded. Skipping this plots')
                new_plot_data=[]
                for option in self.plot_data:
                    if option=='non-excluded':
                        continue
                    new_plot_data.append(option)
                self.plot_data=new_plot_data
            else:
                PlotlyBackend.makeContinuousPlots(self,data_frame_nonexcluded,'non-excluded')
                PlotlyBackend.makeDiscretePlots(self,data_frame_nonexcluded,'non-excluded')
        PlotlyBackend.createIndexHtml(self)
        return True
