#!/usr/bin/env python3

"""
.. module:: xsecResummino
   :synopsis: Computation of reference ("theory") production cross sections.

.. moduleauthor:: Théo Reymermier <theo.reymermier@gmail.com>

"""

from __future__ import print_function
import sys
import os, copy
current = os.getcwd()
sys.path.append(current)


from smodels import installation
from smodels.tools.physicsUnits import pb, TeV, GeV
from smodels.theory import crossSection
from smodels.theory.crossSection import LO, NLO, NLL
from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import subprocess
from concurrent.futures import ProcessPoolExecutor
from smodels.tools.xsecBase import XSecBase, ArgsStandardizer
import tempfile
import argparse
import pyslha
import shutil
import requests
import tarfile
from itertools import combinations

class XSecResummino(XSecBase):
    """ cross section computer class (for resummino), what else? """
    def __init__ ( self, maxOrder,slha_folder_name,sqrt = 13,ncpu=1, maycompile=True, type_writting = None, verbosity = '', json = None, particles = None, xsec_limit = None):
        """
        :param maxOrder: maximum order to compute the cross section, given as an integer
                    if maxOrder == LO, compute only LO resummino xsecs
                    if maxOrder == NLO, compute NLO resummino xsecs
                    if maxOrder == NLL, compute NLO+NLL resummino xsecs
        :param nevents: number of events for pythia run
        :param pythiaVersion: pythia6 or pythia8 (integer)
        :param sqrt: Center of mass energy to consider for the cross section calculation
        :param xsec_limit: Value below which if mode == "check", cross section at NLO order are not calculated
        :param type: If all, put all the order in the slha file, if highest, only the hightest order.
        :param json: Path to the json file with all the relevant informations concerning the resummino calculation
        :param resummino_bin: Path to resummino executable
        :param input_file_original: Path to the template input file of resummino
        :param ncpu: Number of cpu used in parallel for the calculation
        :param verbosity: Type of informations given in the logger file
        """
        self.pwd = installation.installDirectory()
        self.tooldir = os.path.join(self.pwd,"smodels","lib","resummino")
        self.getVersion()
        self.resummino_bin = os.path.join(self.tooldir,f"resummino-{self.version}","bin","resummino")
        self.input_file_original = os.path.join(self.pwd,"smodels","etc","input_resummino.in")
        if json == None:
            self.json_resummino = os.path.join(self.pwd,"smodels","etc","resummino.py")
        else:
            if os.path.isabs(json):
                self.json_resummino = json
            else:
                self.json_resummino = os.path.join(self.pwd, json)
        if particles != None :
            self.particles = list(particles)
        else:
            self.particles = None
        self.slha_folder_name = slha_folder_name
        self.maxOrder = maxOrder
        self.countNoXSecs = 0
        self.countNoNLOXSecs = 0
        self.maycompile = maycompile
        self.ncpu = ncpu
        self.sqrts = sqrt
        self.type_writting = type_writting
        self.verbosity = verbosity
        self.sqrt = sqrt
        
        
        logger.debug('installation directory is '+self.pwd)

        if xsec_limit == None:
            try:    
                with open(self.json_resummino, "r") as f:
                    data = eval(f.read())
                                    
                self.xsec_limit = data["xsec_limit"]
            except KeyError:
                self.xsec_limit = 0.00001
        else:
            self.xsec_limit = xsec_limit
            
    def checkInstallation(self, compile : bool = True ) -> bool:
        """ check if resummino is already compiled.
        :param compile: if true, then attempt a compilation if not installed
        :returns: true if we have a compiled executable
        """
        if os.path.exists ( self.resummino_bin ):
            return True
        if not compile:
            return False
        logger.info ( "Resummino not compiled; trying to build now (this might take a while)" )
        cmd = f"cd {self.tooldir}; make"
        o = subprocess.getoutput ( cmd )
        ret = os.path.exists ( self.resummino_bin )
        if not ret:
            logger.error ( f"attempt at compiling resummino failed:\n{o}" )
        return ret

    def getVersion( self ):
        """ retrieve the version from version_path,
            set self.version
            if it doesnt exist, set to default of 3.1.2 """
        version_path = os.path.join(self.tooldir, 'versions.txt')
        self.version = "?.?.?"
        if os.path.exists(version_path):
            with open(version_path, "r") as f:
                lines = f.readlines()
            for line in lines:
                if line.startswith("#"):
                    continue
                if "resummino_version" in line:
                    self.version = line.split(sep ="=")[1].strip()
        
    def calculate_one_slha(self, particles,input_file, slha_file, output_file, 
            num_try, order, log):
        """
        log file management and launch resummino command. Prepare also the
        cross section list to write the cross section onto the slha file.
        """
        with open(log, 'a') as f:
            f.write(f'{particles} cross-sections written in {slha_file}\n')
        if particles == None:
            with open(slha_file, 'a') as f:
                f.write(' #no_cross-section\n')

        #Use to check the #no_cross_section stuff 
        self.are_crosssection(slha_file, order)
        if particles == None:
            return
        Xsections = crossSection.XSectionList()
        
        for particle_pair in particles:
            self.launch_resummino(input_file, slha_file, output_file, particle_pair[0], particle_pair[1], num_try, order, Xsections, log)
        
        resummino_version = f"Resumminov{self.version}"
        nxsecs = self.addXSecToFile(Xsections, slha_file, comment = f"[pb], {resummino_version}")
        
           
    def launch_command(self,resummino_bin,input_file, output_file, order):
        """
        use resummino at the order asked by the user (order variable).
        """
        if order == 0:
            command = f"{resummino_bin} {input_file} --lo"
        if order == 1:
            command = f"{resummino_bin} {input_file} --nlo"
        if order == 2:
            command = f"{resummino_bin} {input_file}"

        with open(output_file, 'w') as f:
            if self.verbosity == "debug":
                subprocess.run(command, shell=True, stdout=f,stderr=os.sys.stderr, text=True)
            else:
                with open("/dev/null", "w") as errorhandle:
                    subprocess.run(command, shell=True, stdout=f,stderr=errorhandle, text=True)


    def launch_resummino(self, input_file, slha_file, output_file, particle_1, particle_2, num_try, order, Xsections, log):
        """
        Check everything before launching resummino.
        """
        
        #modify_slha_file(input_file, input_file, slha_file)
        self.modify_outgoing_particles(input_file, input_file, particle_1, particle_2)
        
        already_written_channel = self.find_channels(slha_file)

        if self.sqrt > 10:
            _ = str(self.sqrt*10**(-1))+'0E+04'
        else:
            _ = str(self.sqrt)+'.00E+03'
        
        already_written_channel_set = [({x,y},z,w) for (x,y), z,w in already_written_channel]
        
        logger.debug(f"channel, order and cross section {str(particle_1)} {str(particle_2)} {str(order)} {str(_)}" )
        logger.debug('the already written channels are '+ str(already_written_channel))
        if (((particle_1, particle_2), _, order)) in already_written_channel:
            return
        
        if self.mode == "check":
            self.launch_command(self.resummino_bin, input_file, output_file, 0)
            infos = self.search_in_output(output_file)
            infos = infos[0].split(" ")[2][1:]

            logger.debug(str((particle_1,particle_2))+" cross section is "+str(infos)+ " pb at LO order")
            logger.debug("Is "+str((particle_1,particle_2))+ " cross section above the limit ? "+str( float(infos)>self.xsec_limit))
            logger.debug("cross section is "+str(infos)+ " pb at LO order")
            logger.debug("Is cross section above the limit ? "+str( float(infos)>self.xsec_limit))
            if (float(infos))>(self.xsec_limit):
                logger.debug('num try is '+ str(num_try)+ ' for the moment')
                self.launch_command(self.resummino_bin, input_file, output_file, order)
                if num_try == 0:
                    hist = self.write_in_slha(output_file, slha_file, order, particle_1, particle_2, self.type_writting, Xsections, log)
            else:
                hist = self.write_in_slha(output_file, slha_file, 0, particle_1, particle_2, self.type_writting, Xsections, log)
            return
        
        if self.mode == "all":
            self.launch_command(self.resummino_bin, input_file, output_file, order)
            if num_try == 0:
                hist = self.write_in_slha(output_file, slha_file, order, particle_1, particle_2, self.type_writting, Xsections, log)
            return
        
        #:hist: variable to check if there is a problem in the cross section (strange value or no value)
        hist = 0

        if num_try == 0 and self.mode != "check":
            self.launch_command(self.resummino_bin, input_file, output_file, order)

        #here we write in the slha file.
            hist = self.write_in_slha(output_file, slha_file, order, particle_1, particle_2, self.type_writting, Xsections, log)

        #we check if we have written too much cross section
        self.are_crosssection(slha_file, order)

        #if there is an error, we inform the user and launch again resummino
        if hist == 1:
            self.warning("error in resummino, trying again")
            num_try = 0
            self.modify_outgoing_particles(input_file, input_file, particle_1, particle_2)
            self.launch_resummino(input_file, slha_file, output_file, particle_1, particle_2, num_try, order, Xsections, log)


    def search_in_output(self, output_file):
        """
        Search in the .out files of resummino (in tempfiles) to get the cross section asked by the users, 
        then extract the LO,NLO and NLL+NLO.
        If you want to get the incertainties given by resummino, you have everything here in LO, NLO and NLL.
        """
        Infos = []
        with open(output_file, 'r') as f:
            data = f.readlines()
        for i in range(len(data)):
            if "Results:" in data[i]:
                LO = data[i+1][:-1] #[:-1] to get rid of the \n
                NLO = data[i+2][:-1]
                NLL = data[i+3][:-1]
                Infos.append((LO,NLO,NLL))
        if len(Infos) == 0:
            raise RuntimeError("Please check your slha file and that you have resummino correctly installed with the install.sh script in the lib/resummino folder")
        return Infos[0]

    def create_xsection(self, result, particle_1, particle_2, order, Xsections):
        """
        Create cross section list filled with cross section objects, corresponding to all the channels calculated.
        """
        if type(result) == list:
            for i in range(order+1):
                Xsection = crossSection.XSection()
        
                Xsection.value = float(result[i]) * pb
                Xsection.info.order = i
                Xsection.info.sqrts = float(self.sqrt) * TeV
                Xsection.info.label = "WAOUW"
                Xsection.pid = (particle_1, particle_2)
                Xsections.add(Xsection)
            return

        Xsection = crossSection.XSection()
        
        Xsection.value = float(result) * pb
        Xsection.info.order = order
        Xsection.info.sqrts = float(self.sqrt) * TeV
        Xsection.info.label = "WAOUW"
        Xsection.pid = (particle_1, particle_2)
        Xsections.add(Xsection)
        return
            
    def write_in_slha(self, output_file, slha_file, order, particle_1, particle_2, type_writing, Xsections, log):
        """
        Organize here the way cross sections are written into the file (highest,
        all) and then create cross_section object to let smodels take
        care of the writing itself with the create_xsection method.
        """
        results = self.search_in_output(output_file)
        if type_writing == 'highest':
            if order == 0:
                result = results[0].split(" ")[2][1:]
            elif order == 1:
                result = results[1].split(" ")[2][1:]
            elif order == 2:
                result = results[2].split(" ")[2][1:]
            logger.debug('the highest result is'+ str(result) +' pb')
        elif type_writing == "all":
            result = [results[0].split(" ")[2][1:], results[1].split(" ")[2][1:], results[2].split(" ")[2][1:]]
        elif type_writing == None:
            result = [results[0].split(" ")[2][1:], results[1].split(" ")[2][1:], results[2].split(" ")[2][1:]]
            
            _ = ['LO', 'NLO', 'NLO+NLL']
            for i in range(self.maxOrder+1):
                logger.info(f"Cross section is {result[i]} pb for ({particle_1},{particle_2}) channel at {_[i]} in perturbation theory")
        else :
            result_list = [results[0].split(" ")[2][1:], results[1].split(" ")[2][1:], results[2].split(" ")[2][1:]]
            result = 0
            for i in results:
                _ = i.split(" ")[2][1:]
                _ = float(_)
                if _ != 0:
                    result = i.split(" ")[2][1:]     
            return result
        
        if order == 1:
            if float(results[1].split(" ")[2][1:])>2*float(results[0].split(" ")[2][1:]) or float(results[0].split(" ")[2][1:])> 2* float(results[1].split(" ")[2][1:]):
                with open(log, 'a') as f:
                    f.write(f"to much change between LO and NLO for {particle_1} and {particle_2} with {slha_file}")
                return 1
        if order == 2:
            if float(results[2].split(" ")[2][1:])>2*float(results[0].split(" ")[2][1:]) or float(results[0].split(" ")[2][1:])> 2* float(results[1].split(" ")[2][1:]):
                with open(log, 'a') as f:
                    f.write(f"to much change between LO and NLL+NLO for {particle_1} and {particle_2} with {slha_file}")
                return 1
        
        if type_writing is not None:
            self.create_xsection(result, particle_1, particle_2, order, Xsections)
            
        return 0

    def extract_m1_m2_mu(self, file_path : os.PathLike ) -> dict:
        """_summary_
        
        function to extract the breaking term of the electrowikino part (SUSY) in
        an slha file.
        Args:
            file_path (_string_): _path of the slha file_

        Returns: dictionary of:
            _int_: _M1 breaking term in SUSY models_
            _int_: _M2 braking term in SUSY models_
            _int_: _mu breaking term in SUSY models_
        """
        data = pyslha.read(file_path)

        m1 = data.blocks['EXTPAR'][1]
        m2 = data.blocks['EXTPAR'][2]
        mu = data.blocks['EXTPAR'][23]

        result = {
            'M_1(MX)': m1,
            'M_2(MX)': m2,
            'mu(MX)' : mu
        }

        return m1,m2,mu

    def extract_N1_N2_C1(self, file_path):
        """_summary_
        
        function to extract the breaking term of the electrowikino part (SUSY) in
        an slha file.
        Args:
            file_path (_string_): _path of the slha file_

        Returns:
            _float_: _N1 Mass of the neutralino 1
            _float_: _N2 Mass of the neutralino 2
            _float_: _C1 Mass of the chargino 1
            _float_: _C2 Mass of the chargnino 2
        """
        data = pyslha.read(file_path)

        N1 = data.blocks['MASS'][1000022]
        N2 = data.blocks['MASS'][1000023]
        C1 = data.blocks['MASS'][1000024]
        C2 = data.blocks['MASS'][1000037]
        return N1,N2,C1, C2

    def are_crosssection(self, slha_file, order):
        """
        check if the cross sections are already written, and remove the
        cross section written twice.
        """
        with open(slha_file, 'r') as f:
            data = f.readlines()
        test = True
        if data[-1]== " #no_cross-section" or data[-1]==" #no_cross-section\n":
            while test == True:
                if data[-2]== " #no_cross-section" or data[-2]==" #no_cross-section\n" or data[-3]==" #no_cross-section\n":
                    data.pop()
                    logger.info(f'remove from {slha_file} a #no_cross-section"')
                else:
                    test = False
            with open(slha_file, 'w') as f:
                for _ in data:
                    f.write(_)
                return
        channels = {}
        to_delete = []

        for i in range(len(data)):
            line = data[i]
            if line.startswith("XSECTION"):
                _ = [x for x in line.split(" ") if x != '']
                channel = (int(_[5]), int(_[6]),_[1], data[i+1].split(" ")[4])
                
                if channel in channels:
                    start = channels[channel]
                    end = start+1
                    while not data[end].startswith('XSECTION'):
                        end+=1
                    #end = start+order+3
                    to_delete.extend(range(start, end))
                channels[channel] = i
        lines = [line for i, line in enumerate(data) if i not in to_delete]

        with open(slha_file, 'w') as f:
            f.writelines(lines)
        
        
    def find_channels(self, slha_file):
        
        with open(slha_file, 'r') as f:
            data = f.readlines()
            
        channels = []
        for i in range(len(data)):
            line = data[i]
            if line.startswith("XSECTION"):
                _ = [x for x in line.split(" ") if x != '']
                sqrt = _[2]
                channel = (int(_[5]), int(_[6]))
                order = data[i+1].split(" ")[4]
                order = int(order)
                channels.append((channel, sqrt, order))
        return channels

    def create_routine_files(self, order, slha_folder_name):
        """
        Prepare all the paths and everything before turning into parallel task.
        resumino.py is called here to avoid multi-tasking on one file. Create also tempfile to stock all data needed
        by resummino.  
        """

        #You haven't seen anything.
        output_file = "output_file.txt"

        #List of the differents files on which resummino will be launch
        Liste_slha = []
        Liste_resummino_in = []
        Liste_output_file = []
        Liste_particles = []
        Liste = []
        self.resummino_in = tempfile.mkdtemp()
        self.resummino_out = tempfile.mkdtemp()
        self.resummino_log = tempfile.mkdtemp()
        resummino_in = '/tmp/resummino/resummino_in'
        resummino_out = '/tmp/resummino/resummino_out'
        resummino_log = '/tmp/resummino/resummino_log'
        #Check if the folder already exists
        if not os.path.exists(self.resummino_in):
            os.mkdir(self.resummino_in)
        if not os.path.exists(self.resummino_out):
            os.mkdir(self.resummino_out)
        if not os.path.exists(self.resummino_log):
            os.mkdir(self.resummino_log)
        
        #always use absolute path
        slha_folder = os.path.join( os.getcwd(), slha_folder_name)  
        # slha_folder = slha_folder_name
        #Check if the input is a file or a directory
        if not os.path.isfile(slha_folder):
            liste_slha = os.listdir(slha_folder_name)
        else:
            name = os.path.normpath(slha_folder_name).split(os.path.sep)[-1]
            liste_slha = [slha_folder_name]
            
    

        

        '''
        a = number of files
        b = number of files with already a cross section
        c = number of files ignored
        '''
        a,b,c = 0,0,0

        #Loop on all the slha file, in a random order (os.listdir)
        for slha in liste_slha:
            
            if not os.path.isfile(slha_folder):
                slha_path = os.path.join(slha_folder,slha)
            else:
                slha_path = os.path.join(current, slha)
            #Variable used to avoid strange result (huge difference between LO and NLO result)
            num_try = 0
            with open(slha_path, 'r') as f:
                data = f.readlines()
                a+=1
                
            if "SModelSv" in data[-1]:
                b+=1
                # We increase this variable by 1, so if it is > 0 we do not redo the calculation 
                #num_try+=1
            elif data[-1].endswith(" #no_cross-section\n"):
                c+=1
                # We increase this variable by 1, so if it is > 0 we do not redo the calculation 
                num_try+=1

            #remove the .slha
            if not os.path.isfile(slha_folder):
                slha_file_name = slha[6:-5]
            else:
                slha_file_name = name[6:-5]
            

            resummino_in_file = os.path.join(self.resummino_in,f"resummino_{slha_file_name}.in")
            resummino_out_file = os.path.join(self.resummino_out,f"resummino_{slha_file_name}.txt")
            resummino_log_file = os.path.join(self.resummino_log,f"resummino_log_{slha_file_name}.txt")
            
            #On prend le fichier de référence, et on créer une copie dans resummino_in avec le bon fichier slha
            self.modify_slha_file(self.input_file_original, resummino_in_file,slha_path)
            
            #On ajoute les noms à la liste (in, out et slha)
            Liste_resummino_in.append(resummino_in_file)
            Liste_output_file.append(resummino_out_file)
            Liste_slha.append(slha_path)

            # Here we list the channels to use, if scenario excluded then return None
            #particles = self.discrimination_particles(slha_path)
            if self.particles == None:
                self.mode, particles = self.extract_json()
            else:
                self.mode, particles = self.determine_channels()
                
            Liste_particles.append(particles)

            # We could optimize by removing variables that do not change from one iteration to another
            # But it is not very important (negligible in terms of compilation time
            # compared to Resummino)
            Liste.append((particles, resummino_in_file, slha_path, resummino_out_file, num_try, order, resummino_log_file))
        logger.info(f"{a-b-c} files created")
        return Liste


    def modify_slha_file(self, file_before, file_after, slha_file):
        """ _summary_

        Change all the informations in the .in files before launching calculations
        Args:
            file_before (input file for resummino): template
            file_after (input file for resummino): input file ready for resummino 
        """
        with open(file_before, 'r') as f:
            lines = f.readlines()

        try:        
            with open(self.json_resummino, "r") as fi:
                 data = eval(fi.read())
                        
            pdfs = data["pdfs"]
            
            pdf_lo = pdfs['pdf_lo']
            pdfset_lo = pdfs['pdfset_lo']
            pdf_nlo = pdfs['pdf_nlo']
            pdfset_nlo = pdfs['pdfset_nlo']
            
            lhapdf_folder = os.path.join(self.pwd, 'smodels', 'lib', 'resummino', 'lhapdf', 'share', 'LHAPDF')
            
            if not os.path.exists(os.path.join(lhapdf_folder, pdf_lo)):
                url = f"http://lhapdfsets.web.cern.ch/lhapdfsets/current/{pdf_lo}.tar.gz"
                try:
                    response = requests.get(url, stream=True)
                    response.raise_for_status() 

                    with open(f"{pdf_lo}.tar.gz", 'wb') as file:
                        for chunk in response.iter_content(chunk_size=8192): 
                            file.write(chunk)

                    with tarfile.open(f"{pdf_lo}.tar.gz", "r:gz") as tar:
                        tar.extractall(path=lhapdf_folder)

                except requests.exceptions.HTTPError as e:
                    logger.error("pdf name is wrong for LO, see http://lhapdfsets.web.cern.ch/lhapdfsets/current to check")
                    raise RuntimeError(f"Failed during HTTP request for {pdf_lo}, see http://lhapdfsets.web.cern.ch/lhapdfsets/current to check. You may also check if resummino is installed.")
                finally:
                    if os.path.exists(f"{pdf_lo}.tar.gz"):
                        os.remove(f"{pdf_lo}.tar.gz")
            if not os.path.exists(os.path.join(lhapdf_folder, pdf_nlo)):
                url = f"http://lhapdfsets.web.cern.ch/lhapdfsets/current/{pdf_nlo}.tar.gz"
                try:
                    response = requests.get(url, stream=True)
                    response.raise_for_status()  # Vérifie si la requête a réussi

                    # Ouvre un fichier temporaire pour enregistrer le fichier téléchargé
                    with open(f"{pdf_nlo}.tar.gz", 'wb') as file:
                        for chunk in response.iter_content(chunk_size=8192): 
                            file.write(chunk)

                    # Décompression du fichier
                    with tarfile.open(f"{pdf_nlo}.tar.gz", "r:gz") as tar:
                        tar.extractall(path=lhapdf_folder)

                except requests.exceptions.HTTPError as e:
                    logger.error("pdf name is wrong for NLO, see http://lhapdfsets.web.cern.ch/lhapdfsets/current to check")
                    raise RuntimeError(f"Failed during HTTP request for {pdf_nlo}, see http://lhapdfsets.web.cern.ch/lhapdfsets/current to check.")
                finally:
                    # Supprimer le fichier tar.gz temporaire si nécessaire
                    if os.path.exists(f"{pdf_nlo}.tar.gz"):
                        os.remove(f"{pdf_nlo}.tar.gz")
            with open(file_after, 'w') as f:
                for line in lines:
                    if line.startswith("pdf_lo"):
                        line = f"pdf_lo = {pdf_lo}\n"
                    if line.startswith("pdfset_lo"):
                        line = f"pdfset_lo = {pdfset_lo}\n"
                    if line.startswith("pdf_nlo"):
                        line = f"pdf_nlo = {pdf_nlo}\n"
                    if line.startswith("pdfset_nlo"):
                        line = f"pdfset_nlo = {pdfset_nlo}\n" 
                    f.write(line)
                        
        except KeyError:
                logger.info("choosing PDF4LHC21_40 by default")
           
        with open(file_after, 'w') as f:
            for line in lines:
                if line.startswith("slha ="):
                    line = f"slha = {slha_file}\n"
                if self.sqrt == '8.0':
                    if line.startswith("center_of_mass_energy ="):
                        line = f"center_of_mass_energy = 8000"
                if self.sqrt == '13.6':
                    if line.startswith("center_of_mass_energy ="):
                        line = f"center_of_mass_energy = 13600"
                f.write(line)


    def extract_json(self):
        """_summary_
        
        function to extract all the informations in the resummino.py
        file

        Returns: tuple of:
            string: Mode of writting for the slha cross section
            list: list of the daugther particle to consider in the calculation
            of the cross section
        """
        with open(self.json_resummino, "r") as f:
            data = eval(f.read())
        
        try:
            mode = data["mode"]
        except KeyError:
            mode = "check"

        
        channels = data["channels"]
        particles = []
        for key, values in channels.items():
            particles.append((values[0], values[1]))
            
        return mode, particles
    
    def determine_channels(self):
        """_summary_
        
        function to find channels using a set of particles

        Returns: tuple of:
            string: Mode of writting for the slha cross section
            list: list of the daugther particle to consider in the calculation
            of the cross section
        """
        if 1000024 in self.particles:
            self.particles.append(-1000024)
        if 1000037 in self.particles:
            self.particles.append(-1000037)
        channels = list(combinations(self.particles, 2))
        
        mode = 'check'
            
        return mode, channels

    def modify_outgoing_particles( self, input_file, output_file, new_particle1, 
                                   new_particle2 ):
        """
        modify the output particles (mother particles) in the resummino .in
        file. First call the template (input_file),
        then write into the resummino .in file (output_file). 
        Can also write directly onto the output_file.
        """
        with open(input_file, 'r') as f:
            lines = f.readlines()

        with open(output_file, 'w') as f:
            for line in lines:
                if line.startswith("particle1 ="):
                    line = f"particle1 = {new_particle1}\n"
                elif line.startswith("particle2 ="):
                    line = f"particle2 = {new_particle2}\n"
                f.write(line)
            
            
    def launch_all(self):
        """
        Launch all the calculations of the slha files in parallel (limited by
        ncpu), with first the creation of every path needed for the calculations.
        """
        # We create the list
        tasks = self.create_routine_files(self.maxOrder, self.slha_folder_name)
        
        if not os.path.exists(os.path.join(self.pwd, 'smodels', 'lib', 'resummino', 'lhapdf', 'bin', 'lhapdf')):
            logger.error("Warning, lhapdf is not installed, please make resummino, or use ./install.sh in the smodels/lib/resummino folder.")
            return
            
        if not os.path.exists(os.path.join(self.pwd, 'smodels', 'lib', 'resummino', 'resummino_install', 'bin', 'resummino')):
            logger.error("Warning, resummino is not installed, please make resummino, or use ./install.Sh in the smodels/lib/resummino folder.")
            return
        try:
            if len(tasks)==0:
                logger.error("Warning, check your inputfile (.slha), no calculation are available.")
                return
        except TypeError:
            logger.error("Warning, check your inputfile (.slha), no calculation are available.")
            return
        
        if self.ncpu == 1:
            for task in tasks:
                self.calculate_one_slha(*task)
        else:   
            # We launch the program with maximum performance, change if necessary
            with ProcessPoolExecutor(max_workers=self.ncpu) as executor:
                futures = [executor.submit(self.calculate_one_slha, *task) for task in tasks]
                for future in futures:
                    future.result()

        shutil.rmtree(self.resummino_in)
        shutil.rmtree(self.resummino_out)
        shutil.rmtree(self.resummino_log)

def main ( args : argparse.Namespace ):
    """ the central entry point """
    canonizer = ArgsStandardizer()
    setLogLevel ( args.verbosity )
    if not hasattr ( args, "noautocompile" ):
        args.noautocompile = False
        
    sqrtses = canonizer.getSqrtses ( args )
    order = canonizer.getOrder ( args )
    canonizer.checkAllowedSqrtses ( order, sqrtses )
    inputFiles = args.filename.strip()
    ncpus = canonizer.checkNCPUs ( args.ncpus, inputFiles )
    type_writting = canonizer.writeToFile(args)
    json = canonizer.getjson(args)
    particles = canonizer.getParticles( args )
    xsec_limit = canonizer.checkXsec_limit(args)
    #If we get no type_writing then we only print the result 
    if type_writting == None :
        type_writting = None
        
    verbosity = args.verbosity

    ssmultipliers = None
    if hasattr ( args, "ssmultipliers" ):
        ssmultipliers = canonizer.getSSMultipliers ( args.ssmultipliers )
        if ssmultipliers != None:
            for pids,multiplier in ssmultipliers.items():
                if type(pids) != tuple:
                    logger.error ( "keys of ssmultipliers need to be supplied as tuples" )
                    sys.exit()
                if type(multiplier) not in [ int, float ]:
                    logger.error ( "values of ssmultipliers need to be supplied as ints or floats" )
                    sys.exit()

    logger.debug('verbosity is '+ verbosity)
    ssqrtses = ",".join(map(str,sqrtses))
    
    logger.info(f"the calculation will be performed using {ssqrtses} TeV as center of mass energy")
    # logger.debug("The calculation will be performed using : " +str(sqrtses)+ ' TeV as center of mass energy')
    orders_dic = {0:'LO', 1:'NLO', 2:'NLL+NLO'}

    logger.info("the highest perturbation order considered is " + orders_dic[order])
    scpu = "cpus"
    if ncpus == 1:
        scpu = "cpu"
    logger.info(f"we are currently running on {ncpus} {scpu}" )
    logger.debug(f"In this calculation, we'll write out the cross section at "+ str(type_writting) +" order of perturbation theory below "+ orders_dic[order])
    
    for sqrt in sqrtses:
        logger.debug('Current energy considered is '+ str(sqrt)+ ' TeV')
        test = XSecResummino(maxOrder=order, slha_folder_name=inputFiles, sqrt = sqrt, ncpu=ncpus, type_writting = type_writting, verbosity = verbosity, json = json, particles=particles, xsec_limit = xsec_limit)
        test.checkInstallation ( compile = not args.noautocompile )
        test.launch_all()
    return
    
    # test = XSecResummino(maxOrder=order, slha_folder_name=inputFiles, sqrt = sqrtses, ncpu=ncpus, type = type_writting)
    # test.routine_resummino()
 
    # test = XSecResummino(maxOrder=order, slha_folder_name=inputFiles, sqrt = sqrtses, ncpu=ncpus, type = type_writting)
    # test.routine_resummino()
