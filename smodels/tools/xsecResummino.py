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
from smodels.tools import toolBox, runtime
from smodels.tools.physicsUnits import pb, TeV, GeV
from smodels.theory import crossSection
from smodels.theory.crossSection import LO, NLO, NLL
from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import subprocess
from concurrent.futures import ProcessPoolExecutor
from smodels.tools.xsecBasis import XSecBasis, ArgsStandardizer
import pyslha
import math
import json
try:
    import cStringIO as io
except ImportError as e:
    import io


class XSecResummino(XSecBasis):
    """ cross section computer class (for resummino), what else? """
    def __init__ ( self, maxOrder,slha_folder_name,sqrt = 13,ncpu=1, maycompile=True, type = 'all', verbosity = '', json = None):
        """
        :param maxOrder: maximum order to compute the cross section, given as an integer
                    if maxOrder == LO, compute only LO resummino xsecs
                    if maxOrder == NLO, compute NLO resummino xsecs
                    if maxOrder == NLL, compute NLO+NLL resummino xsecs
        :param nevents: number of events for pythia run
        :param pythiaVersion: pythia6 or pythia8 (integer)
        :param sqrt: Center of mass energy to consider for the cross section calculation
        :param mode_limit: Value below which if mode == "check", cross section at NLO order are not calculated
        :param type: If all, put all the order in the slha file, if highest, only the hightest order.
        :param json: Path to the json file with all the relevant informations concerning the resummino calculation
        :param resummino_bin: Path to resummino executable
        :param input_file_original: Path to the template input file of resummino
        :param ncpu: Number of cpu used in parallel for the calculation
        :param verbosity: Type of informations given in the logger file
        """
        self.pwd = installation.installDirectory()
        self.resummino_bin = os.path.join(self.pwd,"./smodels/lib/resummino/resummino-3.1.2/bin/resummino")
        self.input_file_original = os.path.join(self.pwd,"smodels/etc/input_resummino.in")
        if json == None:
            self.json_resummino = os.path.join(self.pwd,"smodels/etc/resummino.py")
        else:
            if os.path.isabs(json):
                self.json_resummino = json
            else:
                self.json_resummino = os.path.join(self.pwd, json)

        
        self.slha_folder_name = slha_folder_name
        self.maxOrder = maxOrder
        self.countNoXSecs = 0
        self.countNoNLOXSecs = 0
        self.maycompile = maycompile
        self.ncpu = ncpu
        self.sqrts = sqrt
        self.type = type
        self.verbosity = verbosity
        self.sqrt = sqrt
        
        
        if verbosity == 'debug':
            logger.info('installation directory is '+self.pwd)
            
        with open(self.json_resummino, "r") as f:
            data = eval(f.read())
                             
        self.mode_limit = data["mode_limit"]
        
    def calculate_one_slha(self, particles,input_file, slha_file, output_file, num_try, order, log):
        """
        log file management and launch resummino command
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
            
        nxsecs = self.addXSecToFile(Xsections, slha_file, comment = "[pb], Resumminov3.1.2")
        
           
    def launch_command(self,resummino_bin,input_file, output_file, order):
        """
        use resummino at the order asked by the user (order variable)
        """
        if order == 0:
            commande = f"{resummino_bin} {input_file} --lo"
        if order == 1:
            commande = f"{resummino_bin} {input_file} --nlo"
        if order == 2:
            commande = f"{resummino_bin} {input_file}"

        with open(output_file, 'w') as f:
            subprocess.run(commande, shell=True, stdout=f, text=True)



    def launch_resummino(self, input_file, slha_file, output_file, particle_1, particle_2, num_try, order, Xsections, log):
        """
        Check everything before launching resummino
        """
        
        #modifie_slha_file(input_file, input_file, slha_file)
        self.modifie_outgoing_particles(input_file, input_file, particle_1, particle_2)
        
        already_written_canal = self.find_channels(slha_file)

        if self.sqrt > 10:
            _ = str(self.sqrt*10**(-1))+'0E+04'
        else:
            _ = str(self.sqrt)+'.00E+03'
        
        already_written_canal_set = [({x,y},z,w) for (x,y), z,w in already_written_canal]
        
        if self.verbosity == 'debug':
            logger.info("channel, order and cross section" + str(particle_1)+str(particle_2)+ str(order)+ str(_))
            logger.info('the already writtent channels are'+ str(already_written_canal))
        if (((particle_1, particle_2), _, order)) in already_written_canal:
            return
        
        if self.mode == "check":
            self.launch_command(self.resummino_bin, input_file, output_file, 0)
            infos = self.search_in_output(output_file)
            infos = infos[0].split(" ")[2][1:]
            if self.verbosity == 'info' or self.verbosity == 'debug':
                logger.info("cross section is "+str(infos)+ " pb at LO order")
                logger.info("Is cross section above the limit ? "+str( float(infos)>10**(-5)))
            if (float(infos))>(10**(-5)):
                if self.verbosity == 'debug':
                    logger.info('num try is '+ str(num_try)+ ' for the moment')
                self.launch_command(self.resummino_bin, input_file, output_file, order)
                if num_try == 0:
                    hist = self.write_in_slha(output_file, slha_file, order, particle_1, particle_2, self.type, Xsections, log)
            else:
                hist = self.write_in_slha(output_file, slha_file, 0, particle_1, particle_2, self.type, Xsections, log)
            return
        
        if self.mode == "all":
            self.launch_command(self.resummino_bin, input_file, output_file, order)
            if num_try == 0:
                hist = self.write_in_slha(output_file, slha_file, order, particle_1, particle_2, self.type, Xsections, log)
            return
        
        #:hist: variable to check if there is a problem in the cross section (strange value or no value)
        hist = 0

        
        if num_try == 0 and self.mode != "check":
            self.launch_command(self.resummino_bin, input_file, output_file, order)


        #here we write in the slha file.
            hist = self.write_in_slha(output_file, slha_file, order, particle_1, particle_2, self.type, Xsections, log)

        #we check if we have written too much cross section
        self.are_crosssection(slha_file, order)

        #if there is an error, we inform the user and launch again resummino
        if hist == 1:
            if self.verbosity == 'warning':
                self.info("error in resummino, trying again")
            num_try = 0
            self.modifie_outgoing_particles(input_file, input_file, particle_1, particle_2)
            self.launch_resummino(input_file, slha_file, output_file, particle_1, particle_2, num_try, order, Xsections, log)


    def search_in_output(self, output_file):
        Infos = []
        with open(output_file, 'r') as f:
            data = f.readlines()
        for i in range(len(data)):
            if "Results:" in data[i]:
                LO = data[i+1][:-1] #[:-1] to get rid of the \n
                NLO = data[i+2][:-1]
                NLL = data[i+3][:-1]
                Infos.append((LO,NLO,NLL))
        return Infos[0]

    # def writing_result(self, value, particle_1, particle_2, slha_file, order, type_writing):
    #     with open(slha_file, 'a') as f:
    #         if type_writing == 'all':
    #             f.write(f"\nXSECTION  1.30E+04  2212 2212 2 {particle_1} {particle_2} # [pb] \n")
    #             for i in range(order+1):
    #                 f.write(f" 0  0  {i}  0  0  0    {value[i]} Resumminov3.1.2\n")
    #         if type_writing == 'max':
    #             f.write(f"\nXSECTION  1.30E+04  2212 2212 2 {particle_1} {particle_2} # [pb] \n 0  0  {order}  0  0  0    {value} Resumminov3.1.2\n")
            
    def create_xsection(self, result, particle_1, particle_2, order, Xsections):
        
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
        results = self.search_in_output(output_file)
        if type_writing == 'highest':
            if order == 0:
                result = results[0].split(" ")[2][1:]
            elif order == 1:
                result = results[1].split(" ")[2][1:]
            elif order == 2:
                result = results[2].split(" ")[2][1:]
            logger.info('the highest result is'+ str(result) +' pb')
        elif type_writing == "all":
            result = [results[0].split(" ")[2][1:], results[1].split(" ")[2][1:], results[2].split(" ")[2][1:]]
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
            
        self.create_xsection(result, particle_1, particle_2, order, Xsections)
            
        return 0

    def extract_m1_m2_mu(self, file_path):
        """_summary_
        
        function to extract the breaking term of the electrowikino part (SUSY) in
        an slha file.
        Args:
            file_path (_string_): _path of the slha file_

        Returns:
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
                channel = (int(line.split(" ")[7]), int(line.split(" ")[8]),line.split(" ")[2], data[i+1].split(" ")[4])
                
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
                sqrt = line.split(" ")[2]
                channel = (int(line.split(" ")[7]), int(line.split(" ")[8]))
                order = data[i+1].split(" ")[4]
                order = int(order)
                channels.append((channel, sqrt, order))
        return channels


    def discrimination_particles(self, slha_file):
        """
        Choice of the differents channel and condtions on the cross section calculation.
        C1<92 or N1>600 or C2>1200 by default.
        """
        m1,m2,mu = self.extract_m1_m2_mu(slha_file)
        N1,N2,C1,C2 = self.extract_N1_N2_C1(slha_file)
        abs_mu = math.fabs(mu)
        if C1<92 or N1>600 or C2>1200:
            return None
        if m2>=abs_mu:
            return [(1000025,1000037),(1000025,-1000037), (-1000037, 1000037)]
        elif m2<abs_mu:
            return [(1000023,1000037), (1000023,-1000037), (1000025, 1000037), (1000025, -1000037), (-1000037,1000037), (1000023,1000025)]
    
      
    def routine_creation(self, order, slha_folder_name):
    #slha_file = "outputM_12000M_20mu100.slha"

        #You haven't seen anything.
        output_file = "output_file.txt"

        #List of the differents files on which resummino will be launch
        Liste_slha = []
        Liste_resummino_in = []
        Liste_output_file = []
        Liste_particles = []
        Liste = []

        resummino_in = os.path.join(self.pwd, 'smodels/lib/resummino/resummino_in')
        resummino_out = os.path.join(self.pwd, 'smodels/lib/resummino/resummino_out')
        resummino_log = os.path.join(self.pwd, 'smodels/lib/resummino/resummino_log')
        #Check if the folder already exists
        if not os.path.exists(resummino_in):
            os.mkdir(resummino_in)
        if not os.path.exists(resummino_out):
            os.mkdir(resummino_out)
        if not os.path.exists(resummino_log):
            os.mkdir(resummino_log)
        
        #always use absolute path
        slha_folder = os.path.join(self.pwd, slha_folder_name)  
        #Check if the input is a file or a directory
        if not os.path.isfile(slha_folder):
            liste_slha = os.listdir(slha_folder_name)
        else:
            liste_slha = [slha_folder_name]
            
    

        

        '''
        a = number of files
        b = number of files with already a cross section
        c = number of files ignored
        '''
        a,b,c = 0,0,0

        #Loop on all the slha file, in a random order (os.listdir)
        for slha in liste_slha:
            
            print(slha, slha_folder)
            if not os.path.isfile(slha_folder):
                slha_path = os.path.join(slha_folder,slha)
            else:
                slha_path = os.path.join(self.pwd, slha)
            #Variable used to avoid strange result (huge difference between LO and NLO result)
            num_try = 0
            with open(slha_path, 'r') as f:
                data = f.readlines()
                a+=1
                
            if data[-1].endswith("SModelSv2.3.0\n"):
                b+=1
                #On augmente cette variable de 1, comme ca si elle est > 0 on ne refait pas le calcul
                #num_try+=1
            elif data[-1].endswith(" #no_cross-section\n"):
                c+=1
                #On augmente cette variable de 1, comme ca si elle est > 0 on ne refait pas le calcul
                num_try+=1

            #remove the .slha
            slha_file_name = slha[6:-5]

            

            resummino_in_file = os.path.join(resummino_in,f"resummino_{slha_file_name}.in")
            resummino_out_file = os.path.join(resummino_out,f"resummino_{slha_file_name}.txt")
            resummino_log_file = os.path.join(resummino_log,f"resummino_log_{slha_file_name}.txt")
            
            #On prend le fichier de référence, et on créer une copie dans resummino_in avec le bon fichier slha
            self.modifie_slha_file(self.input_file_original, resummino_in_file,slha_path)
            
            #On ajoute les noms à la liste (in, out et slha)
            Liste_resummino_in.append(resummino_in_file)
            Liste_output_file.append(resummino_out_file)
            Liste_slha.append(slha_path)

            #On liste ici les canaux à utiliser, si scénario exclu alors renvoi None
            #particles = self.discrimination_particles(slha_path)
            self.mode, particles = self.extract_json()
            Liste_particles.append(particles)

            #On pourrait optimiser en enlevant les variables qui ne changent pas d'une itération à l'autre
            #Mais ce n'est pas très important (négligeable niveau temps de compilation comparé à Resummino)
            Liste.append((particles, resummino_in_file, slha_path, resummino_out_file, num_try, order, resummino_log_file))
        print(f" Number of files created : {a-b-c}")
        return Liste


    def modifie_slha_file(self, file_before, file_after, slha_file):
        """_summary_

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
                comand = f"wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/{pdf_lo}.tar.gz -O- | tar xz -C {lhapdf_folder}"
                subprocess.run(comand, shell=True)
                
            if not os.path.exists(os.path.join(lhapdf_folder, pdf_nlo)):
                comand = f"wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/{pdf_nlo}.tar.gz -O- | tar xz -C {lhapdf_folder}"
                subprocess.run(comand, shell=True)
                
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
        
        function to extract all the informations in the resummino.json
        file

        Returns:
            string: Mode of writting for the slha cross section
            list: liste of the daugther particle to consider in the calculation of the cross section
        """
        with open(self.json_resummino, "r") as f:
            data = eval(f.read())
        
        mode = data["mode"]
        
        channels = data["channels"]
        particles = []
        for clef, valeurs in channels.items():
            particles.append((valeurs[0], valeurs[1]))
            
        return mode, particles

    def modifie_outgoing_particles(self, input_file, output_file, new_particle1, new_particle2):
        with open(input_file, 'r') as f:
            lines = f.readlines()

        with open(output_file, 'w') as f:
            for line in lines:
                if line.startswith("particle1 ="):
                    line = f"particle1 = {new_particle1}\n"
                elif line.startswith("particle2 ="):
                    line = f"particle2 = {new_particle2}\n"
                f.write(line)
            
            
    def routine_resummino(self):
        """
        Launch all the calculations of the slha files in parallel (limited by ncpu), with first the creation of every path needed for the calculations.
        """
        #On créer la liste
        tasks = self.routine_creation(self.maxOrder, self.slha_folder_name)
        #On lance le programme avec les performances maximales, à changer si besoin
        with ProcessPoolExecutor(max_workers=self.ncpu) as executor:
            futures = [executor.submit(self.calculate_one_slha, *task) for task in tasks]
            for future in futures:
                future.result()



def main(args):
    
    
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
    #We choose to select highest by default
    if type_writting == None :
        type_writting = 'highest'
        
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

    logger.info('verbosity is '+ verbosity)
    if verbosity == 'info' or verbosity == 'debug':
        logger.info("The calculation will be done using :" +str(sqrtses)+ ' TeV as center of mass energy')
        logger.info("The max order considered for the calculation is  " + str(order)+ ' (0 = LO, 1 = NLO, 2 = NLL+NLO)')
        logger.info("we are currently running on " + str(ncpus)+ ' cpu')
        logger.info(f"In this calculation, we'll use "+ str(type_writting) +" type of writting for the cross-section")
    
    for sqrt in sqrtses:
        if verbosity == 'info':
            print('Current energy considered is '+ str(sqrt)+ ' TeV')
        test = XSecResummino(maxOrder=order, slha_folder_name=inputFiles, sqrt = sqrt, ncpu=ncpus, type = type_writting, verbosity = verbosity, json = json)
        test.routine_resummino()
    return
    
    # test = XSecResummino(maxOrder=order, slha_folder_name=inputFiles, sqrt = sqrtses, ncpu=ncpus, type = type_writting)
    # test.routine_resummino()
    

