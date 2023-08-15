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
    """ cross section computer class, what else? """
    def __init__ ( self, maxOrder,slha_folder_name,sqrt = 13,ncpu=1, maycompile=True, type = 'all'):
        """
        :param maxOrder: maximum order to compute the cross section, given as an integer
                    if maxOrder == LO, compute only LO pythia xsecs
                    if maxOrder == NLO, apply NLO K-factors from NLLfast (if available)
                    if maxOrder == NLL, apply NLO+NLL K-factors from NLLfast (if available)
        :param nevents: number of events for pythia run
        :param pythiaVersion: pythia6 or pythia8 (integer)
        :param maycompile: if True, then tools can get compiled on-the-fly
        """
        self.resummino_bin = "./smodels/lib/resummino/resummino-3.1.2/bin/resummino"
        self.input_file_original = "smodels/etc/ff1a240db6c1719fe9f299b3390d49d32050c4f1003286d2428411eca45bd50c.in"
        self.json_resummino = "smodels/etc/resummino.json"
        self.slha_folder_name = slha_folder_name
        self.maxOrder = maxOrder
        self.countNoXSecs = 0
        self.countNoNLOXSecs = 0
        self.maycompile = maycompile
        self.ncpu = ncpu
        self.sqrts = sqrt
        self.type = type
        # sqrts = [float(x) for x in self.sqrts]
        # self.sqrt = sqrts[0]
        self.sqrt = sqrt
    
        with open(self.json_resummino, "r") as f:
            data = json.load(f)
                             
        self.mode_limit = data["mode_limit"]
        
    def one_slha_calculation(self, particles,input_file, slha_file, output_file, num_try, order, log):
        """
        log file management and launch resummino command
        """
        with open(log, 'a') as f:
            f.write(f'{particles} cross-sections written in {slha_file}\n')
        if particles == None:
            with open(slha_file, 'a') as f:
                f.write(' #no_cross-section\n')

        #Utilisé ici pour vérifier si on a pas écrit 2 fois le #no_cross-section    
        self.are_crosssection(slha_file, order)
        if particles == None:
            return
        Xsections = crossSection.XSectionList()
        
        for particle_pair in particles:
            self.launcher(input_file, slha_file, output_file, particle_pair[0], particle_pair[1], num_try, order, Xsections, log)
            
        nxsecs = self.addXSecToFile(Xsections, slha_file)
        
           
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



    def launcher(self, input_file, slha_file, output_file, particle_1, particle_2, num_try, order, Xsections, log):
        """
        Check everything before launching resummino
        """
        
        #modifie_slha_file(input_file, input_file, slha_file)
        self.modifie_outgoing_particles(input_file, input_file, particle_1, particle_2)
        
        already_written_canal = self.canaux_finding(slha_file)

        _ = str(self.sqrt*10**(-1))+'0E+04'
        print("channel, order and cross section",particle_1,particle_2, order, _)
        print(already_written_canal)
        if (((particle_1, particle_2), _, order)) in already_written_canal:
            print('wow')
            return
        
        if self.mode == "check":
            self.launch_command(self.resummino_bin, input_file, output_file, 0)
            infos = self.search_in_output(output_file)
            infos = infos[0].split(" ")[2][1:]
            print("cross section is ", infos, "at LO order")
            print("Is cross section above the limit ?", float(infos)>10**(-5))
            if (float(infos))>(10**(-5)):
                print(num_try)
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
            print("error")
            num_try = 0
            self.modifie_outgoing_particles(input_file, input_file, particle_1, particle_2)
            self.launcher(input_file, slha_file, output_file, particle_1, particle_2, num_try, order, Xsections, log)
    #launcher(input_file, slha_file, output_file, particle_1, particle_2)


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
        Xsection.info.sqrts = 13. * TeV
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
            print('the highest result is', result)
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
            #_ = math.fabs(results[0].split(" ")[2][1:]-results[1].split(" ")[2][1:])
            if float(results[1].split(" ")[2][1:])>2*float(results[0].split(" ")[2][1:]) or float(results[0].split(" ")[2][1:])> 2* float(results[1].split(" ")[2][1:]):
                with open(log, 'a') as f:
                    f.write(f"to much change between LO and NLO for {particle_1} and {particle_2} with {slha_file}")
                return 1
        if order == 2:
            #_ = math.fabs(results[0].split(" ")[2][1:]-results[1].split(" ")[2][1:])
            if float(results[2].split(" ")[2][1:])>2*float(results[0].split(" ")[2][1:]) or float(results[0].split(" ")[2][1:])> 2* float(results[1].split(" ")[2][1:]):
                with open(log, 'a') as f:
                    f.write(f"to much change between LO and NLL+NLO for {particle_1} and {particle_2} with {slha_file}")
                return 1
            
        self.create_xsection(result, particle_1, particle_2, order, Xsections)
            
        #self.writing_result(result, particle_1, particle_2, slha_file, order, type_writing)
        return 0

    def extract_m1_m2_mu(self, file_path):
        """_summary_
        
        function to extract the breaking term of the electrowikino part (SUSY) in
        an slha file.
        Args:
            file_path (_string_): _path of the slha file_

        Returns:
            _type_: _description_
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
                    print(f'remove from {slha_file} a #no_cross-section"')
                else:
                    test = False
            with open(slha_file, 'w') as f:
                for _ in data:
                    f.write(_)
                return
        canaux = {}
        to_delete = []

        for i in range(len(data)):
            line = data[i]
            if line.startswith("XSECTION"):
                canal = (int(line.split(" ")[7]), int(line.split(" ")[8]),line.split(" ")[2], data[i+1].split(" ")[4])
                
                if canal in canaux:
                    start = canaux[canal]
                    end = start+1
                    while not data[end].startswith('XSECTION'):
                        end+=1
                    #end = start+order+3
                    to_delete.extend(range(start, end))
                canaux[canal] = i
        lines = [line for i, line in enumerate(data) if i not in to_delete]

        with open(slha_file, 'w') as f:
            f.writelines(lines)
        
        
    def canaux_finding(self, slha_file):
        
        with open(slha_file, 'r') as f:
            data = f.readlines()
            
        canaux = []
        for i in range(len(data)):
            line = data[i]
            if line.startswith("XSECTION"):
                sqrt = line.split(" ")[2]
                canal = (int(line.split(" ")[7]), int(line.split(" ")[8]))
                order = data[i+1].split(" ")[4]
                order = int(order)
                canaux.append((canal, sqrt, order))
        return canaux


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

        #Check if the folder already exists
        if not os.path.exists('smodels/lib/resummino/resummino_in'):
            os.mkdir('smodels/lib/resummino/resummino_in')
        if not os.path.exists('smodels/lib/resummino/resummino_out'):
            os.mkdir('smodels/lib/resummino/resummino_out')
        if not os.path.exists('smodels/lib/resummino/resummino_log'):
            os.mkdir('smodels/lib/resummino/resummino_log')
            
        #Check if the input is a file or a directory
        if not slha_folder_name.endswith(".slha"):
            liste_slha = os.listdir(slha_folder_name)
        else:
            liste_slha = [slha_folder_name]
            
        #current directory
        pwd = os.getcwd()

        #always use absolute path
        slha_folder = os.path.join(pwd, slha_folder_name)

        '''
        a = number of files
        b = number of files with already a cross section
        c = number of files ignored
        '''
        a,b,c = 0,0,0

        #Loop on all the slha file, in a random order (os.listdir)
        for slha in liste_slha:

            slha_path = os.path.join(slha_folder,slha)
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

            #On prend le fichier de référence, et on créer une copie dans resummino_in avec le bon fichier slha
            self.modifie_slha_file(self.input_file_original, f"smodels/lib/resummino/resummino_in/resummino_{slha_file_name}.in",slha_path)

            #On ajoute les noms à la liste (in, out et slha)
            Liste_resummino_in.append(f"smodels/lib/resummino/resummino_in/resummino_{slha_file_name}.in")
            Liste_output_file.append(f"smodels/lib/resummino/resummino_out/resummino_out_{slha_file_name}.txt")
            Liste_slha.append(slha_path)

            #On liste ici les canaux à utiliser, si scénario exclu alors renvoi None
            #particles = self.discrimination_particles(slha_path)
            self.mode, particles = self.json_extraction()
            Liste_particles.append(particles)

            #On pourrait optimiser en enlevant les variables qui ne changent pas d'une itération à l'autre
            #Mais ce n'est pas très important (négligeable niveau temps de compilation comparé à Resummino)
            Liste.append((particles, f"smodels/lib/resummino/resummino_in/resummino_{slha_file_name}.in", slha_path, f"smodels/lib/resummino/resummino_out/resummino_out_{slha_file_name}.txt", num_try, order, f'smodels/lib/resummino/resummino_log/resummino_log_{slha_file_name}.txt'))
        print(f" Number of files created : {a-b-c}")
        return Liste


    def modifie_slha_file(self, file_before, file_after, slha_file):
        """_summary_

        Change all the informations in the .in files before launching calculations
        Args:
            file_before (input file for resummino): template
            file_after (input file for resummino): input file ready for resu;;ino 
        """
        with open(file_before, 'r') as f:
            lines = f.readlines()

        try:        
            with open(self.json_resummino, "r") as fi:
                 data = json.load(fi)
                        
            pdfs = data["pdfs"]
            
            pdf_lo = pdfs['pdf_lo']
            pdfset_lo = pdfs['pdfset_lo']
            pdf_nlo = pdfs['pdf_nlo']
            pdfset_nlo = pdfs['pdfset_nlo']
            
            pwd = os.getcwd()
            lhapdf_folder = os.path.join(pwd, 'smodels', 'lib', 'resummino', 'lhapdf', 'share', 'LHAPDF')
            
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
                print("choosing PDF4LHC21_40 by default")
           
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


    def json_extraction(self):
        """_summary_
        
        function to extract all the informations in the resummino.json
        file

        Returns:
            string: Mode of writting for the slha cross section
            list: liste of the daugther particle to consider in the calculation of the cross section
        """
        with open(self.json_resummino, "r") as f:
            data = json.load(f)
        
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
    
        #On créer la liste
        tasks = self.routine_creation(self.maxOrder, self.slha_folder_name)
        #On lance le programme avec les performances maximales, à changer si besoin
        with ProcessPoolExecutor(max_workers=self.ncpu) as executor:
            futures = [executor.submit(self.one_slha_calculation, *task) for task in tasks]
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

    print("sqrt :" +str(sqrtses))
    print("order" + str(order))
    print("ncpu" + str(ncpus))
    print([float(x) for x in sqrtses])
    
    for sqrt in sqrtses:
        test = XSecResummino(maxOrder=order, slha_folder_name=inputFiles, sqrt = sqrt, ncpu=ncpus, type = type_writting)
        test.routine_resummino()
    return
    
    # test = XSecResummino(maxOrder=order, slha_folder_name=inputFiles, sqrt = sqrtses, ncpu=ncpus, type = type_writting)
    # test.routine_resummino()
    

