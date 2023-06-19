#!/usr/bin/env python3

"""
.. module:: xsecComputer
   :synopsis: Computation of reference ("theory") production cross sections.

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
from __future__ import print_function
from smodels import installation
from smodels.tools import toolBox, runtime
from smodels.tools.physicsUnits import pb, TeV, GeV
from smodels.theory import crossSection
from smodels.theory.crossSection import LO, NLO, NLL
from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import subprocess
from concurrent.futures import ProcessPoolExecutor
import os, copy
import pyslha
import math
try:
    import cStringIO as io
except ImportError as e:
    import io
import sys

class XSecResummino:
    """ cross section computer class, what else? """
    def __init__ ( self, maxOrder,slha_folder_name, maycompile=True):
        """
        :param maxOrder: maximum order to compute the cross section, given as an integer
                    if maxOrder == LO, compute only LO pythia xsecs
                    if maxOrder == NLO, apply NLO K-factors from NLLfast (if available)
                    if maxOrder == NLL, apply NLO+NLL K-factors from NLLfast (if available)
        :param nevents: number of events for pythia run
        :param pythiaVersion: pythia6 or pythia8 (integer)
        :param maycompile: if True, then tools can get compiled on-the-fly
        """
        self.resummino_bin = "./smodels/tools/resummino-3.1.2/bin/resummino"
        self.input_file_original = "smodels/tools/ff1a240db6c1719fe9f299b3390d49d32050c4f1003286d2428411eca45bd50c.in"
        self.slha_folder_name = slha_folder_name
        self.maxOrder = maxOrder
        self.countNoXSecs = 0
        self.countNoNLOXSecs = 0
        self.maycompile = maycompile

    def addXSecToFile( self, xsecs, slhafile, comment=None, complain=True):
        """
        Write cross sections to an SLHA file.

        :param xsecs: a XSectionList object containing the cross sections
        :param slhafile: target file for writing the cross sections in SLHA format
        :param comment: optional comment to be added to each cross section block
        :param complain: complain if there are already cross sections in file

        """

        if not os.path.isfile(slhafile):
            line = f"SLHA file {slhafile} not found."
            logger.error( line )
            raise SModelSError( line )
        if len(xsecs) == 0:
            self.countNoXSecs+=1
            if self.countNoXSecs < 3:
                logger.warning("No cross sections available.")
            if self.countNoXSecs == 3:
                logger.warning("No cross sections available (will quench such warnings in future).")
            return False
        # Check if file already contain cross section blocks
        xSectionList = crossSection.getXsecFromSLHAFile(slhafile)
        if xSectionList and complain:
            logger.info("SLHA file already contains XSECTION blocks. Adding "
                           "only missing cross sections.")

        # Write cross sections to file, if they do not overlap with any cross
        # section in the file
        outfile = open(slhafile, 'r')
        lastline = outfile.readlines()[-1]
        lastline = lastline.strip()
        outfile.close()

        outfile = open(slhafile, 'a')
        if lastline != "":
            outfile.write("\n" )
        nxsecs = 0
        for xsec in xsecs:
            writeXsec = True
            for oldxsec in xSectionList:
                if oldxsec.info == xsec.info and set(oldxsec.pid) == set(xsec.pid):
                    writeXsec = False
                    break
            if writeXsec:
                nxsecs += 1
                outfile.write( self.xsecToBlock(xsec, (2212, 2212), comment) + "\n")
        outfile.close()

        return nxsecs 

    def xsecToBlock( self, xsec, inPDGs=(2212, 2212), comment=None, xsecUnit = pb):
        """
        Generate a string for a XSECTION block in the SLHA format from a XSection
        object.

        :param inPDGs: defines the PDGs of the incoming states
                    (default = 2212,2212)

        :param comment: is added at the end of the header as a comment
        :param xsecUnit: unit of cross sections to be written (default is pb).
                        Must be a Unum unit.

        """
        if type(xsec) != type(crossSection.XSection()):
            logger.error("Wrong input")
            raise SModelSError()
        # Sqrt(s) in GeV
        header = "XSECTION  " + str(xsec.info.sqrts / GeV)
        for pdg in inPDGs:
            # PDGs of incoming states
            header += " " + str(pdg)
        # Number of outgoing states
        header += " " + str(len(xsec.pid))
        for pid in xsec.pid:
            # PDGs of outgoing states
            header += " " + str(pid)
        if comment:
            header += " # " + str(comment)  # Comment
        entry = "  0  " + str(xsec.info.order) + "  0  0  0  0  " + \
                str( "%16.8E" % (xsec.value / xsecUnit) ) + " SModelSv" + \
                    installation.version()

        return "\n" + header + "\n" + entry
    
    def one_slha_calculation(self, particles,input_file, slha_file, output_file, num_try, order):
        """
        Gestion du fichier log et lancement de la commande lancant Resummino
        """
        with open('log.txt', 'a') as f:
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
            self.launcher(input_file, slha_file, output_file, particle_pair[0], particle_pair[1], num_try, order, Xsections)
            
        nxsecs = self.addXSecToFile(Xsections, slha_file)
        
           
    def launch_command(self,resummino_bin,input_file, output_file, order):
        """
        Lance la commande pour resummino, pour modifier l'emplacement de resummino
        il faut préciser le resummino_bin.
        """
        if order == 0:
            commande = f"{resummino_bin} {input_file} --lo"
        if order == 1:
            commande = f"{resummino_bin} {input_file} --nlo"
        if order == 2:
            commande = f"{resummino_bin} {input_file}"

        with open(output_file, 'w') as f:
            subprocess.run(commande, shell=True, stdout=f, text=True)



    def launcher(self, input_file, slha_file, output_file, particle_1, particle_2, num_try, order, Xsections):
        #modifie_slha_file(input_file, input_file, slha_file)
        self.modifie_outgoing_particles(input_file, input_file, particle_1, particle_2)
        #on lance si c'est le premier essai par défaut
        
        hist = 0
        already_written_canal = self.canaux_finding(slha_file)
        
        if (particle_1, particle_2) in already_written_canal:
            return
        
        if num_try == 0:
            self.launch_command(self.resummino_bin, input_file, output_file, order)

        #Ici on écrit dans le fichier slha, la variable hist permet de voir s'il y a eu une erreur
        #Dans le calcul des section efficaces Lo et NLO
            hist = self.write_in_slha(output_file, slha_file, order, particle_1, particle_2, 'all', Xsections)

        #On vérifie si jamais on a écrit trop de choses
        self.are_crosssection(slha_file, order)

        #Si jamais il y a effectivement une erreur, on l'indique et on relance avec cette fois num_try = 0
        if hist == 1:
            print("error")
            num_try = 0
            self.modifie_outgoing_particles(input_file, input_file, particle_1, particle_2)
            self.launcher(input_file, slha_file, output_file, particle_1, particle_2, num_try, order)
    #launcher(input_file, slha_file, output_file, particle_1, particle_2)


    def search_in_output(self, output_file):
        Infos = []
        with open(output_file, 'r') as f: #Chemin a modifier si dossier différent pour les fichiers .txt
            data = f.readlines()
        for i in range(len(data)):
            if "Results:" in data[i]:
                LO = data[i+1][:-1] #[:-1] pour enlever le \n
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
            for i in range(len(result)):
                Xsection = crossSection.XSection()
        
                Xsection.value = float(result[i]) * pb
                Xsection.info.order = i
                Xsection.info.sqrts = 13. * TeV
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
            
    def write_in_slha(self, output_file, slha_file, order, particle_1, particle_2, type_writing, Xsections):
        results = self.search_in_output(output_file)
        if type_writing == 'max':
            if order == 0:
                result = results[0].split(" ")[2][1:]
            elif order == 1:
                result = results[1].split(" ")[2][1:]
            elif order == 2:
                result = results[2].split(" ")[2][1:]
        if type_writing == "all":
            result = [results[0].split(" ")[2][1:], results[1].split(" ")[2][1:], results[2].split(" ")[2][1:]]
        
        if order == 1:
            #_ = math.fabs(results[0].split(" ")[2][1:]-results[1].split(" ")[2][1:])
            if float(results[1].split(" ")[2][1:])>2*float(results[0].split(" ")[2][1:]) or float(results[0].split(" ")[2][1:])> 2* float(results[1].split(" ")[2][1:]):
                with open('log.txt', 'a') as f:
                    f.write(f"to much change between LO and NLO for {particle_1} and {particle_2} with {slha_file}")
                return 1
            
        self.create_xsection(result, particle_1, particle_2, order, Xsections)
            
        #self.writing_result(result, particle_1, particle_2, slha_file, order, type_writing)
        return 0

    def extract_m1_m2_mu(self, file_path):
        
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
        Vérifie si les sections efficaces sont déjà écrites, et supprime
        également les doublons (smodels s'en charge mais c'est toujours mieux)
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
                canal = (int(line.split(" ")[7]), int(line.split(" ")[8]),line.split(" ")[2])
                
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
                if line.split(" ")[2] == '1.30E+04':
                    canal = (int(line.split(" ")[7]), int(line.split(" ")[8]))
                    canaux.append(canal)
        return canaux

    #Choisi ici les canaux que tu souhaites utiliser
    def discrimination_particles(self, slha_file):
        """
        Choix des différents canaux et conditions sur les calculs de sections efficaces
        pour l'instant C1<92 or N1>600 or C2>1200
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
        #Fichier d'input de resummino de référence

        #Ne fait rien
        output_file = "output_file.txt"

        #Listes des differents fichiers sur lesquels faire tourner le logiciel
        Liste_slha = []
        Liste_resummino_in = []
        Liste_output_file = []
        Liste_particles = []
        Liste = []

        #Vérification des dossiers d'entrée et sortie de Resummino
        if not os.path.exists('resummino_in'):
            os.mkdir('resummino_in')
        if not os.path.exists('resummino_out'):
            os.mkdir('resummino_out')

        #On créer la liste des fichiers d'entrée
        liste_slha = os.listdir(slha_folder_name)

        #ancre
        pwd = os.getcwd()

        #Utilisation chemin absolu (à privilégier)
        slha_folder = os.path.join(pwd, slha_folder_name)

        #a = nombre fichiers, b = nombre fichiers déjà écris, c = nombre fichiers déjà ignorés
        a,b,c = 0,0,0

        #Boucle sur tous les fichiers, dans l'ordre de listdir (ordre +/- aléatoire)
        for slha in liste_slha:

            slha_path = os.path.join(slha_folder,slha)
            #Variable utilisée pour éviter les abérations (grosse différence LO/NLO)
            num_try = 0
            with open(slha_path, 'r') as f:
                data = f.readlines()
                a+=1
                
            if data[-1].endswith("Resumminov3.1.2\n"):
                b+=1
                #On augmente cette variable de 1, comme ca si elle est > 0 on ne refait pas le calcul
                num_try+=1
            elif data[-1].endswith(" #no_cross-section\n"):
                c+=1
                #On augmente cette variable de 1, comme ca si elle est > 0 on ne refait pas le calcul
                num_try+=1

            #on enlève le .slha
            slha_file_name = slha[6:-5]

            #On prend le fichier de référence, et on créer une copie dans resummino_in avec le bon fichier slha
            self.modifie_slha_file(self.input_file_original, f"resummino_in/resummino_{slha_file_name}.in",slha_path)

            #On ajoute les noms à la liste (in, out et slha)
            Liste_resummino_in.append(f"resummino_in/resummino_{slha_file_name}.in")
            Liste_output_file.append(f"resummino_out/resummino_out_{slha_file_name}.txt")
            Liste_slha.append(slha_path)

            #On liste ici les canaux à utiliser, si scénario exclu alors renvoi None
            particles = self.discrimination_particles(slha_path)
            Liste_particles.append(particles)

            #On pourrait optimiser en enlevant les variables qui ne changent pas d'une itération à l'autre
            #Mais ce n'est pas très important (négligeable niveau temps de compilation comparé à Resummino)
            Liste.append((particles, f"resummino_in/resummino_{slha_file_name}.in", slha_path, f"resummino_out/resummino_out_{slha_file_name}.txt", num_try, order))
        print(f" Number of files created : {a-b-c}")
        return Liste


    def modifie_slha_file(self, file_before, file_after, slha_file):
        with open(file_before, 'r') as f:
            lines = f.readlines()

        with open(file_after, 'w') as f:
            for line in lines:
                if line.startswith("slha ="):
                    line = f"slha = {slha_file}\n"
                f.write(line)

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
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(self.one_slha_calculation, *task) for task in tasks]
            for future in futures:
                future.result()

class ArgsStandardizer:
    """ simple class to collect all argument manipulators """

    def getSSMultipliers ( self, multipliers ):
        if type ( multipliers ) == str:
            if multipliers in [ "", "None", "none", "no", "{}" ]:
                return None
            if multipliers.count("{") != 1 or multipliers.count("}") != 1:
                logger.error ( "need to pass signal strengh multipliers as dictionary with tuple of pids as keys" )
            if multipliers.count("(") != multipliers.count(")"):
                logger.error ( "need to pass signal strengh multipliers as dictionary with tuple of pids as keys" )
            return eval(multipliers)
        return multipliers

    def getInputFiles ( self, args ):
        """ geth the names of the slha files to run over """
        inputPath  = args.filename.strip()
        if not os.path.exists( inputPath ):
            logger.error( "Path %s does not exist." % inputPath )
            sys.exit(1)
        inputFiles = []
        if os.path.isfile ( inputPath ):
            inputFiles = [ inputPath ]
        else:
            files = os.listdir ( inputPath )
            for f in files:
                inputFiles.append ( os.path.join ( inputPath, f ) )
        import random
        random.shuffle ( inputFiles )
        return inputFiles

    def checkAllowedSqrtses ( self, order, sqrtses ):
        """ check if the sqrtses are 'allowed' """
        if order == 0: return
        allowedsqrtses=[7, 8, 13, 13.6]
        for sqrts in sqrtses:
            if not sqrts in allowedsqrtses:
                logger.error("Cannot compute NLO or NLL xsecs for sqrts = %d "
                        "TeV! Available are: %s TeV." % (sqrts, allowedsqrtses ))
                sys.exit(-2)

    def getOrder ( self, args ):
        """ retrieve the order in perturbation theory from argument list """
        if args.NLL:
            return NLL
        if args.NLO:
            return NLO
        return LO

    def queryCrossSections ( self, filename ):
        if os.path.isdir ( filename ):
            logger.error ( "Cannot query cross sections for a directory." )
            sys.exit(-1)
        xsecsInfile = crossSection.getXsecFromSLHAFile(filename)
        if xsecsInfile:
            print ( "1" )
        else:
            print ( "0" )

    def getSqrtses ( self, args ):
        """ extract the sqrtses from argument list """
        sqrtses = [item for sublist in args.sqrts for item in sublist]
        if len(sqrtses) == 0:
            sqrtses = [8,13]
        sqrtses.sort()
        sqrtses = set(sqrtses)
        return sqrtses

    def checkNCPUs ( self, ncpus, inputFiles ):
        if ncpus < -1 or ncpus == 0:
            logger.error ( "Weird number of CPUs given: %d" % ncpus )
            sys.exit()
        if ncpus == -1:
            ncpus = runtime.nCPUs()
        ncpus = min ( len(inputFiles), ncpus )
        if ncpus == 1:
            logger.info ( "We run on a single cpu" )
        else:
            logger.info ( "We run on %d cpus" % ncpus )
        return ncpus

    def getPythiaVersion ( self, args ):
        pythiaVersion = 8

        if hasattr(args, 'pythia6' ) and args.pythia6 == True:
            pythiaVersion = 6
            if hasattr(args, 'pythia8') and args.pythia8 == True:
                logger.error ( "cannot both use pythia6 and pythia8 for LO xsecs." )
                sys.exit()
        return pythiaVersion

    def writeToFile ( self, args ):
        toFile = None
        if args.tofile:
            toFile="highest"
        if args.alltofile:
            if toFile=="highest":
                logger.warning ( "Specified both --tofile and --alltofile. Will use "\
                              "--alltofile" )
            toFile="all"
        return toFile

def main():
    
    test = XSecResummino(1, 'smodels/tools/exemple')
    test.routine_resummino()
    # canonizer = ArgsStandardizer()
    # setLogLevel ( args.verbosity )
    # if not hasattr ( args, "noautocompile" ):
    #     args.noautocompile = False
    # if args.query:
    #     return canonizer.queryCrossSections ( args.filename )
    # if args.colors:
    #     from smodels.tools.colors import colors
    #     colors.on = True
    # sqrtses = canonizer.getSqrtses ( args )
    # order = canonizer.getOrder ( args )
    # canonizer.checkAllowedSqrtses ( order, sqrtses )
    # inputFiles = canonizer.getInputFiles ( args )
    # ncpus = canonizer.checkNCPUs ( args.ncpus, inputFiles )
    # pythiaVersion = canonizer.getPythiaVersion ( args )
    # ssmultipliers = None
    # if hasattr ( args, "ssmultipliers" ):
    #     ssmultipliers = canonizer.getSSMultipliers ( args.ssmultipliers )
    #     if ssmultipliers != None:
    #         for pids,multiplier in ssmultipliers.items():
    #             if type(pids) != tuple:
    #                 logger.error ( "keys of ssmultipliers need to be supplied as tuples" )
    #                 sys.exit()
    #             if type(multiplier) not in [ int, float ]:
    #                 logger.error ( "values of ssmultipliers need to be supplied as ints or floats" )
    #                 sys.exit()

    # pythiacard = None
    # if hasattr(args, 'pythiacard'):
    #     pythiacard = args.pythiacard

    # children = []
    # for i in range(ncpus):
    #     pid = os.fork()
    #     chunk = inputFiles [ i::ncpus ]
    #     if pid < 0:
    #         logger.error ( "fork did not succeed! Pid=%d" % pid )
    #         sys.exit()
    #     if pid == 0:
    #         logger.debug ( "chunk #%d: pid %d (parent %d)." %
    #                    ( i, os.getpid(), os.getppid() ) )
    #         logger.debug ( " `-> %s" % " ".join ( chunk ) )
    #         computer = XSecResummino( order, args.nevents, pythiaVersion, \
    #                                  not args.noautocompile )
    #         toFile = canonizer.writeToFile ( args )
    #         computer.computeForBunch (  sqrtses, chunk, not args.keep,
    #                       args.LOfromSLHA, toFile, pythiacard=pythiacard, \
    #                     ssmultipliers = ssmultipliers )
    #         os._exit ( 0 )
    #     if pid > 0:
    #         children.append ( pid )
    # for child in children:
    #     r = os.waitpid ( child, 0 )
    #     logger.debug ( "child %d terminated: %s" % (child,r) )
    # logger.debug ( "all children terminated." )

main()