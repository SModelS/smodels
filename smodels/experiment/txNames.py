#!/usr/bin/env python

"""
.. module:: txNames
   :synopsis: Holds a dictionary with decays for every txName and a small object to read them.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Michael Traub <michael.traub@gmx.at>

"""

import logging

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.ERROR)


class TxDecay(object):

    def __init__(self,txName):
        
        self._name = txName
        
    def __str__(self):
        
        return self._name
        
    def __nonzero__(self):
        
        return self._name in decays

        
    @property
    def shortdecay(self):
        return self._shortdecay
        
    @property
    def decay(self):
        return self._decay
        
    @property
    def motherParticle(self):
       return self._motherParticle
       
    @property
    def intermediateParticles(self):
        return self._intermediateParticles
        
        
    def _slackExpTopologyName(self):
        """Bypassing case sensitivity
        # ### FIX ME: doesn't know much at the moment. Is this still needed?
        """
        if any(c in self._name for c in ['w', 'W', 'z', 'Z']):
            return self._name.replace("W","w").replace("Z","z" )
        return self._name
        
    def _searchDecayDict(self):
        """Searches for topology name in descriptions.decay
        :returns: dictionary entry without formatting 
        
        """
        if decays.has_key(self._name):
            logger.info('found decay for topology %s' %self._name)
            return decays[self._name]
        if decays.has_key(self._slackExpTopologyName()):
            logger.info('found decay for topology %s with \
            slack name %s' %(self._name, self._slackExpTopologyName()))
            return decays[self._slackExpTopologyName()]
        logger.warning('no decay found for topology %s' %self._name)
        return None        

    @property
    def _decay(self):
        """:returns: decay as string, formatted for ROOT.TLatex
        
        """
        
        decay = self._searchDecayDict()
        if isinstance(decay,str): return self._latexDecay(decay)
        if isinstance(decay,list):
            i = 1
            lenght = len(decay)
            decayString =''
            for line in decay:
                if i != 1: decayString = decayString + '{'
                if i != lenght:
                    decayString = decayString +  '#splitline{' + self._latexDecay(line) 
                if i == lenght:
                    decayString = decayString + self._latexDecay(line)
                decayString = decayString + '}'
                if i == lenght:
                    decayString = decayString + '}'*(i-2)
                i += 1
            return decayString
            
    def _latexDecay(self, decayString):
        """Translates decay description as given in decays dictionary
        to a string readable by ROOT.TLatex object.
        
        """
        for key, value in prettySUSYParticle.items():
            decayString = self._latexParticle(decayString,key,value)
        for key, value in prettySMParticle.items():
            decayString = self._latexParticle(decayString,key,value)
        for key, value in highstrings.items():
            decayString = decayString.replace(key,value)
        for key, value in lowstrings.items():
            decayString = decayString.replace(key,value)
        decayString = decayString.replace('-->','#rightarrow')
        return decayString
        
    def _latexParticle(self,decayString,key,value):
        """Translates particle description as given in decays dictionary
        to a string readable by ROOT.TLatex object.
        
        """
        decayString = decayString.replace('anti' + key + ' ','#bar{' + value + '}')
        decayString = decayString.replace(key + ' ',value)
        decayString = decayString.replace(key + '_',value + '_')
        decayString = decayString.replace(key + '^',value + '^')
        return decayString
        
    @property
    def _motherParticle(self):
        """ :returns: mother particle in simple format as string or None
        
        """
        # ### FIX ME: This is not done yet! 
        decay = self._searchDecayDict()
        if isinstance(decay,list): decay = decay[0]
        motherPart = decay.split('-->')[0]
        motherPart = motherPart.strip()
        if motherPart == 'gluino': return 'g'
        if motherPart == 'squark': return 'q'
        if motherPart == 'stop': return 't'
        if motherPart == 'sbottom': return 'b'
        if motherPart == 'slepton': return 'l'
        if motherPart == 'stop_2': return 't2'
        if 'gluino' in motherPart and 'squark' in motherPart:
            return 'gq'
        if 'chargino' in motherPart and 'neutralino' in motherPart:
            return 'c0cpm'
        if 'chargino' in motherPart and not 'neutralino' in motherPart:
            return 'cpm'
        if not 'chargino' in motherPart and 'neutralino' in motherPart \
        and not 'neutralino_2' in motherPart and not 'neutralino_3' in motherPart:
            return 'c0'
        if  'neutralino_2' in motherPart and 'neutralino_3' in motherPart:#
            return 'c02c03'
        logger.error('could not identify mother particle for  %s' %self._name)
        return None
        
    @property
    def shortdecay(self):
        """:returns: short version of decay as string
        
        """
        
        decay = self._searchDecayDict()
        if isinstance(decay,list): decay = decay[0]
        decaySteps = decay.split('-->')
        if len(decaySteps) == 2: return self._latexDecay(decay)
        decay = decaySteps[0] + '--> '
        lsp = 'lsp '*(len(decay.split())-1)
        decaySteps = decaySteps[1:]
        for decayStep in decaySteps:
            #decayStep.replace('(','')
            #decayStep.replace(')','')
            #decayStep.replace('|','')
            for particle in prettySMParticle:
                if particle in decayStep: decay = decay + particle + ' '
        decay = decay + lsp
        return self._latexDecay(decay)
        
    @property
    def _intermediateParticles(self):
        """:returns: dictionary with intermediate particles
        
        """
        particles = []
        decays = self._searchDecayDict()
        if isinstance(decays,str): decays = [decays]
        for decay in decays:
            decay = decay.split('-->')
            decay = decay[1:-1]
            if not decay: continue
            for expression in decay:
                expression = expression.replace('(','')
                expression = expression.replace(')','')
                expression = expression.replace('|','')
                expression = expression.replace('lsp','')
                [particles.append(particle.strip()) for particle in expression.split(' ')]
        if not particles: return
        interParticles = []
        for particle in particles:
            for sparticle in prettySUSYParticle:
                if sparticle in particle and not particle in interParticles: 
                    interParticles.append(particle)
        return interParticles
        
# dictionary containing all decays for the different topologies.
# special format expected: 
# -supported particle names can be found in the dictionaries prettySMParticels and 
#  prettySUSYParticle
# -supported postfixes can be found in the dictionaries highstrings and lowstrings
# -prefix "anti" for anti-particles supported
# -space is expected to end a particle description (eg: chargino^pm_2 )
# -if there is more then one possible decay for one topology a list with decays is expected
# -neutrino, lepton, neutralion, slepton means first 2 generations
# -Neutrino, Lepton, Neutralion, sLepton means all 3 generations

decays = { 
    'T1': 'gluino  --> quark antiquark  lsp ' ,
    'T1bbbb': 'gluino  --> bottom antibottom  lsp ', 
    'T1tttt': 'gluino  --> top antitop  lsp ',
    'T1ttttoff': 'gluino  --> top antitop  lsp ',
    'T1gg':'gluino  --> quark antiquark (neutralino_2 --> photon lsp )', 
    'T1lg':'gluino  --> quark antiquark (neutralino_2  --> photon lsp |chargino^pm  --> w lsp )', 
    'T1lnu':'gluino  --> quark antiquark (chargino^pm --> lepton^pm neutrino  lsp )', 
    'T1lh':'gluino  --> quark antiquark  neutralino_2 neutralino_2  --> lepton^p lepton^m lsp ', 
    'T2':'squark  --> quark lsp ',
    'T2FVttcc': 'stop  --> charm lsp ',
    'T2llnunubb': 'stop  --> lepton neutrino bottom lsp ',
    'T2bb':'sbottom  --> bottom lsp ', 
    'T2bw':'stop  --> bottom w lsp ',
    'T2ttww': 'sbottom  --> top w lsp ',
    'T2bbWW': 'stop  --> bottom w lsp ',
    'T2tt': 'stop  --> top lsp ', 
    'T3w': 'gluino --> quark antiquark (chargino^pm_1 --> w lsp | lsp )' ,
    'T3wb':'gluino  --> bottom antibottom (w )lsp ', 
    'T3lh':'gluino  --> quark antiquark (neutralino_2 --> lepton^p lepton^m lsp | lsp )',
    'T3tauh':'gluino  --> quark antiquark (neutralino_2 --> tau tau lsp | lsp )', 
    'T5WW':'gluino  --> quark antiquark (chargino^pm_1 --> w lsp )',
    'T5wg':'gluino  --> quark antiquark (neutralino_2 --> photon lsp | chargino^pm_1 --> w lsp )',
    'T5WH':'gluino  --> quark antiquark (neutralino_2 --> higgs lsp | chargino^pm_1 --> w lsp )',
    'T5gg':'gluino  --> quark antiquark (neutralino_2 --> photon lsp )',
    'T5lnu':'gluino  --> quark antiquark (chargino^pm --> lepton^pm neutrino lsp )',
    'T5ZZ':'gluino  --> quark antiquark (neutralino_2 --> z lsp )',
    'T5ZZInc':'neutralino_2 --> z lsp ',
    'T5zzgmsb':'gluino --> quark antiquark (neutralino_2 --> z lsp )', 
    'T5tttt':'gluino  --> top (stop --> top antitop lsp )',
    'T6ttww': 'sbottom  --> top (chargino^pm_1 --> w lsp )',
    'T6WW': 'squark  --> quark (chargino^pm_1 --> w lsp )',
    'T6ttHH': 'stop  --> top higgs lsp ',
    'T6ZZtt': 'stop_2  --> z (stop_1 --> top lsp ) ',
    'T6bbWW':'stop  --> bottom (chargino^p --> w lsp )',
    'T6bbWWoff':'stop  --> bottom (chargino^p --> w lsp )',
    'T6bbZZ':'sbottom  -->  bottom (neutralino_2 --> z lsp )',
    'T7btW':'gluino  --> bottom top w lsp ',
    'T7btbtWW':'gluino  --> bottom (sbottom --> top (chargino^pm --> w lsp ))',
    'TGQ':'gluino squark --> quark antiquark lsp | quark lsp ',
    'TChizz':'neutralino_3 neutralino_2  --> z z lsp lsp ',
    'TChiSlep':'neutralino_2 chargino^pm_1  --> lepton lepton lepton neutrino lsp lsp ',
    'TChiNuSlep':'neutralino_2 chargino^pm_1  --> lepton lepton lepton neutrino lsp lsp ',
    'TChizz':'neutralino_3 neutralino_2  --> z z lsp lsp ',
    'TChiwz':'chargino^pm neutralino_2  --> w z lsp lsp ',
    'TChiWZon':'chargino^pm neutralino_2  --> w z lsp lsp ',
    'TChiWH':'chargino^pm neutralino_2  --> w higgs lsp lsp ',
    'TChiWW':'chargino^pm  --> w lsp ',
    'TChiWZoff':'chargino^pm neutralino_2  --> w z lsp lsp ',
    'TChiChipmSlepSlep':'neutralino_2 chargino^pm_1  --> lepton (slepton --> lepton lsp ) | neutrino (slepton --> lepton lsp )',
    'TChiChipmStauStau':'neutralino_2 chargino^pm_1  --> tau (stau --> tau lsp ) | neutrino (stau --> tau lsp )', 
    'TChiChipmSlepL':[
    'neutralino_2 chargino^pm_1  --> Lepton (sLepton --> Lepton lsp ) | Neutrino (sLepton --> Lepton lsp )' ,
    'neutralino_2 chargino^pm_1  --> Lepton (sLepton --> Lepton lsp ) | Lepton (sNeutrino --> Neutrino lsp )',
    'neutralino_2 chargino^pm_1  --> Neutrino (sNeutrino --> Neutrino lsp ) | Lepton (sNeutrino --> Neutrino lsp )' ,
    'neutralino_2 chargino^pm_1  --> Neutrino (sNeutrino --> Neutrino lsp ) | Neutrino (sLepton --> Lepton lsp )'
    ], 
    'TChiChiSlepSlep': 'neutralino_2 neutralino_3  --> lepton (slepton --> lepton lsp )',
    'TChiChipmStauL':[
    'neutralino_2 chargino^pm_1  --> tau (stau --> tau lsp ) | tau-neutrino (stau --> tau lsp )',
    'neutralino_2 chargino^pm_1  --> tau (stau --> tau lsp ) | tau (stau --> tau-neutrino lsp )',
    'neutralino_2 chargino^pm_1  --> neutrino (tau-sneutrino --> neutrino lsp ) | tau (tau-sneutrino --> tau-neutrino lsp )',
    'neutralino_2 chargino^pm_1  --> neutrino (tau-sneutrino --> neutrino lsp ) | tau-neutrino (stau -->  lsp )'
    ], 
    'TChiChipmHW':'neutralino_2 chargino^pm_1  --> w lsp higgs lsp ', 
    'TChiChipmSlepStau':'neutralino_2 chargino^pm_1  --> Lepton (slepton --> Lepton lsp ) | neutrino (stau --> tau lsp )', 
    'TChiChipmStauStau':'neutralino_2 chargino^pm_1  --> tau (stau --> tau lsp ) | neutrino (stau --> tau lsp ) ',
    'TChipChimSlepSnu':[
    'chargino^pm chargino^mp  --> Lepton (sneutrino --> neutrino lsp ) | neutrino (slepton --> Lepton lsp )',
    'chargino^pm chargino^mp  --> Lepton (sneutrino --> neutrino lsp ) | Lepton (sneutrino --> neutrino lsp )',
    'chargino^pm chargino^mp  --> neutralino (slepton --> Lepton lsp ) | neutrino (slepton --> Lepton lsp )'
    ], 
    'TSlepSlep':'slepton  --> lepton lsp '
}

prettySMParticle = {
    'graviton':'#tilde{G}',         #graviton
    'photon': '#gamma',             #photon
    'gluon':'g',                    #gluon
    'w' : 'W',                  #W
    'z' : 'Z',                  #Z
    'higgs' : 'H',                  #higgs
    
    'quark': 'q',           #quark
    'up': 'u',           #up
    'down': 'd',           #down
    'charm': 'c',           #charm
    'strange': 's',           #strange
    'top': 't',           #top
    'bottom': 'b',           #bottom
    
    'lepton' : 'l',             #lepton (first 2 generations)
    'Lepton' : 'l',             #lepton (all 3 generations)
    'electron' : 'e',               #electron
    'mu' : '#mu',            #mu
    'tau' : '#tau',  #tau
    
    'neutrino' : '#nu',                     #neutrino (first 2 generations)
    'Neutrino' : '#nu',                     #neutrino (all 3 generations)
    'electron-neutrino' : '#nu_{e}',               #electron-neutrino
    'mu-neutrino' : '#nu_{#mu}',            #mu-neutrino
    'tau-neutrino' : '#nu_{#tau}',          #tau-neutrino
}

prettySUSYParticle = {
    'lsp' : '#tilde{#chi}^{0}_{1}',  # lightesd SUSY particle
    'neutralino' : '#chi^{0}',      #neutralino
    'chargino' : '#chi',            #Chargino
    'gravitino':'G',              #gravitino
    'photino':'#tilde{#gamma}',   #photino
    'gluino': '#tilde{g}',        #gluino
    'wino' : '#tilde{W}',       #Wino
    'zino' : '#ti:lde{Z}',       #Zino
    'higgsino' : '#tilde{H}',       #higgsino
    
    'squark': '#tilde{q}',  #squarkfound in
    'sup': '#tilde{u}',  #sup
    'sdown': '#tilde{d}',  #sdown
    'scharm': '#tilde{c}',  #scharm
    'sstrange': '#tilde{s}',  #sstrange
    'stop': '#tilde{t}',  #stop
    'sbottom': '#tilde{b}',  #sbottom
    
    'slepton' : '#tilde{l}',    #slepton (first 2 generations)
    'sLepton' : '#tilde{l}',    #sLepton (all 3 generations)
    'selectron' : '#tilde{e}',      #selectron
    'smu' : '#tilde{#mu}',   #smu
    'stau' : '#tilde{#tau}', #stau
    
    'sneutrino' : '#tilde{#nu}',            #sneutrino (first 2 generations)
    'sNeutrino' : '#tilde{#nu}',            #sneutrino (all 3 generations)
    'electron-sneutrino' : '#tilde{#nu}_{e}',      #electron-sneutrino
    'mu-sneutrino' : '#tilde{#nu}_{#mu}',   #mu-sneutrino
    'tau-sneutrino' : '#tilde{#nu}_{#tau}', #tau-sneutrino  
}

highstrings = {
    '^0' : '^{0}',
    '^pm' : '^{#pm}',
    '^mp' : '^{#mp}',
    '^p' : '^{+}',
    '^m' : '^{-}',
}

lowstrings = {
    '_1' : '_{1}',
    '_2' : '_{2}',
    '_3' : '_{3}',
    '_4' : '_{4}',
}
            
