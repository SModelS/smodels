"""
.. module:: particles
   :synopsis: Defines the particles to be used.
              All particles appearing in the model as well as the SM particles
              must be defined (or imported) here.
              Only the particles stored in these lists are used by
              external methods.
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
      
   HOW TO ADD NEW PARTICLES: simply add a new Particle object (or ParticleList object).
   Define all known properties and set other properties to None.
   
   Properties defined by the LHE or SLHA input file (such as masses, width and BRs)
   are automatically added later.
"""


from smodels import SMparticleDefinitions, MSSMparticleDefinitions
from smodels.theory.particle import Particle,ParticleList,ParticleWildcard
from smodels.tools.wildCards import ValueWildcard

def add(ptcList,obj):
    """
    If obj is a Particle, ParticleList or list of particle objects,
    add the particle objects to ptcList.
    """
    
    if isinstance(obj,Particle) or isinstance(obj,ParticleWildcard):
        ptcList.append(obj)
    elif isinstance(obj,ParticleList):
        ptcList.append(obj)        
        for ptc in obj.particles:
            add(ptcList,ptc)
    elif isinstance(obj,list):
        for ptc in obj:
            add(ptcList,ptc)
        

# all particles and particle lists from the SM
SM = []
for pname in dir(SMparticleDefinitions):
    obj =  getattr(SMparticleDefinitions,pname)
    add(SM,obj)

    
BSM = []
for pname in dir(MSSMparticleDefinitions):
    obj =  getattr(MSSMparticleDefinitions,pname)
    add(BSM,obj)
    

#Used to construct generic allParticles (z2-odd) and SM (z2-even) particles:
anyBSM = ParticleWildcard(label='anyBSM',Z2parity='odd', mass =ValueWildcard())
anySM = ParticleWildcard(label='*',Z2parity='even')

#Used to construct final states:

MET = ParticleList(label='MET', particles = [sparticle for sparticle in BSM 
                                if isinstance(sparticle,Particle) and sparticle.eCharge == 0 and sparticle.colordim == 0
                                and sparticle.Z2parity == 'odd'],
                                mass = ValueWildcard())
HSCP = ParticleList(label='HSCP', particles = [sparticle for sparticle in BSM 
                                if isinstance(sparticle,Particle) and abs(sparticle.eCharge) == 1 and sparticle.colordim == 0
                                and sparticle.Z2parity == 'odd'],
                                mass = ValueWildcard())
RHadronG = ParticleList(label='RHadronG', particles = [sparticle for sparticle in BSM 
                                if isinstance(sparticle,Particle) and sparticle.eCharge == 0 and sparticle.colordim == 8
                                and sparticle.Z2parity == 'odd'],
                                mass = ValueWildcard())
RHadronQ = ParticleList(label='RHadronQ', particles = [sparticle for sparticle in BSM 
                                if isinstance(sparticle,Particle) and abs(sparticle.eCharge) in [2./3.,1./3.] and sparticle.colordim == 3
                                and sparticle.Z2parity == 'odd'],
                                mass = ValueWildcard())


#Store all particles:
allParticles = BSM + SM + [anyBSM,anySM,MET,HSCP,RHadronG,RHadronQ]

