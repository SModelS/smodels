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


from smodels.share.models import MSSMparticles
from smodels.share.models import SMparticles
from smodels.theory.particle import Particle,ParticleList,ParticleWildcard
from smodels.experiment import finalStateParticles

def add(ptcList,obj):
    """
    If obj is a Particle, ParticleList or list of particle objects,
    add the particle objects to ptcList.
    """
    
    if isinstance(obj,Particle) or isinstance(obj,ParticleWildcard):
        if not any(obj is ptc for ptc in ptcList):
            ptcList.append(obj)
    elif isinstance(obj,ParticleList):
        if not any(obj is ptc for ptc in ptcList):
            ptcList.append(obj)        
        for ptc in obj.particles:
            add(ptcList,ptc)
    elif isinstance(obj,list):
        for ptc in obj:
            add(ptcList,ptc)
        

# all particles and particle lists from the SM
SM = []
for pname in dir(SMparticles):
    obj =  getattr(SMparticles,pname)
    add(SM,obj)

    
BSM = []
for pname in dir(MSSMparticles):
    obj =  getattr(MSSMparticles,pname)
    add(BSM,obj)
    
finalStates = []
for pname in dir(finalStateParticles):
    obj =  getattr(finalStateParticles,pname)
    add(finalStates,obj)

#Store all particles:
allParticles = []
for particle in SM+BSM+finalStates:
    if not any(particle is ptc for ptc in allParticles):
        allParticles.append(particle)

