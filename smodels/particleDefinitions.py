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


from smodels.SMparticleDefinitions import SMList, SMpdgs, SMnames, particleLists
from smodels.MSSMparticleDefinitions import BSMList, BSMpdgs, BSMnames


# all particles and particle lists from the SM
SM = SMList + particleLists

