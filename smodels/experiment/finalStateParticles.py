"""
.. module:: finalStateParticles
   :synopsis: Defines the final state particles used in the experimental results.
              It also makes use of the particles defined in smodels.share.models.SMparticles.py
   
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory.particle import InclusiveParticle,Particle,ParticleList
from smodels.share.models.SMparticles import *
# from smodels.share.models.SMparticles import leptons,quarks,quarksC,leptonsC,gauge,gaugeC,pi
from smodels.theory.model import Model

#Particle groups
eList = ParticleList('e' , [leptons[0], leptonsC[0]])
muList = ParticleList('mu' , [leptons[1], leptonsC[1]])
taList = ParticleList('ta' , [leptons[2], leptonsC[2]])
lpList = ParticleList('l+' , [leptonsC[0],leptonsC[1]])
lmList = ParticleList('l-'  , [leptons[0],leptons[1]]) 
lList = ParticleList('l' , lpList.particles + lmList.particles ) 
nuList = ParticleList('nu' , [leptons[3],leptons[4],leptons[5],leptonsC[3],leptonsC[4],leptonsC[5]])
WList = ParticleList('W'  , [gauge[2],gaugeC[2]])
tList = ParticleList('t'  , [quarks[4],quarksC[4]])
LpList = ParticleList('L+' , lpList.particles + [ leptonsC[2] ])
LmList = ParticleList('L-' , lmList.particles + [ leptons[2] ])
LList = ParticleList('L'  , LpList.particles + LmList.particles )
jetList = ParticleList('jet' ,  quarks[0:4] + [gauge[0]] + [pip,piz] + quarksC[0:4] + [ gaugeC[0] ] 
                       + [piz.chargeConjugate('pi'),pip.chargeConjugate('pi')])



#Used to construct generic allParticles (z2-odd) and SM (z2-even) particles:
anyBSM = InclusiveParticle(label='anyBSM',Z2parity='odd')
anySM = InclusiveParticle(label='*',Z2parity='even')


#Used to construct BSM final states:
MET = Particle(label='MET', Z2parity = 'odd', eCharge = 0, colordim = 1)
HSCPp = Particle(label='HSCP+', Z2parity = 'odd', eCharge = +1, colordim = 1)
HSCPm = Particle(label='HSCP-', Z2parity = 'odd', eCharge = -1, colordim = 1)
HSCP = ParticleList(label='HSCP', particles = [HSCPp,HSCPm])

RHadronG = Particle(label='RHadronG', Z2parity = 'odd', eCharge = 0, colordim = 8)
RHadronU = Particle(label='RHadronU', Z2parity = 'odd', eCharge = 2./3., colordim = 3)
RHadronD = Particle(label='RHadronD', Z2parity = 'odd', eCharge = -1./3., colordim = 3)
RHadronQ = ParticleList(label='RHadronQ', particles = [RHadronU,RHadronU.chargeConjugate(),
                                                                         RHadronD,RHadronD.chargeConjugate()])

#Get all objects defined so far:
objects = list(locals().values())[:]
    
allFinalStates = []
#Get all particles defined here:
for obj in objects:
    if isinstance(obj,list):
        allFinalStates += [ptc for ptc in obj if isinstance(ptc,(Particle,ParticleList,InclusiveParticle)) and 
                           not any(obj is x for x in allFinalStates)]
    elif isinstance(obj,(Particle,ParticleList,InclusiveParticle)):
        if not any(obj is x for x in allFinalStates):
            allFinalStates.append(obj)
            
#Protect all final state properties:
for ptc in allFinalStates:
    ptc._static = True
    
finalStates = Model(SMparticles = allFinalStates, BSMparticles=[])

