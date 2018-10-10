"""
.. module:: finalStateParticles
   :synopsis: Defines the final state particles used in the experimental results.
              It also makes use of the particles defined in smodels.share.models.SMparticles.py
   
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory.particle import InclusiveParticle,Particle,MultiParticle
from smodels.share.models.SMparticles import *
# from smodels.share.models.SMparticles import leptons,quarks,quarksC,leptonsC,gauge,gaugeC,pi
from smodels.theory.model import Model

#Particle groups
eList = MultiParticle('e' , [leptons[0], leptonsC[0]])
muList = MultiParticle('mu' , [leptons[1], leptonsC[1]])
taList = MultiParticle('ta' , [leptons[2], leptonsC[2]])
lpList = MultiParticle('l+' , [leptonsC[0],leptonsC[1]])
lmList = MultiParticle('l-'  , [leptons[0],leptons[1]]) 
lList = MultiParticle('l' , lpList.particles + lmList.particles ) 
nuList = MultiParticle('nu' , [leptons[3],leptons[4],leptons[5],leptonsC[3],leptonsC[4],leptonsC[5]])
WList = MultiParticle('W'  , [gauge[2],gaugeC[2]])
tList = MultiParticle('t'  , [quarks[4],quarksC[4]])
LpList = MultiParticle('L+' , lpList.particles + [ leptonsC[2] ])
LmList = MultiParticle('L-' , lmList.particles + [ leptons[2] ])
LList = MultiParticle('L'  , LpList.particles + LmList.particles )
jetList = MultiParticle('jet' ,  quarks[0:4] + [gauge[0]] + [pip,piz] + quarksC[0:4] + [ gaugeC[0] ] 
                       + [piz.chargeConjugate('pi'),pip.chargeConjugate('pi')])



#Used to construct generic z2-odd and z2-even particles:
anyOdd = InclusiveParticle(label='anyOdd',Z2parity='odd')
anyEven = InclusiveParticle(label='*',Z2parity='even')


#Used to construct BSM final states:
MET = Particle(label='MET', Z2parity = 'odd', eCharge = 0, colordim = 1)
HSCPp = Particle(label='HSCP+', Z2parity = 'odd', eCharge = +1, colordim = 1)
HSCPm = Particle(label='HSCP-', Z2parity = 'odd', eCharge = -1, colordim = 1)
HSCP = MultiParticle(label='HSCP', particles = [HSCPp,HSCPm])

RHadronG = Particle(label='RHadronG', Z2parity = 'odd', eCharge = 0, colordim = 8)
RHadronU = Particle(label='RHadronU', Z2parity = 'odd', eCharge = 2./3., colordim = 3)
RHadronD = Particle(label='RHadronD', Z2parity = 'odd', eCharge = -1./3., colordim = 3)
RHadronQ = MultiParticle(label='RHadronQ', particles = [RHadronU,RHadronU.chargeConjugate(),
                                                                         RHadronD,RHadronD.chargeConjugate()])

#Get all objects defined so far:
objects = list(locals().values())[:]
    
allFinalStates = []
#Get all particles defined here:
for obj in objects:
    if isinstance(obj,list):
        allFinalStates += [ptc for ptc in obj if isinstance(ptc,(Particle,MultiParticle,InclusiveParticle)) and 
                           not any(obj is x for x in allFinalStates)]
    elif isinstance(obj,(Particle,MultiParticle,InclusiveParticle)):
        if not any(obj is x for x in allFinalStates):
            allFinalStates.append(obj)
            
#Protect all final state properties:
for ptc in allFinalStates:
    ptc._static = True
    
finalStates = Model(SMparticles = allFinalStates, BSMparticles=[])

