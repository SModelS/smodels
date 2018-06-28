"""
.. module:: finalStateParticles
   :synopsis: Defines the final state particles used in the experimental results.
              It also makes use of the particles defined in smodels.share.models.SMparticles.py
   
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory.particle import ParticleWildcard,Particle,ParticleList
from smodels.share.models.SMparticles import leptons,quarks,leptonsC,quarksC,gauge,gaugeC,pi



#Particle groups
eList = ParticleList('e' , [leptons[0], leptonsC[0]])
muList = ParticleList('mu' , [leptons[2], leptonsC[2]])
taList = ParticleList('ta' , [leptons[4], leptonsC[4]])
lpList = ParticleList('l+' , [leptonsC[0],leptonsC[2]])
lmList = ParticleList('l-'  , [leptons[0],leptons[2]]) 
lList = ParticleList('l' , lpList.particles + lmList.particles ) 
nuList = ParticleList('nu' , [leptons[1],leptons[3],leptons[5],leptonsC[1],leptonsC[3],leptonsC[5]])
WList = ParticleList('W'  , [gauge[2],gaugeC[2]])
tList = ParticleList('t'  , [quarks[4],quarksC[4]])
LpList = ParticleList('L+' , lpList.particles + [ leptonsC[4] ])
LmList = ParticleList('L-' , lmList.particles + [ leptons[4] ])
LList = ParticleList('L'  , LpList.particles + LmList.particles )
jetList = ParticleList('jet' ,  quarks[0:4] + [ gauge[0] ] + [ pi ] + quarksC[0:4] + [ gaugeC[0] ] + [ pi.chargeConjugate() ])



#Used to construct generic allParticles (z2-odd) and SM (z2-even) particles:
anyBSM = ParticleWildcard(label='anyBSM',Z2parity='odd')
anySM = ParticleWildcard(label='*',Z2parity='even')


#Used to construct BSM final states:
MET = Particle(label='MET', Z2parity = 'odd', eCharge = 0, colordim = 0)
HSCPp = Particle(label='HSCP+', Z2parity = 'odd', eCharge = +1, colordim = 0)
HSCPm = Particle(label='HSCP-', Z2parity = 'odd', eCharge = -1, colordim = 0)
HSCP = ParticleList(label='HSCP', particles = [HSCPp,HSCPm])

RHadronG = Particle(label='RHadronG', Z2parity = 'odd', eCharge = 0, colordim = 8)
RHadronU = Particle(label='RHadronU', Z2parity = 'odd', eCharge = 2./3., colordim = 3)
RHadronD = Particle(label='RHadronD', Z2parity = 'odd', eCharge = -1./3., colordim = 3)
RHadronQ = ParticleList(label='RHadronQ', particles = [RHadronU,RHadronU.chargeConjugate(),
                                                                         RHadronD,RHadronD.chargeConjugate()])

