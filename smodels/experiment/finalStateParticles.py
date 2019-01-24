"""
.. module:: finalStateParticles
   :synopsis: Defines the final state particles used in the experimental results.
              It also makes use of the particles defined in smodels.share.models.SMparticles.py
   
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory.particle import Particle,MultiParticle
from smodels.share.models.SMparticles import SMList,e,eC,mu,muC,ta,taC,pion,t,tC,W,WC,g,q,c,nu
from smodels.theory.model import Model

###Below we define all particles that are relevant for the database constraints.
###One must avoid defining non-unique labels, so there is a one-to-one correspondence
###between the particles defined in finalStates and the labels used to describe the
###particles in the simplified model.


#Particle groups
eList = MultiParticle('e' , [e,eC])
muList = MultiParticle('mu' , [mu,muC])
taList = MultiParticle('ta' , [ta,taC])
lpList = MultiParticle('l+' , [eC,muC])
lmList = MultiParticle('l-'  , [e,mu])
lList = MultiParticle('l' , [e,mu,eC,muC])
WList = MultiParticle('W'  , [W,WC])
tList = MultiParticle('t'  , [t,tC])
LpList = MultiParticle('L+' , [eC,muC,taC])
LmList = MultiParticle('L-' , [e,mu,ta])
LList = MultiParticle('L'  , [e,mu,ta,eC,muC,taC] )
jetList = MultiParticle('jet' ,[q,c,g,pion])
nuList  = nu

#Used to construct generic z2-odd and z2-even particles:
anyOdd = Particle(label='anyOdd',Z2parity='odd')
anyEven = Particle(label='*',Z2parity='even')


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



#Define list of inclusive final states:
finalStates = [eList,muList,taList,lpList,lmList,lList,WList,
               tList,LpList,LmList,LList,jetList,anyEven]
#Include list of exclusive final states:
finalStates +=  SMList
#Define list of BSM final states:
BSMfinalStates = [MET,HSCP,RHadronG,RHadronQ,anyOdd]

#Avoid double counting:
for i,ptc in enumerate(finalStates):
    if any((ptc is p and i != j) for j,p in enumerate(finalStates)):
        finalStates.remove(ptc)
for i,ptc in enumerate(BSMfinalStates):
    if any((ptc is p and i != j) for j,p in enumerate(BSMfinalStates)):
        BSMfinalStates.remove(ptc)
            
#Protect all final state properties:
for ptc in finalStates+BSMfinalStates:
    ptc._static = True

finalStates = Model(SMparticles = finalStates,
                    BSMparticles=BSMfinalStates)

