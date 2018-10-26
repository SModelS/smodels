"""
.. module:: finalStateParticles
   :synopsis: Defines the final state particles used in the experimental results.
              It also makes use of the particles defined in smodels.share.models.SMparticles.py
   
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory.particle import Particle,MultiParticle
from smodels.share.models.SMparticles import (chargedLeptons,chargedLeptonsC,
                                              leptons,leptonsC,gauge,gaugeC,
                                            quarks,quarksC,pip,piz,pim,pizz,higgs,higgsC)
from smodels.theory.model import Model
from smodels.theory.exceptions import SModelSTheoryError as SModelSError

#Particle groups
eList = MultiParticle('e' , [leptons[0], leptonsC[0]])
muList = MultiParticle('mu' , [leptons[1], leptonsC[1]])
taList = MultiParticle('ta' , [leptons[2], leptonsC[2]])
lpList = MultiParticle('l+' , [leptonsC[0],leptonsC[1]])
lmList = MultiParticle('l-'  , [leptons[0],leptons[1]]) 
lList = MultiParticle('l' , lpList.particles + lmList.particles ) 
nuList = MultiParticle('nu' , [leptons[3],leptons[4],leptons[5],leptonsC[3],leptonsC[4],leptonsC[5]])
WList = MultiParticle('W'  , [gauge[2],gaugeC[2]])
uList = MultiParticle('u'  , [quarks[0],quarksC[0]])
dList = MultiParticle('d'  , [quarks[1],quarksC[1]])
cList = MultiParticle('c'  , [quarks[2],quarksC[2]])
sList = MultiParticle('s'  , [quarks[3],quarksC[3]])
tList = MultiParticle('t'  , [quarks[4],quarksC[4]])
bList = MultiParticle('b'  , [quarks[5],quarksC[5]])
LpList = MultiParticle('L+' , lpList.particles + [ leptonsC[2] ])
LmList = MultiParticle('L-' , lmList.particles + [ leptons[2] ])
LList = MultiParticle('L'  , LpList.particles + LmList.particles )
jetList = MultiParticle('jet' ,  quarks[0:4] + [gauge[0]] + [pip,piz,pim,pizz] + quarksC[0:4] + [ gaugeC[0] ])
qList = MultiParticle('q' ,  quarks[0:4] +quarksC[0:4])
higgsList = MultiParticle('higgs' , [higgs,higgsC])



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


###Below we define all particles that are relevant for the database constraints.
###One must avoid defining non-unique labels, so there is a one-to-one correspondence
###between the particles defined in finalStates and the labels used to describe the
###particles in the simplified model.
###For instance, we define 'b' as a MultiParticle corresponding to [b+,b-],
###instead of defining it twice with repeated labels, since in the database b means any of the
###charge assignments.

#Define list of inclusive final states:
finalStates = [eList,muList,taList,lpList,lmList,lList,nuList,WList,uList,
               dList,sList,bList,cList,
               tList,LpList,LmList,LList,jetList,qList,higgsList,anyEven]
#Include list of exclusive final states:
finalStates += chargedLeptons + chargedLeptonsC + [quarks[4],quarksC[4]] + gauge + gaugeC
#Define list of BSM final states:
BSMfinalStates = [MET,HSCPp,HSCPm,HSCP,RHadronG,RHadronQ,anyOdd]

#Avoid double counting:
for i,ptc in enumerate(finalStates):
    if any((ptc is p and i != j) for j,p in enumerate(finalStates)):
        finalStates.remove(ptc)
for i,ptc in enumerate(BSMfinalStates):
    if any((ptc is p and i != j) for j,p in enumerate(BSMfinalStates)):
        BSMfinalStates.remove(ptc)

#Check for ambiguous label definitions:
allLabels = [p.label for p in finalStates+BSMfinalStates]
for p in finalStates+BSMfinalStates:
    counter = allLabels.count(p.label)
    if counter > 1:
        raise SModelSError("Particle with label %s has been defined %i times in finalStates. Make sure to use unique labels for final states."
                           %(p.label,counter))
            
#Protect all final state properties:
for ptc in finalStates+BSMfinalStates:
    ptc._static = True

finalStates = Model(SMparticles = finalStates,
                    BSMparticles=BSMfinalStates)

