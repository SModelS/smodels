"""
.. module:: finalStateParticles
   :synopsis: Defines the final state particles used in the experimental results.
              It also makes use of the particles defined in smodels.share.models.SMparticles.py

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.base.particle import Particle, MultiParticle
from smodels.share.models.SMparticles import SMList, e, eC, mu, muC, ta, taC, pion, t, tC, W, WC, g, q, c, nu
from smodels.base.model import Model

# Below we define all particles that are relevant for the database constraints.
# One must avoid defining non-unique labels, so there is a one-to-one correspondence
# between the particles defined in finalStates and the labels used to describe the
# particles in the simplified model.


# Particle groups
eList = MultiParticle('e', [e, eC])
muList = MultiParticle('mu', [mu, muC])
taList = MultiParticle('ta', [ta, taC])
lpList = MultiParticle('l+', [eC, muC])
lmList = MultiParticle('l-', [e, mu])
lList = MultiParticle('l', [e, mu, eC, muC])
WList = MultiParticle('W', [W, WC])
tList = MultiParticle('t', [t, tC])
LpList = MultiParticle('L+', [eC, muC, taC])
LmList = MultiParticle('L-', [e, mu, ta])
LList = MultiParticle('L', [e, mu, ta, eC, muC, taC])
jetList = MultiParticle('jet', [q, c, g, pion])
nuList = nu

# Used to construct generic z2-odd and z2-even particles:
anyBSM = Particle(label='anyBSM', isSM=False)
anySM = Particle(label='anySM', isSM=True)
anyParticle = Particle(label='*')

# Used to construct BSM final states:
MET = Particle(label='MET', isSM=False, eCharge=0, colordim=1)
HSCPp = Particle(label='HSCP+', isSM=False, eCharge=+1, colordim=1)
HSCPm = Particle(label='HSCP-', isSM=False, eCharge=-1, colordim=1)
HSCP = MultiParticle(label='HSCP', particles=[HSCPp, HSCPm])

RHadronG = Particle(label='RHadronG', isSM=False, eCharge=0, colordim=8)
RHadronU = Particle(label='RHadronU', isSM=False, eCharge=2./3., colordim=3)
RHadronD = Particle(label='RHadronD', isSM=False, eCharge=-1./3., colordim=3)
RHadronQ = MultiParticle(label='RHadronQ', particles=[RHadronU, RHadronU.chargeConjugate(),
                                                      RHadronD, RHadronD.chargeConjugate()])

gluino = Particle(label='gluino', isSM=False, eCharge=0, colordim=8, spin=1./2.)
chargino = Particle(label='C1+', isSM=False, eCharge=1, colordim=1, spin=1./2.)
charginoBar = Particle(label='C1-', isSM=False, eCharge=-1, colordim=1, spin=1./2.)
C1 = MultiParticle(label='C1', particles=[chargino, charginoBar])
Hp = Particle(label='H+', isSM=False, eCharge=1, colordim=1, spin=0.)
Hm = Particle(label='H-', isSM=False, eCharge=-1, colordim=1, spin=0.)
Hpm = MultiParticle(label='Hpm', particles=[Hp, Hm])


# Define list of inclusive final states:
finalStates = [eList, muList, taList, lpList, lmList, lList, WList,
               tList, LpList, LmList, LList, jetList, anySM, anyParticle]
# Include list of exclusive final states:
finalStates += SMList
# Define list of BSM final states:
BSMfinalStates = [MET, HSCP, RHadronU, RHadronD, RHadronG, RHadronQ, anyBSM,
                  gluino, chargino, charginoBar, C1, Hp, Hm, Hpm]


# Avoid double counting:
for i, ptc in enumerate(finalStates):
    if any((ptc is p and i != j) for j, p in enumerate(finalStates)):
        finalStates.remove(ptc)
for i, ptc in enumerate(BSMfinalStates):
    if any((ptc is p and i != j) for j, p in enumerate(BSMfinalStates)):
        BSMfinalStates.remove(ptc)

# Protect all final state properties:
for ptc in finalStates+BSMfinalStates:
    ptc._static = True

finalStates = Model(SMparticles=finalStates,
                    BSMparticles=BSMfinalStates, label='DB Final States (default)')
