"""
.. module:: databaseParticles
   :synopsis: Defines the final state particles used in the experimental results.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.base.particle import Particle, MultiParticle
from smodels.base.model import Model
from smodels.base.physicsUnits import MeV, GeV
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError

# Below we define all particles that are relevant for the database constraints.
# One must avoid defining non-unique labels, so there is a one-to-one correspondence
# between the particles defined in finalStates and the labels used to describe the
# particles in the simplified model.

#  SM particles

# Charged leptons:
e = Particle(isSM=True, label='e-', pdg=11, mass=0.5*MeV, eCharge=-1, colordim=1, spin=1./2, totalwidth=0.*GeV)
mu = Particle(isSM=True, label='mu-', pdg=13, mass=106.*MeV, eCharge=-1, colordim=1, spin=1./2, totalwidth=0.*GeV)
ta = Particle(isSM=True, label='ta-', pdg=15, mass=1777.*MeV, eCharge=-1, colordim=1, spin=1./2, totalwidth=0.*GeV)
eC = e.chargeConjugate()
muC = mu.chargeConjugate()
taC = ta.chargeConjugate()
# -----------------------------------------------------------------------------------------------------------------------------
# Neutrinos
nue = Particle(isSM=True, label='nue', pdg=12, mass=0.*MeV, eCharge=0, colordim=1, spin=1./2, totalwidth=0.*GeV)
numu = Particle(isSM=True, label='numu', pdg=14, mass=0.*MeV, eCharge=0, colordim=1, spin=1./2, totalwidth=0.*GeV)
nuta = Particle(isSM=True, label='nuta', pdg=16, mass=0.*MeV, eCharge=0, colordim=1, spin=1./2, totalwidth=0.*GeV)
# Group all neutrinos into a single particle:
nu = MultiParticle('nu', [nue, numu, nuta, nue.chargeConjugate(), numu.chargeConjugate(), nuta.chargeConjugate()])
# -----------------------------------------------------------------------------------------------------------------------------

# Light quarks:
d = Particle(isSM=True, label='d', pdg=1, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, totalwidth=0.*GeV)
u = Particle(isSM=True, label='u', pdg=2, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, totalwidth=0.*GeV)
s = Particle(isSM=True, label='s', pdg=3, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, totalwidth=0.*GeV)
# Group all light quarks in a single particle
q = MultiParticle('q', [u, d, s, u.chargeConjugate(), d.chargeConjugate(), s.chargeConjugate()])
c = Particle(isSM=True, label='c', pdg=4, mass=1.28*GeV, eCharge=(2./3.), colordim=3, spin=1./2, totalwidth=0.*GeV)
# Group c and c* in a single particle:
c = MultiParticle('c', [c, c.chargeConjugate('c')])
# -----------------------------------------------------------------------------------------------------------------------------

# Heavy quarks:
b = Particle(isSM=True, label='b', pdg=5, mass=4.7*GeV, eCharge=(-1./3.), colordim=3, spin=1./2, totalwidth=0.*GeV)
# Group b and b* in a single particle:
b = MultiParticle('b', [b, b.chargeConjugate('b')])
# (We want to be able to differentiate between t+ and t-, so we define both particles:
t = Particle(isSM=True, label='t+', pdg=6, mass=173.1*GeV, eCharge=(2./3.), colordim=3, spin=1./2, totalwidth=1.4*GeV)
tC = t.chargeConjugate()
# -----------------------------------------------------------------------------------------------------------------------------

# Gauge bosons:
g = Particle(isSM=True, label='g', pdg=21, mass=0.*MeV, eCharge=0, colordim=8, spin=1, totalwidth=0.*GeV)
photon = Particle(isSM=True, label='photon', pdg=22, mass=0.*MeV, eCharge=0, colordim=1, spin=1, totalwidth=0.*GeV, _isInvisible=False)
Z = Particle(isSM=True, label='Z', pdg=23, mass=91.*GeV, eCharge=0, colordim=1, spin=1, totalwidth=2.5*GeV, _isInvisible=False)
# We group each electrically neutral gauge boson with its antiparticle:
g = MultiParticle('g', [g, g.chargeConjugate('g')])
photon = MultiParticle('photon', [photon, photon.chargeConjugate('photon')])
Z = MultiParticle('Z', [Z, Z.chargeConjugate('Z')])
W = Particle(isSM=True, label='W+', pdg=24, mass=80.*GeV, eCharge=1, colordim=1, spin=1, totalwidth=2.0*GeV)
WC = W.chargeConjugate()
# -----------------------------------------------------------------------------------------------------------------------------

# Higgs:
higgs = Particle(isSM=True, label='higgs', pdg=25, mass=125.*GeV, eCharge=0, colordim=1, spin=0, totalwidth=0.*GeV, _isInvisible=False)
# Group higgs and conjugate in single particle:
higgs = MultiParticle('higgs', [higgs, higgs.chargeConjugate('higgs')])
# -----------------------------------------------------------------------------------------------------------------------------

#  A dummy particle to represent the primary vertex (production mode)
pv = Particle(isSM=True, label='PV', pdg=0)


# Pions
pip = Particle(isSM=True, label='pi+', pdg=211, mass=140.*MeV, eCharge=+1, colordim=1, spin=0, totalwidth=0.*GeV)
piz = Particle(isSM=True, label='pi0', pdg=111, mass=140.*MeV, eCharge=0, colordim=1, spin=0, totalwidth=0.*GeV)
# Group all pions in a single particle:
pion = MultiParticle('pion', [pip, piz, pip.chargeConjugate(), piz.chargeConjugate('pi0')])

leptons = [e, mu, ta, eC, muC, taC, nu]
gauge = [g, photon, Z, W, WC]
quarks = [q, c, b, t, tC]

SMList = leptons + gauge + quarks + [higgs, pion, pv]


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
jetbList = MultiParticle('jetb', [q, c, g, b])

# Used to construct generic z2-odd and z2-even particles:
anyBSM = Particle(label='anyBSM', isSM=False)
anySM = Particle(label='anySM', isSM=True)
anyParticle = Particle(label='*')

# Used to construct BSM final states:
MET = Particle(label='MET', isSM=False, eCharge=0, colordim=1)
HSCPp = Particle(label='HSCP+', isSM=False, eCharge=+1, colordim=1)
HSCPm = Particle(label='HSCP-', isSM=False, eCharge=-1, colordim=1)
HSCP = MultiParticle(label='HSCP', particles=[HSCPp, HSCPm])

RHadronG = Particle(label='RHadronG', isSM=False, eCharge = 0, colordim = 8)
RHadronUp = Particle(label='RHadronU+', isSM=False, eCharge = 2./3., colordim = 3)
RHadronDm = Particle(label='RHadronD-', isSM=False, eCharge = -1./3., colordim = 3)
RHadronU = MultiParticle(label='RHadronU', particles = [RHadronUp,RHadronUp.chargeConjugate()])
RHadronD = MultiParticle(label='RHadronD', particles = [RHadronDm,RHadronDm.chargeConjugate()])
RHadronQ = MultiParticle(label='RHadronQ', particles = [RHadronUp,RHadronUp.chargeConjugate(),
                                                                         RHadronDm,RHadronDm.chargeConjugate()])

gluino = Particle(label='gluino', isSM=False, eCharge = 0, colordim = 8, spin = 1./2.)
chargino = Particle(label='C1+', isSM=False, eCharge = 1, colordim = 1, spin = 1./2.)
charginoBar = Particle(label='C1-', isSM=False, eCharge = -1, colordim = 1, spin = 1./2.)
C1 = MultiParticle(label='C1', particles = [chargino,charginoBar])
Hp = Particle(label='H+', isSM=False, eCharge = 1, colordim = 1, spin = 0.)
Hm = Particle(label='H-', isSM=False, eCharge = -1, colordim = 1, spin = 0.)
Hpm = MultiParticle(label='Hpm', particles = [Hp,Hm])

Zprime = Particle(label='Zprime', isSM=False, eCharge = 0, colordim = 1, spin = 1)
Vprime = Particle(label='Vprime', isSM=False, colordim = 1, spin = 1)
H0 = Particle(label='H0', isSM=False, eCharge = 0, colordim = 1, spin = 0)

#Define list of inclusive final states:
SMfinalStates = [eList,muList,taList,lpList,lmList,lList,WList,
               tList,LpList,LmList,LList,jetList,jetbList,anySM]
#Include list of exclusive final states:
SMfinalStates +=  SMList
#Define list of BSM final states:
BSMfinalStates = [MET,HSCP,RHadronU,RHadronD,RHadronG,RHadronQ,anyBSM,
                  gluino,chargino,charginoBar,C1,Hp,Hm,Hpm,Zprime,H0]


# Avoid double counting:
for i, ptc in enumerate(SMfinalStates):
    if any((ptc is p and i != j) for j, p in enumerate(SMfinalStates)):
        SMfinalStates.remove(ptc)

for i, ptc in enumerate(BSMfinalStates):
    if any((ptc is p and i != j) for j, p in enumerate(BSMfinalStates)):
        BSMfinalStates.remove(ptc)

# Protect all final state properties:
for ptc in SMfinalStates+BSMfinalStates:
    ptc._static = True

finalStates = Model(SMparticles=SMfinalStates,
                    BSMparticles=BSMfinalStates)
