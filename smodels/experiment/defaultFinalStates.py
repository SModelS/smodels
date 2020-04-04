"""
.. module:: databaseParticles
   :synopsis: Defines the final state particles used in the experimental results.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory.particle import Particle,MultiParticle
from smodels.theory.model import Model
from smodels.tools.physicsUnits import MeV , GeV
from smodels.theory.exceptions import SModelSTheoryError as SModelSError

###Below we define all particles that are relevant for the database constraints.
###One must avoid defining non-unique labels, so there is a one-to-one correspondence
###between the particles defined in finalStates and the labels used to describe the
###particles in the simplified model.

# SM particles

##Charged leptons:
e = Particle(Z2parity=1, label='e-', pdg=11, mass=0.5*MeV, eCharge=-1, colordim=1, spin=1./2, totalwidth = 0.*GeV)
mu = Particle(Z2parity=1, label='mu-', pdg=13, mass=106.*MeV, eCharge=-1, colordim=1, spin=1./2, totalwidth = 0.*GeV)
ta = Particle(Z2parity=1, label='ta-', pdg=15, mass=1777.*MeV, eCharge=-1, colordim=1, spin=1./2, totalwidth = 0.*GeV)
eC = e.chargeConjugate()
muC = mu.chargeConjugate()
taC = ta.chargeConjugate()
##-----------------------------------------------------------------------------------------------------------------------------
##Neutrinos
nue = Particle(Z2parity=1, label='nue', pdg=12, mass=0.*MeV, eCharge=0, colordim=1, spin=1./2, totalwidth = 0.*GeV)
numu = Particle(Z2parity=1, label='numu', pdg=14, mass=0.*MeV, eCharge=0, colordim=1, spin=1./2, totalwidth = 0.*GeV)
nuta = Particle(Z2parity=1, label='nuta', pdg=16, mass=0.*MeV, eCharge=0, colordim=1, spin=1./2, totalwidth = 0.*GeV)
###Group all neutrinos into a single particle:
nu = MultiParticle('nu',[nue,numu,nuta,nue.chargeConjugate(),numu.chargeConjugate(),nuta.chargeConjugate()])
##-----------------------------------------------------------------------------------------------------------------------------

##Light quarks:
d = Particle(Z2parity=1, label='d', pdg=1, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, totalwidth = 0.*GeV)
u = Particle(Z2parity=1, label='u', pdg=2, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, totalwidth = 0.*GeV)
s = Particle(Z2parity=1, label='s', pdg=3, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, totalwidth = 0.*GeV)
#Group all light quarks in a single particle
q = MultiParticle('q', [u,d,s,u.chargeConjugate(),d.chargeConjugate(),s.chargeConjugate()])
c = Particle(Z2parity=1, label='c', pdg=4, mass=1.28*GeV, eCharge=(2./3.), colordim=3, spin=1./2, totalwidth = 0.*GeV)
###Group c and c* in a single particle:
c = MultiParticle('c',[c,c.chargeConjugate('c')])
##-----------------------------------------------------------------------------------------------------------------------------

##Heavy quarks:
b = Particle(Z2parity=1, label='b', pdg=5, mass=4.7*GeV, eCharge=(-1./3.), colordim=3, spin=1./2, totalwidth = 0.*GeV)
###Group b and b* in a single particle:
b = MultiParticle('b',[b,b.chargeConjugate('b')])
###(We want to be able to differentiate between t+ and t-, so we define both particles:
t = Particle(Z2parity=1, label='t+', pdg=6, mass=173.1*GeV, eCharge=(2./3.), colordim=3, spin=1./2, totalwidth=1.4*GeV)
tC = t.chargeConjugate()
##-----------------------------------------------------------------------------------------------------------------------------

##Gauge bosons:
g = Particle(Z2parity=1, label='g', pdg=21, mass=0.*MeV, eCharge=0, colordim=8, spin=1, totalwidth = 0.*GeV)
photon = Particle(Z2parity=1, label='photon',pdg=22, mass=0.*MeV, eCharge=0, colordim=1, spin=1, totalwidth = 0.*GeV, _isInvisible=False)
Z = Particle(Z2parity=1, label='Z', pdg=23, mass=91.*GeV , eCharge=0, colordim=1, spin=1, totalwidth=2.5*GeV, _isInvisible=False)
#We group each electrically neutral gauge boson with its antiparticle:
g = MultiParticle('g',[g,g.chargeConjugate('g')])
photon = MultiParticle('photon',[photon,photon.chargeConjugate('photon')])
Z = MultiParticle('Z',[Z,Z.chargeConjugate('Z')])
W = Particle(Z2parity=1, label='W+', pdg=24, mass=80.*GeV, eCharge=1, colordim=1, spin=1, totalwidth=2.0*GeV)
WC = W.chargeConjugate()
##-----------------------------------------------------------------------------------------------------------------------------

##Higgs:
higgs = Particle(Z2parity=1, label='higgs', pdg=25, mass=125.*GeV, eCharge=0, colordim=1, spin=0, totalwidth = 0.*GeV, _isInvisible=False)
###Group higgs and conjugate in single particle:
higgs = MultiParticle('higgs',[higgs,higgs.chargeConjugate('higgs')])
##-----------------------------------------------------------------------------------------------------------------------------

##Pions
pip = Particle(Z2parity=1, label='pi+', pdg=211, mass=140.*MeV, eCharge=+1, colordim=1, spin=0, totalwidth = 0.*GeV)
piz = Particle(Z2parity=1, label='pi0', pdg=111, mass=140.*MeV, eCharge=0, colordim=1, spin=0, totalwidth = 0.*GeV)
###Group all pions in a single particle:
pion = MultiParticle('pion',[pip,piz,pip.chargeConjugate(),piz.chargeConjugate('pi0')])

leptons = [e,mu,ta,eC,muC,taC,nu]
gauge = [g,photon,Z,W,WC]
quarks = [q,c,b,t,tC]

SMList = leptons + gauge + quarks + [higgs,pion]



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
anyOdd = Particle(label='anyOdd',Z2parity=-1)
anyEven = Particle(label='*',Z2parity=1)


#Used to construct BSM final states:
MET = Particle(label='MET', Z2parity = -1, eCharge = 0, colordim = 1)
HSCPp = Particle(label='HSCP+', Z2parity = -1, eCharge = +1, colordim = 1)
HSCPm = Particle(label='HSCP-', Z2parity = -1, eCharge = -1, colordim = 1)
HSCP = MultiParticle(label='HSCP', particles = [HSCPp,HSCPm])

RHadronG = Particle(label='RHadronG', Z2parity = -1, eCharge = 0, colordim = 8)
RHadronU = Particle(label='RHadronU', Z2parity = -1, eCharge = 2./3., colordim = 3)
RHadronD = Particle(label='RHadronD', Z2parity = -1, eCharge = -1./3., colordim = 3)
RHadronQ = MultiParticle(label='RHadronQ', particles = [RHadronU,RHadronU.chargeConjugate(),
                                                                         RHadronD,RHadronD.chargeConjugate()])

gluino = Particle(label='gluino', Z2parity = -1, eCharge = 0, colordim = 8, spin = 1./2.)
chargino = Particle(label='C1+', Z2parity = -1, eCharge = 1, colordim = 0, spin = 1./2.)


#Define list of inclusive final states:
SMfinalStates = [eList,muList,taList,lpList,lmList,lList,WList,
               tList,LpList,LmList,LList,jetList,anyEven]
#Include list of exclusive final states:
SMfinalStates +=  SMList
#Define list of BSM final states:
BSMfinalStates = [MET,HSCP,RHadronG,RHadronQ,anyOdd,gluino,chargino]

allFinalStates = SMfinalStates + BSMfinalStates
#Avoid double counting:
for i,ptc in enumerate(allFinalStates):
    if any((ptc is p and i != j) for j,p in enumerate(allFinalStates)):
        allFinalStates.remove(ptc)

#Define a dummy model just to use the facilities for filtering particles
finalStates = Model(SMparticles = allFinalStates,
                    BSMparticles = [], label='DB Final States (default)')

#Check consistency:
for label in finalStates.getValuesFor('label'):
    particles = finalStates.getParticlesWith(label=label)
    if len(particles) != 1:
        raise SModelSError("%i particles defined with label %s. Particles defined in databaseParticles must have unique labels."
                           %(len(particles),label))
