"""
.. module:: SMparticleDefinitions
   :synopsis: Defines the SM particles.
   
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.tools.physicsUnits import MeV , GeV
from smodels.theory.particle import Particle, MultiParticle

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
SMparticleList = MultiParticle('SM',SMList)
