"""
.. module:: SMparticleDefinitions
   :synopsis: Defines the SM particles.
   
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.tools.physicsUnits import MeV , GeV
from smodels.theory.particle import Particle, MultiParticle




# SM particles
e = Particle(Z2parity='even', label='e-', pdg=11, mass=0.5*MeV, eCharge=-1, colordim=0, spin=1./2, totalwidth = 0.*GeV, decays=[]) 
mu = Particle(Z2parity='even', label='mu-', pdg=13, mass=106.*MeV, eCharge=-1, colordim=0, spin=1./2, totalwidth = 0.*GeV, decays=[])
ta = Particle(Z2parity='even', label='ta-', pdg=15, mass=1777.*MeV, eCharge=-1, colordim=0, spin=1./2, totalwidth = 0.*GeV, decays=[])

nue = Particle(Z2parity='even', label='nu', pdg=12, mass=0.*MeV, eCharge=0, colordim=0, spin=1./2, totalwidth = 0.*GeV, decays=[], _isMET=True)
numu = Particle(Z2parity='even', label='nu', pdg=14, mass=0.*MeV, eCharge=0, colordim=0, spin=1./2, totalwidth = 0.*GeV, decays=[], _isMET=True)
nuta = Particle(Z2parity='even', label='nu', pdg=16, mass=0.*MeV, eCharge=0, colordim=0, spin=1./2, totalwidth = 0.*GeV, decays=[], _isMET=True)

d = Particle(Z2parity='even', label='q', pdg=1, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, totalwidth = 0.*GeV, decays=[])
u = Particle(Z2parity='even', label='q', pdg=2, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, totalwidth = 0.*GeV, decays=[])
s = Particle(Z2parity='even', label='q', pdg=3, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, totalwidth = 0.*GeV, decays=[])
c = Particle(Z2parity='even', label='c', pdg=4, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, totalwidth = 0.*GeV, decays=[])
b = Particle(Z2parity='even', label='b', pdg=5, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, totalwidth = 0.*GeV, decays=[])
t = Particle(Z2parity='even', label='t+', pdg=6, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, totalwidth=1.4*GeV, decays=[])

g = Particle(Z2parity='even', label='g', pdg=21, mass=0.*MeV, eCharge=0, colordim=8, spin=1, totalwidth = 0.*GeV, decays=[])
photon = Particle(Z2parity='even', label='photon',pdg=22, mass=0.*MeV, eCharge=0, colordim=0, spin=1, totalwidth = 0.*GeV, decays=[], _isMET=False)
Z = Particle(Z2parity='even', label='Z', pdg=23, mass=91.*GeV , eCharge=0, colordim=0, spin=1, totalwidth=2.5*GeV, decays=[], _isMET=False)
W = Particle(Z2parity='even', label='W+', pdg=24, mass=80.*GeV, eCharge=1, colordim=0, spin=1, totalwidth=2.0*GeV, decays=[])
higgs = Particle(Z2parity='even', label='higgs', pdg=25, mass=125.*GeV, eCharge=0, colordim=0, spin=0, totalwidth = 0.*GeV, decays=[], _isMET=False)

pip = Particle(Z2parity='even', label='pi', pdg=211, mass=140.*MeV, eCharge=+1, colordim=0, spin=0, totalwidth = 0.*GeV, decays=[])
piz = Particle(Z2parity='even', label='pi', pdg=111, mass=140.*MeV, eCharge=+1, colordim=0, spin=0, totalwidth = 0.*GeV, decays=[])


quarks = [u,d] + [c,s] + [t,b]
quarksC = [p.chargeConjugate(p.label) if p.label != 't+' else p.chargeConjugate() for p in quarks]
chargedLeptons = [e,mu,ta]
neutrinos = [nue,numu,nuta]
leptons = chargedLeptons + neutrinos
leptonsC = [p.chargeConjugate() for p in chargedLeptons] + [p.chargeConjugate('nu') for p in neutrinos] 
gauge = [g,photon,W,Z]
gaugeC = [p.chargeConjugate() for p in gauge]

SMparticles = quarks + leptons + gauge + [higgs] + [pip,piz]
SMparticlesC = quarksC + leptonsC + gaugeC + [higgs.chargeConjugate()] + [pip.chargeConjugate('pi'),piz.chargeConjugate('pi')]

SMList = SMparticles + SMparticlesC
#Protect all particles properties:
for ptc in SMList:
    ptc._static = True

SMparticleList = MultiParticle('SM', SMList)
          
