"""
.. module:: SMparticleDefinitions
   :synopsis: Defines the particles to be used.
              All particles appearing in the model as well as the SM particles
              must be defined here.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.tools.physicsUnits import MeV , GeV
from smodels.theory.particle import Particle, ParticleList




# SM particles

anything = Particle(Z2parity=None, label='*', pdg=None, mass=None, eCharge=None, colordim=None, spin=None, width=None, decays=None)


e = Particle(Z2parity='even', label='e-', pdg=11, mass=0.5*MeV, eCharge=-1, colordim=0, spin=1./2, width=0, decays=None) 
mu = Particle(Z2parity='even', label='mu-', pdg=13, mass=106.*MeV, eCharge=-1, colordim=0, spin=1./2, width=0, decays=None)
ta = Particle(Z2parity='even', label='ta-', pdg=15, mass=1777.*MeV, eCharge=-1, colordim=0, spin=1./2, width=0, decays=None)

nue = Particle(Z2parity='even', label='neu', pdg=12, mass=0.*MeV, eCharge=0, colordim=0, spin=1./2, width=0, decays=None)
numu = Particle(Z2parity='even', label='neu', pdg=14, mass=0.*MeV, eCharge=0, colordim=0, spin=1./2, width=0, decays=None)
nuta = Particle(Z2parity='even', label='neu', pdg=16, mass=0.*MeV, eCharge=0, colordim=0, spin=1./2, width=0, decays=None)

d = Particle(Z2parity='even', label='q', pdg=1, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, width=0, decays=None)
u = Particle(Z2parity='even', label='q', pdg=2, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, width=0, decays=None)
s = Particle(Z2parity='even', label='q', pdg=3, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, width=0, decays=None)
c = Particle(Z2parity='even', label='c', pdg=4, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, width=0, decays=None)
b = Particle(Z2parity='even', label='b', pdg=5, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, width=0, decays=None)
t = Particle(Z2parity='even', label='t+', pdg=6, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, width=1.4*GeV, decays=None)

g = Particle(Z2parity='even', label='g', pdg=21, mass=0.*MeV, eCharge=0, colordim=8, spin=1, width=0, decays=None)
photon = Particle(Z2parity='even', label='photon',pdg=22, mass=0.*MeV, eCharge=0, colordim=0, spin=1, width=0, decays=None)
Z = Particle(Z2parity='even', label='Z', pdg=23, mass=91.*GeV , eCharge=0, colordim=0, spin=1, width=2.5*GeV, decays=None)
W = Particle(Z2parity='even', label='W+', pdg=24, mass=80.*GeV, eCharge=1, colordim=0, spin=1, width=2.0*GeV, decays=None)
higgs = Particle(Z2parity='even', label='higgs', pdg=25, mass=125.*GeV, eCharge=0, colordim=0, spin=0, width=0, decays=None)

pi = Particle(Z2parity='even', label='pi', pdg=211, mass=140.*MeV, eCharge=+1, colordim=0, spin=0, width=0, decays=None)


quarks = [u,d] + [c,s] + [t,b]
quarksC = [p.chargeConjugate() for p in quarks]
leptons = [e,nue] + [mu,numu] + [ta,nuta]
leptonsC = [p.chargeConjugate() for p in leptons]
gauge = [g,photon,W,Z]
gaugeC = [p.chargeConjugate() for p in gauge]

SMparticles = quarks + leptons + gauge + [higgs] + [pi] + [anything]
SMparticlesC = quarksC + leptonsC + gaugeC + [higgs.chargeConjugate()] + [pi.chargeConjugate()]

SMList = SMparticles + SMparticlesC

SMparticleList = ParticleList( 'SM', SMList)
SMpdgs = SMparticleList.getPdgs()
SMLabels = SMparticleList.getLabels()





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
jetList = ParticleList('jet' ,  quarks[0:4] + [ gauge[0] ] + [ pi ])
jetbarList =  ParticleList('jet~' , quarksC[0:4] + [ gaugeC[0] ] + [ pi.chargeConjugate() ])
allParticles = ParticleList('all' , SMList)  
          
particleLists = [eList, muList,taList,lpList,lmList,lList,nuList,WList,tList,
LpList,LmList,LList,jetList,jetbarList,allParticles]

