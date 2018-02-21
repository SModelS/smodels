"""
.. module:: particles
   :synopsis: Defines the particles to be used.
              All particles appearing in the model as well as the SM particles
              must be defined here.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.tools.physicsUnits import MeV , GeV
from smodels.theory.particleClass import Particles, ParticleList




# SM particles


e = Particles(Z2parity='even', label='e-', pdg=11, mass=0.5*MeV, eCharge=-1, colordim=0, spin=1./2, width=0, branches=None) 
mu = Particles(Z2parity='even', label='mu-', pdg=13, mass=106.*MeV, eCharge=-1, colordim=0, spin=1./2, width=0, branches=None)
ta = Particles(Z2parity='even', label='ta-', pdg=15, mass=1777.*MeV, eCharge=-1, colordim=0, spin=1./2, width=0, branches=None)

enu = Particles(Z2parity='even', label='enu', pdg=12, mass=0.*MeV, eCharge=0, colordim=0, spin=1./2, width=0, branches=None)
munu = Particles(Z2parity='even', label='munu', pdg=14, mass=0.*MeV, eCharge=0, colordim=0, spin=1./2, width=0, branches=None)
tanu = Particles(Z2parity='even', label='tanu', pdg=16, mass=0.*MeV, eCharge=0, colordim=0, spin=1./2, width=0, branches=None)


d = Particles(Z2parity='even', label='d', pdg=1, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, width=0, branches=None)
u = Particles(Z2parity='even', label='u', pdg=2, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, width=0, branches=None)
s = Particles(Z2parity='even', label='s', pdg=3, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, width=0, branches=None)
c = Particles(Z2parity='even', label='c', pdg=4, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, width=0, branches=None)
b = Particles(Z2parity='even', label='b', pdg=5, mass=0.*MeV, eCharge=(-1./3.), colordim=3, spin=1./2, width=0, branches=None)
t = Particles(Z2parity='even', label='t+', pdg=6, mass=0.*MeV, eCharge=(2./3.), colordim=3, spin=1./2, width=1.4*GeV, branches=None)

g = Particles(Z2parity='even', label='g', pdg=21, mass=0.*MeV, eCharge=0, colordim=8, spin=1, width=0, branches=None)
photon = Particles(Z2parity='even', label='photon',pdg=22, mass=0.*MeV, eCharge=0, colordim=0, spin=1, width=0, branches=None)
Z = Particles(Z2parity='even', label='Z', pdg=23, mass=91.*GeV , eCharge=0, colordim=0, spin=1, width=2.5*GeV, branches=None)
W = Particles(Z2parity='even', label='W+', pdg=24, mass=80.*GeV, eCharge=1, colordim=0, spin=1, width=2.0*GeV, branches=None)
higgs = Particles(Z2parity='even', label='higgs', pdg=25, mass=125.*GeV, eCharge=0, colordim=0, spin=0, width=0, branches=None)





quarks = [u,d] + [c,s] + [t,b]
quarksC = [p.chargeConjugate() for p in quarks]
leptons = [e,enu] + [mu,munu] + [ta,tanu]
leptonsC = [p.chargeConjugate() for p in leptons]
gauge = [g,photon,W,Z]
gaugeC = [p.chargeConjugate() for p in gauge]

SMparticles = quarks + leptons + gauge + [higgs]
SMparticlesC = quarksC + leptonsC + gaugeC + [higgs.chargeConjugate()]

SMList = SMparticles + SMparticlesC

SMparticleList = ParticleList( 'SM', SMList)
SMpdgs = SMparticleList.getPdgs()
SMparticles = SMparticleList.getNames()





#Particle groups 

eList = ParticleList('e' , [leptons[0], leptonsC[0]])
muList = ParticleList('mu' , [leptons[2], leptonsC[2]])
taList = ParticleList('ta' , [leptons[4], leptonsC[4]])
lpList = ParticleList('l+' , [leptonsC[0],leptonsC[2]])
lmList = ParticleList('l-'  , [leptons[0],leptons[2]]) 
lList = ParticleList('l' , [lpList,lmList])
nuList = ParticleList('nu' , [leptons[1],leptons[3],leptons[5],
                                 leptonsC[1],leptonsC[3],leptonsC[5]])
WList = ParticleList('W'  , [gauge[2],gaugeC[2]])
tList = ParticleList('t'  , [quarks[4],quarksC[4]])
LpList = ParticleList('L+' , [lpList,leptonsC[4]])
LmList = ParticleList('L-' , [lmList,leptons[4]])
LList = ParticleList('L'  , [LpList,LmList])
jet = ParticleList('jet' ,  quarks[0:4] + [ gauge[0] ])
jetbar =  ParticleList('jetbar' , quarksC[0:4] + [ gaugeC[0] ])
allParticles = ParticleList('all' , SMList)  
          
particleLists = [eList, muList,taList,lpList,lmList,lList,nuList,WList,tList,
LpList,LmList,LList,jet,jetbar,allParticles]



