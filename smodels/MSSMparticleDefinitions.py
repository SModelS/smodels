"""
.. module:: particles
   :synopsis: Defines the particles to be used.
              All particles appearing in the model as well as the SM particles
              must be defined here.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>


   Properties not defined here and defined by the LHE or SLHA input file 
   (such as masses, width and BRs) are automatically added later.
"""

from smodels.tools.physicsUnits import MeV, GeV
from smodels.theory.particle import Particle, ParticleList


# MSSM particles

####  R-odd   ##########
#1st generation squarks and its conjugates:
sdl = Particle(Z2parity='odd', label='sd_L', pdg=1000001, mass=None, eCharge=-1./3, colordim=3, spin=0, width=None, decays=None)
sul = Particle(Z2parity='odd', label='su_L', pdg=1000002, mass=None, eCharge=2./3, colordim=3, spin=0, width=None, decays=None)
sdr = Particle(Z2parity='odd', label='sd_R', pdg=2000001, mass=None, eCharge=-1./3, colordim=3, spin=0, width=None, decays=None)
sur = Particle(Z2parity='odd', label='su_R', pdg=2000002, mass=None, eCharge=2./3, colordim=3, spin=0, width=None, decays=None)

#2nd generation squarks and its conjugates:
ssl = Particle(Z2parity='odd', label='ss_L', pdg=1000003, mass=None, eCharge=-1./3, colordim=3, spin=0, width=None, decays=None)
scl = Particle(Z2parity='odd', label='sc_L', pdg=1000004, mass=None, eCharge=2./3, colordim=3, spin=0, width=None, decays=None)
ssr = Particle(Z2parity='odd', label='ss_R', pdg=2000003, mass=None, eCharge=-1./3, colordim=3, spin=0, width=None, decays=None)
scr = Particle(Z2parity='odd', label='sc_R', pdg=2000004, mass=None, eCharge=2./3, colordim=3, spin=0, width=None, decays=None)

#3rd generation squarks and its conjugates:
sb1 = Particle(Z2parity='odd', label='sb_1', pdg=1000005, mass=None, eCharge=2./3, colordim=3, spin=0, width=None, decays=None)
st1 = Particle(Z2parity='odd', label='st_1', pdg=1000006, mass=None, eCharge=-1./3, colordim=3, spin=0, width=None, decays=None)
sb2 = Particle(Z2parity='odd', label='sb_2', pdg=2000005, mass=None, eCharge=2./3, colordim=3, spin=0, width=None, decays=None)
st2 = Particle(Z2parity='odd', label='st_2', pdg=2000006, mass=None, eCharge=-1./3, colordim=3, spin=0, width=None, decays=None)

#1st generation sleptons and its conjugates:
sel = Particle(Z2parity='odd', label='se_L', pdg=1000011, mass=None, eCharge=-1, colordim=0, spin=0, width=None, decays=None)
snel = Particle(Z2parity='odd', label='sne_L', pdg=1000012, mass=None, eCharge=0, colordim=0, spin=0, width=None, decays=None)
ser = Particle(Z2parity='odd', label='se_R', pdg=2000011, mass=None, eCharge=0, colordim=0, spin=0, width=None, decays=None)

#2nd generation sleptons and its conjugates:
smul = Particle(Z2parity='odd', label='smu_L', pdg=1000013, mass=None, eCharge=-1, colordim=0, spin=0, width=None, decays=None)
snmul = Particle(Z2parity='odd', label='snmu_L', pdg=1000014, mass=None, eCharge=0, colordim=0, spin=0, width=None, decays=None)
smur = Particle(Z2parity='odd', label='smu_R', pdg=2000013, mass=None, eCharge=-1, colordim=0, spin=0, width=None, decays=None)

#3rd generation sleptons and its conjugates:
sta1 = Particle(Z2parity='odd', label='sta_1', pdg=1000015, mass=None, eCharge=-1, colordim=0, spin=0, width=None, decays=None)
sntal = Particle(Z2parity='odd', label='snta_L', pdg=1000016, mass=None, eCharge=0, colordim=0, spin=0, width=None, decays=None)
sta2 = Particle(Z2parity='odd', label='sta_2', pdg=2000015, mass=None, eCharge=-1, colordim=0, spin=0, width=None, decays=None)

#Gluino:
gluino = Particle(Z2parity='odd', label='gluino', pdg=1000021, mass=None, eCharge=0, colordim=8, spin=1./2, width=None, decays=None)
#Neutralinos
n1 = Particle(Z2parity='odd', label='N1', pdg=1000022, mass=None, eCharge=0, colordim=0, spin=1./2, width=None, decays=None)  
n2 = Particle(Z2parity='odd', label='N2', pdg=1000023, mass=None, eCharge=0, colordim=0, spin=1./2, width=None, decays=None)  
n3 = Particle(Z2parity='odd', label='N3', pdg=1000025, mass=None, eCharge=0, colordim=0, spin=1./2, width=None, decays=None)  
n4 = Particle(Z2parity='odd', label='N4', pdg=1000035, mass=None, eCharge=0, colordim=0, spin=1./2, width=None, decays=None)  

#Charginos
c1 = Particle(Z2parity='odd', label='C1+', pdg=1000024, mass=None, eCharge=1, colordim=0, spin=1./2, width=None, decays=None)  
c2 = Particle(Z2parity='odd', label='C2+', pdg=1000037, mass=None, eCharge=1, colordim=0, spin=1./2, width=None, decays=None)  

##### R-even  ###############
#Higgs
H = Particle(Z2parity='even', label='H+', pdg=37, mass=None, eCharge=+1, colordim=0, spin=0, width=None, decays=None)  
A0 = Particle(Z2parity='even', label='A0', pdg=36, mass=None, eCharge=0, colordim=0, spin=0, width=None, decays=None)  
H0 = Particle(Z2parity='even', label='H0', pdg=35, mass=None, eCharge=0, colordim=0, spin=0, width=None, decays=None)  


squarks = [sdl,sul,sdr,sur] + [ssl,scl,ssr,scr] + [sb1,st1,sb2,st2]
sleptons = [sel,snel,ser] + [smul,snmul,smur] + [sta1,sntal,sta2]
inos = [gluino] + [n1,n2,n3,n4] + [c1,c2]
higgs = [H,A0,H0]

sparticles = squarks + sleptons + inos + higgs
sparticlesC = [p.chargeConjugate() for p in sparticles]  #Define the charge conjugates



BSMList = sparticles + sparticlesC

BSMparticleList = ParticleList('BSM', BSMList)
BSMpdgs = BSMparticleList.getPdgs()
BSMLabels = BSMparticleList.getLabels()



