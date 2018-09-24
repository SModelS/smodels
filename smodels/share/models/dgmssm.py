"""
.. module:: dgmssm 
   :synopsis: Defines the BSM particles to be used
   .
.. moduleauthor:: Sabine Kraml <sabine.kraml@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>

   Properties not defined here and defined by the LHE or SLHA input file 
   (such as masses, width and BRs) are automatically added later.
"""

from smodels.theory.particle import Particle, ParticleList

####  R-odd   ##########
#1st generation squarks and its conjugates:
sdl = Particle(Z2parity='odd', label='sd_L', pdg=1000001, eCharge=-1./3, colordim=3, spin=0)
sul = Particle(Z2parity='odd', label='su_L', pdg=1000002, eCharge=2./3, colordim=3, spin=0)
sdr = Particle(Z2parity='odd', label='sd_R', pdg=2000001, eCharge=-1./3, colordim=3, spin=0)
sur = Particle(Z2parity='odd', label='su_R', pdg=2000002, eCharge=2./3, colordim=3, spin=0)

#2nd generation squarks and its conjugates:
ssl = Particle(Z2parity='odd', label='ss_L', pdg=1000003, eCharge=-1./3, colordim=3, spin=0)
scl = Particle(Z2parity='odd', label='sc_L', pdg=1000004, eCharge=2./3, colordim=3, spin=0)
ssr = Particle(Z2parity='odd', label='ss_R', pdg=2000003, eCharge=-1./3, colordim=3, spin=0)
scr = Particle(Z2parity='odd', label='sc_R', pdg=2000004, eCharge=2./3, colordim=3, spin=0)

#3rd generation squarks and its conjugates:
sb1 = Particle(Z2parity='odd', label='sb_1', pdg=1000005, eCharge=2./3, colordim=3, spin=0)
st1 = Particle(Z2parity='odd', label='st_1', pdg=1000006, eCharge=-1./3, colordim=3, spin=0)
sb2 = Particle(Z2parity='odd', label='sb_2', pdg=2000005, eCharge=2./3, colordim=3, spin=0)
st2 = Particle(Z2parity='odd', label='st_2', pdg=2000006, eCharge=-1./3, colordim=3, spin=0)

#1st generation sleptons and its conjugates:
sel = Particle(Z2parity='odd', label='se_L', pdg=1000011, eCharge=-1, colordim=1, spin=0)
snel = Particle(Z2parity='odd', label='sne_L', pdg=1000012, eCharge=0, colordim=1, spin=0)
ser = Particle(Z2parity='odd', label='se_R', pdg=2000011, eCharge=0, colordim=1, spin=0)
sner = Particle(Z2parity='odd', label='sne_R', pdg=2000012, eCharge=0, colordim=1, spin=0)

#2nd generation sleptons and its conjugates:
smul = Particle(Z2parity='odd', label='smu_L', pdg=1000013, eCharge=-1, colordim=1, spin=0)
snmul = Particle(Z2parity='odd', label='snmu_L', pdg=1000014, eCharge=0, colordim=1, spin=0)
smur = Particle(Z2parity='odd', label='smu_R', pdg=2000013, eCharge=-1, colordim=1, spin=0)
snmur = Particle(Z2parity='odd', label='snmu_R', pdg=2000014, eCharge=0, colordim=1, spin=0)

#3rd generation sleptons and its conjugates:
sta1 = Particle(Z2parity='odd', label='sta_1', pdg=1000015, eCharge=-1, colordim=1, spin=0)
sntal = Particle(Z2parity='odd', label='snta_L', pdg=1000016, eCharge=0, colordim=1, spin=0)
sta2 = Particle(Z2parity='odd', label='sta_2', pdg=2000015, eCharge=-1, colordim=1, spin=0)
sntar = Particle(Z2parity='odd', label='snta_R', pdg=2000016, eCharge=0, colordim=1, spin=0)

#Gluino:
gluino1 = Particle(Z2parity='odd', label='gluino1', pdg=1000021, eCharge=0, colordim=8, spin=1./2)
gluino2 = Particle(Z2parity='odd', label='gluino2', pdg=2000021, eCharge=0, colordim=8, spin=1./2)

#Neutralinos
n1 = Particle(Z2parity='odd', label='N1', pdg=1000022, eCharge=0, colordim=1, spin=1./2)  
n2 = Particle(Z2parity='odd', label='N2', pdg=1000023, eCharge=0, colordim=1, spin=1./2)  
n3 = Particle(Z2parity='odd', label='N3', pdg=1000025, eCharge=0, colordim=1, spin=1./2)  
n4 = Particle(Z2parity='odd', label='N4', pdg=1000035, eCharge=0, colordim=1, spin=1./2)  
n5 = Particle(Z2parity='odd', label='N5', pdg=1000045, eCharge=0, colordim=1, spin=1./2)  
n6 = Particle(Z2parity='odd', label='N6', pdg=1000055, eCharge=0, colordim=1, spin=1./2)  

#Charginos
c1 = Particle(Z2parity='odd', label='C1+', pdg=1000024, eCharge=1, colordim=1, spin=1./2)  
c2 = Particle(Z2parity='odd', label='C2+', pdg=1000037, eCharge=1, colordim=1, spin=1./2)  
c3 = Particle(Z2parity='odd', label='C3+', pdg=1000047, eCharge=1, colordim=1, spin=1./2)  

#Gravitino
gravitino = Particle(Z2parity='odd', label='G', pdg=1000039, eCharge=1, colordim=1, spin=1./2)

##### R-even  ###############
#Higgs
H1p = Particle(Z2parity='even', label='H1+', pdg=37, eCharge=+1, colordim=1, spin=0)    
H2p = Particle(Z2parity='even', label='H2+', pdg=47, eCharge=+1, colordim=1, spin=0)    
H3p = Particle(Z2parity='even', label='H3+', pdg=57, eCharge=+1, colordim=1, spin=0)    
H2 = Particle(Z2parity='even', label='H2', pdg=35, eCharge=0, colordim=1, spin=0)  
H3 = Particle(Z2parity='even', label='H3', pdg=45, eCharge=0, colordim=1, spin=0)  
H4 = Particle(Z2parity='even', label='H4', pdg=55, eCharge=0, colordim=1, spin=0)  
A1 = Particle(Z2parity='even', label='A1', pdg=36, eCharge=0, colordim=1, spin=0)
A2 = Particle(Z2parity='even', label='A2', pdg=46, eCharge=0, colordim=1, spin=0)
A3 = Particle(Z2parity='even', label='A3', pdg=56, eCharge=0, colordim=1, spin=0)

#Sgluons
sgluon1 = Particle(Z2parity='odd', label='sgluon1', pdg=3000021, eCharge=0, colordim=8, spin=1./2)
sgluon1 = Particle(Z2parity='odd', label='sgluon2', pdg=3000022, eCharge=0, colordim=8, spin=1./2)


squarks = [sdl,sul,sdr,sur] + [ssl,scl,ssr,scr] + [sb1,st1,sb2,st2]
sleptons = [sel,snel,ser, snelr] + [smul,snmul,smur,snmur] + [sta1,sntal,sta2,sntar]
inos = [gluino1,gluino2] + [n1,n2,n3,n4,n5,n6] + [c1,c2,c3] + [gravitino]

rOdd = squarks + sleptons + inos 
rOddC = [p.chargeConjugate() for p in rOdd]  #Define the charge conjugates

higgs = [H1p,H2p,H3p,H2,H3,H4,A1,A2,A3]
higgsC = [p.chargeConjugate() for p in higgs]
from smodels.share.models.SMparticles import SMList
rEven = SMList + higgs + higgsC

#Generic BSM particles:

BSMList = rOdd + rOddC
BSMparticleList = ParticleList('BSM', BSMList)

