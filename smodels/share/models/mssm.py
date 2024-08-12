"""
.. module:: mssm
   :synopsis: Defines the BSM particles to be used. Properties not defined here and defined by the LHE or SLHA input file (such as masses, width and BRs) are automatically added later.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.base.particle import Particle

####  R-odd   ##########
#1st generation squarks and its conjugates:
sdl = Particle(isSM=False, label='sd_L', pdg=1000001, eCharge=-1./3, colordim=3, spin=0)
sul = Particle(isSM=False, label='su_L', pdg=1000002, eCharge=2./3, colordim=3, spin=0)
sdr = Particle(isSM=False, label='sd_R', pdg=2000001, eCharge=-1./3, colordim=3, spin=0)
sur = Particle(isSM=False, label='su_R', pdg=2000002, eCharge=2./3, colordim=3, spin=0)

#2nd generation squarks and its conjugates:
ssl = Particle(isSM=False, label='ss_L', pdg=1000003, eCharge=-1./3, colordim=3, spin=0)
scl = Particle(isSM=False, label='sc_L', pdg=1000004, eCharge=2./3, colordim=3, spin=0)
ssr = Particle(isSM=False, label='ss_R', pdg=2000003, eCharge=-1./3, colordim=3, spin=0)
scr = Particle(isSM=False, label='sc_R', pdg=2000004, eCharge=2./3, colordim=3, spin=0)

#3rd generation squarks and its conjugates:
sb1 = Particle(isSM=False, label='sb_1', pdg=1000005, eCharge=-1./3, colordim=3, spin=0)
st1 = Particle(isSM=False, label='st_1', pdg=1000006, eCharge=2./3, colordim=3, spin=0)
sb2 = Particle(isSM=False, label='sb_2', pdg=2000005, eCharge=-1./3, colordim=3, spin=0)
st2 = Particle(isSM=False, label='st_2', pdg=2000006, eCharge=2./3, colordim=3, spin=0)

#1st generation sleptons and its conjugates:
sel = Particle(isSM=False, label='se_L', pdg=1000011, eCharge=-1, colordim=1, spin=0)
snel = Particle(isSM=False, label='sne_L', pdg=1000012, eCharge=0, colordim=1, spin=0)
ser = Particle(isSM=False, label='se_R', pdg=2000011, eCharge=-1, colordim=1, spin=0)

#2nd generation sleptons and its conjugates:
smul = Particle(isSM=False, label='smu_L', pdg=1000013, eCharge=-1, colordim=1, spin=0)
snmul = Particle(isSM=False, label='snmu_L', pdg=1000014, eCharge=0, colordim=1, spin=0)
smur = Particle(isSM=False, label='smu_R', pdg=2000013, eCharge=-1, colordim=1, spin=0)

#3rd generation sleptons and its conjugates:
sta1 = Particle(isSM=False, label='sta_1', pdg=1000015, eCharge=-1, colordim=1, spin=0)
sntal = Particle(isSM=False, label='snta_L', pdg=1000016, eCharge=0, colordim=1, spin=0)
sta2 = Particle(isSM=False, label='sta_2', pdg=2000015, eCharge=-1, colordim=1, spin=0)

#Gluino:
gluino = Particle(isSM=False, label='gluino', pdg=1000021, eCharge=0, colordim=8, spin=1./2)
#Neutralinos
n1 = Particle(isSM=False, label='N1', pdg=1000022, eCharge=0, colordim=1, spin=1./2)
n2 = Particle(isSM=False, label='N2', pdg=1000023, eCharge=0, colordim=1, spin=1./2)
n3 = Particle(isSM=False, label='N3', pdg=1000025, eCharge=0, colordim=1, spin=1./2)
n4 = Particle(isSM=False, label='N4', pdg=1000035, eCharge=0, colordim=1, spin=1./2)

#Gravitino:
# g = Particle(isSM=False, label='G', pdg=1000039, eCharge=0, colordim=1, spin=3./2)

#Charginos
c1 = Particle(isSM=False, label='C1+', pdg=1000024, eCharge=1, colordim=1, spin=1./2)
c2 = Particle(isSM=False, label='C2+', pdg=1000037, eCharge=1, colordim=1, spin=1./2)

##### R-even  ###############
#Higgs
H = Particle(isSM=False, label='H+', pdg=37, eCharge=+1, colordim=1, spin=0)
A0 = Particle(isSM=False, label='A0', pdg=36, eCharge=0, colordim=1, spin=0, _isInvisible=False)
H0 = Particle(isSM=False, label='H0', pdg=35, eCharge=0, colordim=1, spin=0, _isInvisible=False)


squarks = [sdl,sul,sdr,sur] + [ssl,scl,ssr,scr] + [sb1,st1,sb2,st2]
sleptons = [sel,snel,ser] + [smul,snmul,smur] + [sta1,sntal,sta2]
inos = [gluino] + [n1,n2,n3,n4] + [c1,c2] # + [g]

rOdd = squarks + sleptons + inos
rOddC = [p.chargeConjugate() for p in rOdd]  #Define the charge conjugates

bsm_higgs = [H,A0,H0]
bsm_higgsC = [p.chargeConjugate() for p in bsm_higgs]

#Generic BSM particles:

BSMList = rOdd + rOddC + bsm_higgs + bsm_higgsC
