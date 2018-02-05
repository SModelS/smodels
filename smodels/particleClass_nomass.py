from smodels.tools.physicsUnits import MeV, GeV

class Particles(object):
	
	""" A SM particle with the following properties:
		
		label: string
		pdg: number in pdg
		mass: in MeV
		charge: as multiples of the unit charge
		color charge: 
		width: in GeV  """

	SMList = []
	BSMList = []

	def __init__(self, label, Z2parity, pdg, mass, charge, colordim, width, branches):
		
		self.Z2parity = Z2parity
		self.label = label
		self.pdg = pdg
		self.mass = mass
		self.charge = charge
		self.colordim = colordim
		self.width = width
		self.branches = branches
		
		if Z2parity=='odd': Particles.BSMList.append(self)
		elif Z2parity=='even': Particles.SMList.append(self)
		else: print("Z2parity must be 'even' or 'odd' ")


	def __eq__(self, other): 
		return self.__dict__ == other.__dict__
	
	def __str__(self): 
		return str(self.__dict__ )




# SM particles

e_min = Particles('e-', 'even',11, 0.5*MeV, -1, 0, 0, None) 
mu_min = Particles('mu-', 'even',13, 106.*MeV, -1, 0, 0, None)
ta_min = Particles('ta-', 'even',15, 1777.*MeV, -1, 0, 0, None)
enu = Particles('enu', 'even',12, 0.*MeV, 0, 0, 0, None)
munu = Particles('munu', 'even',14, 0.*MeV, 0, 0, 0, None)
tanu = Particles('tanu', 'even',16, 0.*MeV, 0, 0, 0, None)

e_plus = Particles('e+', 'even',-11, 0.5*MeV, 1, 0, 0, None) 
mu_plus = Particles('mu+', 'even',-13, 106.*MeV, 1, 0, 0, None)
ta_plus = Particles('ta+', 'even',-15, 1777.*MeV, 1, 0, 0, None)


d = Particles('d', 'even',1, 0.*MeV, (-1./3.), 3, 0, None)
u = Particles('u', 'even',2, 0.*MeV, (2./3.), 3, 0, None)
s = Particles('s', 'even',3, 0.*MeV, (-1./3.), 3, 0, None)
c = Particles('c', 'even',4, 0.*MeV, (2./3.), 3, 0, None)
b = Particles('b', 'even',5, 0.*MeV, (-1./3.), 3, 0, None)
t = Particles('t', 'even',6, 0.*MeV, (2./3.), 3, 1.4*GeV, None)

dbar = Particles('dbar', 'even',-1, 0.*MeV, (1./3.), 3, 0, None)
ubar = Particles('ubar', 'even',-2, 0.*MeV, (-2./3.), 3, 0, None)
sbar = Particles('sbar', 'even',-3, 0.*MeV, (1./3.), 3, 0, None)
cbar = Particles('cbar', 'even',-4, 0.*MeV, (-2./3.), 3, 0, None)
bbar = Particles('bbar', 'even',-5, 0.*MeV, (1./3.), 3, 0, None)
tbar = Particles('tbar', 'even',-6, 0.*MeV, (-2./3.), 3, 1.4*GeV, None)

g = Particles('g', 'even',21, 0.*MeV, 0, 8, 0, None)
photon = Particles('photon', 'even',22, 0.*MeV, 0, 0, 0, None)
Z = Particles('Z', 'even',23, 91.*GeV , 0, 0, 2.5*GeV, None)
W_plus = Particles('W+', 'even',24, 80.*GeV, 1, 0, 2.0*GeV, None)
W_min = Particles('W-', 'even',-24, 80.*GeV, -1, 0, 2.0*GeV, None)
higgs = Particles('higgs', 'even',25, 125.*GeV, 0, 0, 0, None)

SMList = Particles.SMList
SMpdgs = [particle.pdg for particle in SMList]
SMparticles = [particle.label for particle in SMList]



# BSM particles

snu = Particles('sneutrino', 'odd', 1000012, None, 0, 0, None, None)
gluino = Particles('gluino', 'odd', 1000021, None, 0, 8, None, None)
N1 = Particles('N1', 'odd', 1000022, None, 0, 0, None, None)  


BSMList = Particles.BSMList
BSMpdgs = [particle.pdg for particle in BSMList]
BSMparticles = [particle.label for particle in BSMList]



#Particle groups 

ptcDic = {"e"  : [e_min, e_plus],
          "mu" : [mu_min, mu_plus],
          "ta" : [ta_min, ta_plus],
          "l+" : [e_plus, mu_plus],
          "l-" : [e_min,  mu_min, e_plus, mu_plus],
		  "l"  : [e_min,  mu_min], 
		  "nu" : [enu, munu, tanu],
          "W"  : [W_plus, W_min],
          "t+" : [t],
		  "t-" : [tbar],
		  "t"  : [t,tbar],
          "L+" : [e_plus, mu_plus, ta_plus],
          "L-" : [e_min,  mu_min, ta_min],
		  "L"  : [e_plus, mu_plus, ta_plus, e_min,  mu_min, ta_min],
          "jet" : [ d, u, s, g, c ],
		  "jetbar" : [ dbar, ubar, sbar, g, cbar ],
          "all" : [e_min,  mu_min, ta_min, e_plus, mu_plus, ta_plus, W_plus, W_min, Z, photon, higgs, d, u, s, t, b, c, g]}



