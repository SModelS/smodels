.. index:: Introduction

Introduction
============

(From the original SModelS paper, `arXiv:1312.4175 <http://arxiv.org/abs/arXiv:1312.4175>`_ )

Searches at the ATLAS and CMS experiments at the LHC show no signs of physics
beyond the Standard Model (BSM).  After the first phase of LHC operation at
centre-of-mass energies of 7–8 TeV in 2010–2012, the limits for the masses of
supersymmetric particles, in particular of 1st/2nd generation squarks and
gluinos, have been pushed well into the TeV range [1,2]. Likewise, precision
measurements in the flavor sector, in particular in B-physics, are well
consistent with Standard Model (SM) expectations [3,4] and show no sign, or
need, of new physics. At the same time the recent discovery [5,6] of a
Higgs-like particle with mass around 125 GeV makes the question of stability of
the electroweak scale—the infamous gauge hierarchy problem—even more imminent.
Indeed, supersymmetry (SUSY) is arguably the best-motivated theory to solve the
gauge hierarchy problem and to explain a light SM-like Higgs boson. So, the
Higgs has very likely been discovered—but where is supersymmetry?  

Looking closely [7–12] one soon realizes that many of the current limits on SUSY
particles are based on severe model assumptions, which impose particular
relations between particle masses, decay branching ratios, etc. The prime
example is the interpretation of the search results within the Constrained
Minimal Supersymmetric Standard Model (CMSSM). The interpretation of the search
results within a much more general realization of the MSSM is perfectly
feasible, see [7,8,13], but computationally very demanding and certainly not
suitable for a “quick” survey.  

An approach which has therefore been adopted systematically by the ATLAS and
CMS collaborations, is to interpret the results within so-called Simplified
Model Spectra [14,15]. Simplified Model Spectra, or SMS for short, are
effective-Lagrangian descriptions involving just a small number of new
particles. They were designed as a useful tool for the characterization of new
physics, see e.g. [16,17]. A large variety of results on searches in many
different channels are available from both ATLAS and CMS, providing general
cross section limits for SMS topologies. However, using these results to
constrain complex SUSY (or general BSM) scenarios is not straightforward.  In
this paper, we present a method to decompose the signal of an arbitrary SUSY
spectrum into simplified model topologies and test it against all the existing
LHC bounds in the SMS context. 
 
This document describes the computer package that does all this.
SModelS can be used just like an application, running :ref:`runSModelS.py <runSModelS>`.
In addition, SModelS is also usable as a :ref:`library <exampleCode>`, providing functionality to

 * decompose models into simplified models,
 * confront input models with LHC constraints,
 * compute LO, NLO, NLL SUSY cross sections,
 * identify missing topologies,
 * browse the SModelS database of SMS results
 * and a few more tasks.

For example code for various tasks, see :ref:`Examples`

**References**:

 [1] https://twiki.cern.ch/twiki/bin/view/AtlasPublic/SupersymmetryPublicResults

 [2] https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS

 [3] Heavy Flavor Averaging Group Collaboration, Y. Amhis et al.,
 Averages of B-Hadron, C-Hadron, and tau-lepton properties as of
 early 2012, 1207.1158. and HFAG2013 update at 
 http://www.slac.stanford.edu/xorg/hfag/

 [4] LHCb Collaboration, R. Aaij et al., First evidence for the decay
 B s → μ+ μ- , Phys. Rev. Lett. 110, 021801 (2013), 1211.2674

 [5] ATLAS Collaboration, G. Aad et al., Observation of a new particle
 in the search for the Standard Model Higgs boson with the ATLAS
 detector at the LHC. Phys. Lett. B 716, 1–29 (2012), 1207.7214

 [6] CMS Collaboration, S. Chatrchyan et al., Observation of a new
 boson at a mass of 125 GeV with the CMS experiment at the LHC.
 Phys. Lett. B 716, 30–61 (2012), 1207.7235

 [7] S. Sekmen, S. Kraml, J. Lykken, F. Moortgat, S. Padhi, et al.,
 Interpreting LHC SUSY searches in the phenomenological MSSM.
 JHEP 1202, 075 (2012), 1109.5119

 [8] A. Arbey, M. Battaglia, F. Mahmoudi, Implications of LHC
 searches on SUSY particle spectra: the pMSSM parameter space
 with neutralino dark matter. Eur. Phys. J. C 72, 1847 (2012),
 1110.3726

 [9] M. Papucci, J.T. Ruderman, A. Weiler, Natural SUSY endures.
 JHEP 1209, 035 (2012), 1110.6926

 [10] M.W. Cahill-Rowley, J.L. Hewett, A. Ismail, T.G. Rizzo, More
 energy, more searches, but the pMSSM lives on, Phys. Rev. D 88,
 035002 (2013), 1211.1981

 [11] H.K. Dreiner, M. Krämer, J. Tattersall, How low can SUSY go?
 Matching, monojets and compressed spectra. Europhys. Lett. 99,
 61001 (2012), 1207.1613

 [12] R. Mahbubani, M. Papucci, G. Perez, J.T. Ruderman, A. Weiler,
 Light non-degenerate squarks at the LHC. Phys. Rev. Lett. 110,
 151804 (2013), 1212.3328

 [13] CMS Collaboration, Phenomenological MSSM interpretation of
 the CMS 2011 5fb-1 results, Tech. Rep. CMS-PAS-SUS-12-030,
 CERN, Geneva (2013)

 [14] ATLAS Collaboration, H. Okawa, Interpretations of SUSY
 searches in ATLAS with simplified models, 1110.0282
 Page 21 of 23 2868

 [15] CMS Collaboration, S. Chatrchyan et al., Interpretation of searches
 for supersymmetry with simplified models. Phys. Rev. D 88,
 052017 (2013), 1301.2175

 [16] J. Alwall, P. Schuster, N. Toro, Simplified models for a first characterization of new physics at the LHC. Phys. Rev. D 79, 075020 (2009), 0810.3921

 [17] LHC New Physics Working Group Collaboration, D. Alves et al.,
 Simplified models for LHC new physics searches, J. Phys. G 39,
 105005 (2012), 1105.2838
