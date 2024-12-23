"""
.. module:: SMparticleDefinitions
   :synopsis: Defines the SM particles.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.share.models.SMparticles import *
from smodels.base.vertexGraph import VertexGraph



# SM vertices
SMvertices = []

# Gauge vertices
for p in [e,mu,ta,u,c,t,d,s,b,W]:
    v_photon = VertexGraph(incoming=[photon],outgoing=[p,p.chargeConjugate()])
    SMvertices.append(v_photon)

for p in [e,mu,ta,u,c,t,d,s,b,nue,numu,nuta,W]:    
    v_Z = VertexGraph(incoming=[Z],outgoing=[p,p.chargeConjugate()])
    SMvertices.append(v_Z)
    
for p1,p2 in zip([e,mu,ta,d,s,b],[nue,numu,nuta,u,c,t]):
    v_W = VertexGraph(incoming=[W],outgoing=[p1.chargeConjugate(),p2])
    SMvertices.append(v_W)

for p in [u,c,t,d,s,b,g]:
    v_gluon = VertexGraph(incoming=[g],outgoing=[p,p.chargeConjugate()])
    SMvertices.append(v_gluon)

for p in [t,b,g,photon,W,Z]:
    v_higgs = VertexGraph(incoming=[higgs],outgoing=[p,p.chargeConjugate()])
    SMvertices.append(v_higgs)
