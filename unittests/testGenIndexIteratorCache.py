#!/usr/bin/env python3

"""
.. module:: testGenIndexIteratorCache
   :synopsis: Tests cache invalidation in GenericSMS.genIndexIterator

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest

from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS
from smodels.base.particle import Particle
from smodels.base.particleNode import ParticleNode


class TestGenIndexIteratorCache(unittest.TestCase):

    @staticmethod
    def _node(label, isInclusive=False):
        particle = Particle(label=label, isSM=False)
        return ParticleNode(particle=particle, isInclusive=isInclusive)

    @staticmethod
    def _build_sms(sms_cls):
        sms = sms_cls()
        pv = sms.add_node(TestGenIndexIteratorCache._node("PV"))
        a = sms.add_node(TestGenIndexIteratorCache._node("A"))
        b = sms.add_node(TestGenIndexIteratorCache._node("B"))
        c = sms.add_node(TestGenIndexIteratorCache._node("C"))
        sms.add_edges_from([(pv, a), (pv, b), (a, c)])
        return sms, pv, a, b, c

    @staticmethod
    def _mothers(sms, includeLeaves=True, ignoreInclusiveNodes=False):
        return [mom for mom, _ in sms.genIndexIterator(includeLeaves=includeLeaves,
                                                       ignoreInclusiveNodes=ignoreInclusiveNodes)]

    def test_cache_invalidation_on_topology_and_node_updates(self):
        for sms_cls in (TheorySMS, ExpSMS):
            with self.subTest(smsClass=sms_cls.__name__):
                sms, pv, a, b, c = self._build_sms(sms_cls)

                mothers = self._mothers(sms, includeLeaves=True)
                self.assertEqual(mothers, [pv, a, b, c])

                sms.remove_edge(a, c)
                mothers_after_remove = self._mothers(sms, includeLeaves=True)
                self.assertEqual(mothers_after_remove, [pv, a, b])

                sms.add_edge(b, c)
                mothers_after_add = self._mothers(sms, includeLeaves=True)
                self.assertEqual(mothers_after_add, [pv, a, b, c])

                mothers_non_inclusive = self._mothers(sms, includeLeaves=True,
                                                      ignoreInclusiveNodes=True)
                self.assertEqual(mothers_non_inclusive, [pv, a, b, c])

                sms.updateNodeObjects({a: self._node("Ainc", isInclusive=True)})
                mothers_ignore_inclusive = self._mothers(sms, includeLeaves=True,
                                                         ignoreInclusiveNodes=True)
                self.assertEqual(mothers_ignore_inclusive, [pv, b, c])


if __name__ == "__main__":
    unittest.main()
