#!/usr/bin/env python
# module TEST_RUN_SIMULATION

import unittest
import numpy as np
from pylocus.point_set import HeterogenousSet

from run_simulation import get_Om_from_abs_angles


class TestSimulation(unittest.TestCase):
    def setUp(self):
        N = 4
        d = 2
        self.points = HeterogenousSet(N, d)
        self.points.set_points('normal')

    def test_multiple(self):
        for i in range(100):
            print('running test number', i)
            self.points.set_points('normal')
            self.test_Om_from_abs_angles()
            self.test_constraints()

    def test_Om_from_abs_angles(self):
        Om = get_Om_from_abs_angles(self.points.abs_angles, self.points.m)
        self.assertTrue(np.allclose(Om, self.points.Om),
                        'Om from all methods does not give the same result.')

    def test_constraints(self):
        C, b = self.points.get_KE_constraints()
        self.assertTrue(np.allclose(C.dot(self.points.KE), b))


if __name__ == "__main__":
    unittest.main()
