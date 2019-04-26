######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

from tutils import target, assertsame
import numpy
import pytest
import gdspy


def test_robustpath1(target):
    cell = gdspy.Cell('test', True)
    rp = gdspy.RobustPath((0, 0), 0.1, layer=[1], gdsii_path=True)
    rp.segment((1, 1))
    cell.add(rp)
    rp = gdspy.RobustPath((1, 0), 0.1, [-0.1, 0.1], tolerance=1e-5,
                          ends=['round', 'extended'], layer=[2, 3], max_points=6)
    rp.segment((2, 1))
    cell.add(rp)
    rp = gdspy.RobustPath((2, 0), [0.1, 0.2], 0.2, ends=(0.2, 0.1),
                          layer=4, datatype=[1, 1])
    rp.segment((3, 1))
    cell.add(rp)
    rp = gdspy.RobustPath((3, 0), [0.1, 0.2, 0.1], [-0.2, 0, 0.2],
                          ends=[(0.2, 0.1), 'smooth', 'flush'], datatype=5)
    rp.segment((4, 1))
    cell.add(rp)
    assertsame(target['RobustPath1'], cell)

def test_robustpath2(target):
    cell = gdspy.Cell('test', True)
    rp = gdspy.RobustPath((0, 0), [0.1, 0.2, 0.1], 0.15, layer=[1, 2, 3])
    rp.segment((1, 0))
    rp.segment((1, 1), 0.1, 0.05)
    rp.segment((1, 1), [0.2, 0.1, 0.1], -0.05, True)
    rp.segment((-1, 1), 0.2, [-0.2, 0, 0.3], True)
    rp.arc(2, 0, 0.5 * numpy.pi)
    rp.arc(3, 0.7 * numpy.pi, numpy.pi, 0.1, 0)
    rp.arc(2, 0.4 * numpy.pi, -0.4 * numpy.pi, [0.1, 0.2, 0.1], [0.2, 0, -0.2])
    rp.turn(1, -0.3 * numpy.pi)
    rp.turn(1, 'rr', 0.15)
    rp.turn(0.5, 'l', [0.05, 0.1, 0.05], [0.15, 0, -0.15])
    cell.add(rp)
    rp = gdspy.RobustPath((-5, 6), 0.8, layer=20, ends='round', tolerance=1e-4)
    rp.segment((1, 1), 0.1, relative=True)
    cell.add(rp)
    rp = gdspy.RobustPath((-5, 6), 0.8, layer=21, ends='extended', tolerance=1e-4)
    rp.segment((1, 1), 0.1, relative=True)
    cell.add(rp)
    rp = gdspy.RobustPath((-5, 6), 0.8, layer=22, ends=(0.1, 0.2), tolerance=1e-4)
    rp.segment((1, 1), 0.1, relative=True)
    cell.add(rp)
    rp = gdspy.RobustPath((-5, 6), 0.8, layer=23, ends='smooth', tolerance=1e-4)
    rp.segment((1, 1), 0.1, relative=True)
    cell.add(rp)
    rp = gdspy.RobustPath((-3, 6), 0.8, layer=10, ends='round', tolerance=1e-5)
    rp.segment((1, 0), 0.1, relative=True)
    rp.segment((0, 1), 0.8, relative=True)
    cell.add(rp)
    rp = gdspy.RobustPath((-3, 6), 0.8, layer=11, ends='extended', tolerance=1e-5)
    rp.segment((1, 0), 0.1, relative=True)
    rp.segment((0, 1), 0.8, relative=True)
    cell.add(rp)
    rp = gdspy.RobustPath((-3, 6), 0.8, layer=12, ends='smooth', tolerance=1e-5)
    rp.segment((1, 0), 0.1, relative=True)
    rp.segment((0, 1), 0.8, relative=True)
    cell.add(rp)
    rp = gdspy.RobustPath((-3, 8), 0.1, layer=13, ends='round', tolerance=1e-5)
    rp.segment((1, 0), 0.8, relative=True)
    rp.segment((0, 1), 0.1, relative=True)
    cell.add(rp)
    rp = gdspy.RobustPath((-3, 8), 0.1, layer=14, ends=(0.2, 0.2), tolerance=1e-5)
    rp.segment((1, 0), 0.8, relative=True)
    rp.segment((0, 1), 0.1, relative=True)
    cell.add(rp)
    rp = gdspy.RobustPath((-3, 8), 0.1, layer=15, ends='smooth', tolerance=1e-5)
    rp.segment((1, 0), 0.8, relative=True)
    rp.segment((0, 1), 0.1, relative=True)
    cell.add(rp)
    rp = gdspy.RobustPath((5, 2), [0.05, 0.1, 0.2], [-0.2, 0, 0.4], layer=[4, 5, 6])
    rp.parametric(lambda u: numpy.array((5.5 + 3 * u, 2 + 3 * u**2)), relative=False)
    rp.segment((0, 1), relative=True)
    rp.parametric(lambda u: numpy.array((2 * numpy.cos(0.5 * numpy.pi * u) - 2,
                                         3 * numpy.sin(0.5 * numpy.pi * u))),
                  width=[0.2, 0.1, 0.05], offset=[-0.3, 0, 0.3])
    rp.parametric(lambda u: numpy.array((-2*u, 0)), width=0.1, offset=0.2)
    rp.bezier([(-3, 0), (-2, -3), (0, -4), (0, -5)], offset=[-0.2, 0, 0.2])
    rp.bezier([(4.5, 0), (1, -1),  (1, 5), (3, 2), (5, 2)], width=[0.05, 0.1, 0.2],
              offset=[-0.2, 0, 0.4], relative=False)
    cell.add(rp)
    rp = gdspy.RobustPath((2, -1), 0.1, layer=7, tolerance=1e-4, max_points=0)
    rp.smooth([(1, 0), (1, -1), (0, -1)], angles=[numpy.pi / 3, None, -2 / 3 * numpy.pi, None],
              cycle=True)
    cell.add(rp)
    rp = gdspy.RobustPath((2.5, -1.5), 0.1, layer=8)
    rp.smooth([(3, -1.5), (4, -2), (5, -1), (6, -2), (7, -1.5), (7.5, -1.5)], relative=False,
              width=0.2)
    cell.add(rp)
    assertsame(target['RobustPath2'], cell)
