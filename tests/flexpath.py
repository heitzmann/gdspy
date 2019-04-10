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


def broken(p0, v0, p1, v1, p2, w):
    den = v1[1] * v0[0] - v1[0] * v0[1]
    lim = 1e-12 * (v0[0]**2 + v0[1]**2) * (v1[0]**2 + v1[1]**2)
    if den**2 < lim:
        u0 = u1 = 0
        p = 0.5 * (p0 + p1)
    else:
        dx = p1[0] - p0[0]
        dy = p1[1] - p0[1]
        u0 = (v1[1] * dx - v1[0] * dy) / den
        u1 = (v0[1] * dx - v0[0] * dy) / den
        p = 0.5 * (p0 + v0 * u0 + p1 + v1 * u1)
    if u0 <= 0 and u1 >= 0:
        return [p]
    return [p0, p2, p1]

def pointy(p0, v0, p1, v1):
    r = 0.5 * numpy.sqrt(numpy.sum((p0 - p1)**2))
    v0 /= numpy.sqrt(numpy.sum(v0**2))
    v1 /= numpy.sqrt(numpy.sum(v1**2))
    return [p0, 0.5 * (p0 + p1) + 0.5 * (v0 - v1) * r, p1]

def test_flexpath_warnings():
    f = lambda *args: None
    for ends in ['smooth', f]:
        with pytest.warns(UserWarning):
            gdspy.FlexPath([(0, 0)], 1, ends=ends, gdsii_path=True)
    for corners in ['miter', 'bevel', 'round', 'smooth', f]:
        with pytest.warns(UserWarning):
            gdspy.FlexPath([(0, 0)], 1, corners=corners, gdsii_path=True)

def test_flexpath1(target):
    cell = gdspy.Cell('test', True)
    fp = gdspy.FlexPath([(0, 0), (1, 1)], 0.1, layer=[1], gdsii_path=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(1, 0), (2, 1)], 0.1, [-0.1, 0.1], tolerance=1e-5,
                        ends=['round', 'extended'], layer=[2, 3], max_points=6)
    cell.add(fp)
    fp = gdspy.FlexPath([(2, 0), (3, 1)], [0.1, 0.2], 0.2, ends=(0.2, 0.1),
                        layer=4, datatype=[1, 1])
    cell.add(fp)
    fp = gdspy.FlexPath([(3, 0), (4, 1)], [0.1, 0.2, 0.1], [-0.2, 0, 0.2],
                        ends=[(0.2, 0.1), 'smooth', pointy], datatype=5)
    cell.add(fp)
    assertsame(target['FlexPath1'], cell)

def test_flexpath2(target):
    cell = gdspy.Cell('test', True)
    fp = gdspy.FlexPath([(0, 0), (1, 0), (1, 1), (0, 1), (-1, -2), (-2, 0)], 0.05, [0, -0.1, 0, 0.1],
                        corners=['natural', 'circular bend', 'circular bend', 'circular bend'],
                        ends=['flush', 'extended', (0.1, 0.2), 'round'], tolerance=1e-4,
                        layer=[0, 1, 1, 2], bend_radius=[0, 0.3, 0.3, 0.2], max_points=10)
    cell.add(fp)
    assertsame(target['FlexPath2'], cell)
    cell = gdspy.Cell('test2', True)

def test_flexpath3(target):
    cell = gdspy.Cell('test', True)
    pts = numpy.array([(0, 0), (1, 0), (1, 2), (3, 0), (2, -1), (2, -2), (0, -1),
                       (1, -2), (1, -3)])
    fp = gdspy.FlexPath(pts + numpy.array((0, 5)), [0.1, 0.1, 0.1], 0.15, layer=[1, 2, 3],
                        corners=['natural', 'miter', 'bevel'], ends=(0.5, 0))
    cell.add(fp)
    fp = gdspy.FlexPath(pts + numpy.array((5, 0)), [0.1, 0.1, 0.1], 0.15, layer=[4, 5, 6],
                        corners=['round', 'smooth', broken],
                        ends=[pointy, 'smooth', (0, 0.5)])
    cell.add(fp)
    assertsame(target['FlexPath3'], cell)

def test_flexpath_gdsiipath():
    cells = []
    for gdsii_path in [True, False]:
        cells.append(gdspy.Cell(str(gdsii_path), True))
        fp = gdspy.FlexPath([(0, 0), (1, 0), (1, 1), (0, 1), (-1, -2), (-2, 0)], 0.05, [-0.1, 0.1],
                            corners=['natural', 'circular bend'], ends=['extended', (0.1, 0.2)],
                            layer=[0, 1,], bend_radius=[0, 0.3], gdsii_path=gdsii_path)
        cells[-1].add(fp)
    assertsame(*cells)
