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
    fp = gdspy.FlexPath([(0, 0), (0.5, 0), (1, 0), (1, 1), (0, 1), (-1, -2), (-2, 0)], 0.05, [0, -0.1, 0, 0.1],
                        corners=['natural', 'circular bend', 'circular bend', 'circular bend'],
                        ends=['flush', 'extended', (0.1, 0.2), 'round'], tolerance=1e-4,
                        layer=[0, 1, 1, 2], bend_radius=[0, 0.3, 0.3, 0.2], max_points=10)
    cell.add(fp)
    assertsame(target['FlexPath2'], cell)
    cell = gdspy.Cell('test2', True)

def test_flexpath3(target):
    cell = gdspy.Cell('test', True)
    pts = numpy.array([(0, 0), (0.5, 0), (1, 0), (1, 2), (3, 0), (2, -1), (2, -2), (0, -1),
                       (1, -2), (1, -3)])
    fp = gdspy.FlexPath(pts + numpy.array((0, 5)), [0.1, 0.1, 0.1], 0.15, layer=[1, 2, 3],
                        corners=['natural', 'miter', 'bevel'], ends=(0.5, 0))
    cell.add(fp)
    fp = gdspy.FlexPath(pts + numpy.array((5, 0)), [0.1, 0.1, 0.1], 0.15, layer=[4, 5, 6],
                        corners=['round', 'smooth', broken],
                        ends=[pointy, 'smooth', (0, 0.5)])
    cell.add(fp)
    assertsame(target['FlexPath3'], cell)

def test_flexpath4(target):
    cell = gdspy.Cell('test', True)
    fp = gdspy.FlexPath([(0, 0)], [0.1, 0.2, 0.1], 0.15, layer=[1, 2, 3],
                        corners=['natural', 'miter', 'bevel'])
    fp.segment((1, 0))
    fp.segment((1, 1), 0.1, 0.05)
    fp.segment((1, 1), [0.2, 0.1, 0.1], -0.05, True)
    fp.segment((-1, 1), 0.2, [-0.2, 0, 0.3], True)
    fp.arc(2, 0, 0.5 * numpy.pi)
    fp.arc(3, 0.5 * numpy.pi, numpy.pi, 0.1, 0)
    fp.arc(1, 0.4 * numpy.pi, -0.4 * numpy.pi, [0.1, 0.2, 0.1], [0.2, 0, -0.2])
    fp.turn(1, 0.4 * numpy.pi)
    fp.turn(1, 'll', 0.15, 0)
    fp.turn(0.5, 'r', [0.1, 0.05, 0.1], [0.15, 0, -0.15])
    cell.add(fp)
    fp = gdspy.FlexPath([(-5, 6)], 0.8, layer=20, ends='round', tolerance=1e-4)
    fp.segment((1, 1), 0.1, relative=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(-5, 6)], 0.8, layer=21, ends='extended', tolerance=1e-4)
    fp.segment((1, 1), 0.1, relative=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(-5, 6)], 0.8, layer=22, ends=(0.1, 0.2), tolerance=1e-4)
    fp.segment((1, 1), 0.1, relative=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(-5, 6)], 0.8, layer=23, ends='smooth', tolerance=1e-4)
    fp.segment((1, 1), 0.1, relative=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(-3, 6)], 0.8, layer=10, corners='round', ends='round', tolerance=1e-5)
    fp.segment((1, 0), 0.1, relative=True)
    fp.segment((0, 1), 0.8, relative=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(-3, 6)], 0.8, layer=11, corners='smooth', ends='extended', tolerance=1e-5)
    fp.segment((1, 0), 0.1, relative=True)
    fp.segment((0, 1), 0.8, relative=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(-3, 6)], 0.8, layer=12, corners='smooth', ends='smooth', tolerance=1e-5)
    fp.segment((1, 0), 0.1, relative=True)
    fp.segment((0, 1), 0.8, relative=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(-3, 8)], 0.1, layer=13, corners='round', ends='round', tolerance=1e-5)
    fp.segment((1, 0), 0.8, relative=True)
    fp.segment((0, 1), 0.1, relative=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(-3, 8)], 0.1, layer=14, corners='smooth', ends=(0.2, 0.2), tolerance=1e-5)
    fp.segment((1, 0), 0.8, relative=True)
    fp.segment((0, 1), 0.1, relative=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(-3, 8)], 0.1, layer=15, corners='round', ends='smooth', tolerance=1e-5)
    fp.segment((1, 0), 0.8, relative=True)
    fp.segment((0, 1), 0.1, relative=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(5, 2)], [0.05, 0.1, 0.2], [-0.2, 0, 0.4], layer=[4, 5, 6])
    fp.parametric(lambda u: numpy.array((5.5 + 3 * u, 2 + 3 * u**2)), relative=False)
    fp.segment((0, 1), relative=True)
    fp.parametric(lambda u: numpy.array((2 * numpy.cos(0.5 * numpy.pi * u) - 2,
                                         3 * numpy.sin(0.5 * numpy.pi * u))),
                  [0.2, 0.1, 0.05], [-0.3, 0, 0.3])
    fp.parametric(lambda u: numpy.array((-2*u, 0)), 0.1, 0.2)
    fp.bezier([(-3, 0), (-2, -3), (0, -4), (0, -5)], offset=[-0.2, 0, 0.2])
    fp.bezier([(5, 0), (1, -1),  (1, 5), (3, 2), (5, 2)], [0.05, 0.1, 0.2],
              [-0.2, 0, 0.4], relative=False)
    cell.add(fp)
    fp = gdspy.FlexPath([(2, -1)], 0.1, layer=7, tolerance=1e-5, max_points=0)
    fp.smooth([(1, 0), (1, -1), (0, -1)], angles=[numpy.pi / 3, None, -2 / 3.0 * numpy.pi, None],
              cycle=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(2.5, -1.5)], 0.1, layer=8)
    fp.smooth([(3, -1.5), (4, -2), (5, -1), (6, -2), (7, -1.5), (7.5, -1.5)], relative=False,
              width=0.2)
    cell.add(fp)
    assertsame(target['FlexPath4'], cell)

def test_flexpath_gdsiipath():
    cells = []
    for gdsii_path in [True, False]:
        cells.append(gdspy.Cell(str(gdsii_path), True))
        fp = gdspy.FlexPath([(0, 0), (0.5, 0), (1, 0), (1, 1), (0, 1), (-1, -2), (-2, 0)], 0.05, [-0.1, 0.1],
                            corners=['natural', 'circular bend'], ends=['extended', (0.1, 0.2)],
                            layer=[0, 1,], bend_radius=[0, 0.3], gdsii_path=gdsii_path)
        cells[-1].add(fp)
    assertsame(*cells)

def test_flexpath_getpolygons():
    fp = gdspy.FlexPath([(0, 0), (0.5, 0), (1, 0), (1, 1), (0, 1), (-1, -2), (-2, 0)],
                        0.05, [-0.1, 0.1], corners=['natural', 'circular bend'],
                        ends=['extended', (0.1, 0.2)], bend_radius=[0, 0.3],
                        layer=[0, 1], datatype=[1, 0])
    d = fp.get_polygons(True)
    l = fp.get_polygons()
    assert len(d) == 2
    assert (1, 0) in d
    assert (0, 1) in d
    assert sum(len(p) for p in d.values()) == len(l)
    assert sum(sum(len(x) for x in p) for p in d.values()) == sum(len(x) for x in l)
    ps = fp.to_polygonset()
    assert len(ps.layers) == len(ps.datatypes) == len(ps.polygons)
    assert gdspy.FlexPath([(0, 0)], 1).to_polygonset() == None

def test_flexpath_togds(tmpdir):
    cell = gdspy.Cell('flexpath')
    fp = gdspy.FlexPath([(0, 0), (0.5, 0.5), (1, 1)], 0.1, layer=[1], gdsii_path=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(3, 0), (3.5, 0.5), (4, 1)], 0.1, layer=[21])
    cell.add(fp)
    fp = gdspy.FlexPath([(1, 0), (2, 1)], 0.1, [-0.1, 0.1], max_points=6,
                        ends=['round', 'extended'], layer=[2, 3], gdsii_path=True)
    cell.add(fp)
    fp = gdspy.FlexPath([(2, 0), (3, 1)], [0.1, 0.2], 0.2, ends=(0.2, 0.1),
                        layer=4, datatype=[10, 10], gdsii_path=True)
    cell.add(fp)

    fp = gdspy.FlexPath([(0, 0), (0.5, 0), (1, 0), (1, 1), (0, 1), (-1, -2), (-2, 0)],
                        0.05, [0, -0.1, 0, 0.1],
                        corners=['natural', 'circular bend', 'circular bend',
                        'circular bend'], tolerance=1e-5,
                        ends=['flush', 'extended', (0.1, 0.2), 'flush'],
                        layer=[10, 11, 11, 12], bend_radius=[0, 0.3, 0.3, 0.2],
                        gdsii_path=True).translate(-5, 0)
    cell.add(fp)
    fp = gdspy.FlexPath([(i, 2 + i**2) for i in numpy.linspace(0, 1, 8192)],
                        0.01, gdsii_path=True)
    cell.add(fp)
    fname = str(tmpdir.join('test.gds'))
    with pytest.warns(UserWarning):
        gdspy.write_gds(fname, unit=1, precision=1e-7)
    lib = gdspy.GdsLibrary(infile=fname, rename={'flexpath': 'file'})
    assertsame(lib.cell_dict['file'], cell, tolerance=1e-3)
    gdspy.current_library = gdspy.GdsLibrary()
