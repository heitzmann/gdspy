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


def test_polygonset(target):
    ps = gdspy.PolygonSet([
        [(10, 0), (11, 0), (10, 1)],
        numpy.array([(11.0, 0), (10, 1), (11, 1)]),
        [numpy.array((11, 1)), numpy.array((12.0, 1.0)), (11, 2)],
    ], 1, 2)
    assert str(ps) == "PolygonSet (3 polygons, 9 vertices, layers [1], datatypes [2])"
    assert ps.area() == 1.5
    assert ps.area(True) == {(1,2): 1.5}
    bb = ps.get_bounding_box()
    assert bb.shape == (2, 2)
    assert numpy.max(numpy.abs(bb - numpy.array(((10, 0), (12, 2))))) == 0
    assert gdspy.PolygonSet([]).get_bounding_box() == None
    cell = gdspy.Cell('test', True).add(ps)
    assertsame(cell, target['PolygonSet'])

def test_translate():
    ps = gdspy.PolygonSet([[(0, 0), (1, 0), (1, 2), (0, 2)]])
    ps.translate(-1, -2)
    tgt = gdspy.PolygonSet([[(-1, -2), (0, -2), (0, 0), (-1, 0)]])
    assertsame(gdspy.Cell('test', True).add(ps), gdspy.Cell('TGT', True).add(tgt))

def test_rotation():
    ps = gdspy.PolygonSet([[(0, 0), (1, 0), (1, 2), (0, 2)]])
    ps.rotate(-numpy.pi / 6, center=(1, 0))
    x = 3**0.5 / 2
    a = numpy.arctan(2) + numpy.pi / 6
    l = 5**0.5
    tgt = gdspy.PolygonSet([[(1 - x, 0.5), (1, 0), (2, 2 * x),
                             (1 - l * numpy.cos(a), l * numpy.sin(a))]])
    assertsame(gdspy.Cell('test', True).add(ps), gdspy.Cell('TGT', True).add(tgt))

def test_mirror():
    ps = gdspy.PolygonSet([[(0, 0), (1, 0), (1, 2), (0, 2)]])
    ps.mirror((-1, 1), (1, -1))
    tgt = gdspy.PolygonSet([[(0, 0), (-2, 0), (-2, -1), (0, -1)]])
    assertsame(gdspy.Cell('test', True).add(ps), gdspy.Cell('TGT', True).add(tgt))

def test_scale():
    ps = gdspy.PolygonSet([[(0, 0), (1, 0), (1, 2), (0, 2)]])
    cell = gdspy.Cell('test', True).add(ps)
    ps.scale(0.5)
    tgt = gdspy.PolygonSet([[(0, 0), (0.5, 0), (0.5, 1), (0, 1)]])
    assertsame(cell, gdspy.Cell('TGT', True).add(tgt))
    ps = gdspy.PolygonSet([[(0, 0), (1, 0), (1, 2), (0, 2)]])
    cell = gdspy.Cell('test', True).add(ps)
    ps.scale(0.2, 2, center=(1, 2))
    tgt = gdspy.PolygonSet([[(0.8, -2), (1, -2), (1, 2), (0.8, 2)]])
    assertsame(cell, gdspy.Cell('TGT', True).add(tgt))

def test_togds(tmpdir):
    ps = gdspy.PolygonSet([
        [(10 + i * 1e-3, i**2 * 1e-6) for i in range(8191)],
    ])
    cell = gdspy.Cell('LARGE', True).add(ps)
    fname = str(tmpdir.join('large.gds'))
    lib =gdspy.GdsLibrary(unit=1, precision=1e-7).add(cell)
    with pytest.warns(UserWarning):
        lib.write_gds(fname)
    lib = gdspy.GdsLibrary(infile=fname)
    assertsame(lib.cell_dict['LARGE'], cell)

def test_fracture():
    ps1 = gdspy.PolygonSet([[(-2, -1), (-2, 1), (-1, 2), (1, 2),
                             (2, 1), (2, -1), (1, -2), (-1, -2)]])
    ps2 = gdspy.PolygonSet([[(-2, -1), (-2, 1), (-1, 2), (1, 2),
                             (2, 1), (2, -1), (1, -2), (-1, -2)]])
    ps2.fracture(5)
    assert all(len(p) <= 5 for p in ps2.polygons)
    assertsame(gdspy.Cell('1', True).add(ps1), gdspy.Cell('2', True).add(ps2))

    def spiral(u):
        r = 4 - 3 * u
        theta = 5 * u * numpy.pi
        x = r * numpy.cos(theta) - 4
        y = r * numpy.sin(theta)
        return (x, y)

    def dspiral_dt(u):
        theta = 5 * u * numpy.pi
        dx_dt = -numpy.sin(theta)
        dy_dt = numpy.cos(theta)
        return (dx_dt, dy_dt)

    pt1 = gdspy.Path(0.5, (0, 0))
    pt1.parametric(spiral, dspiral_dt, number_of_evaluations=512, tolerance=10, max_points=0)
    assert len(pt1.polygons) == 1
    assert len(pt1.polygons[0]) == 1024
    pt2 = gdspy.Path(0.5, (0, 0))
    pt2.parametric(spiral, dspiral_dt, number_of_evaluations=512, tolerance=10, max_points=0)
    pt2.fracture(199, precision=1e-6)
    assert all(len(p) <= 199 for p in pt2.polygons)
    assertsame(gdspy.Cell('3', True).add(pt1), gdspy.Cell('4', True).add(pt2))

def test_fillet(target):
    cell = gdspy.Cell('test')
    orig = gdspy.PolygonSet([
        [(0, 0), (-1, 0), (0, -1), (0.5, -0.5), (1, 0), (1, 1), (4, -1), (1, 3), (1, 2), (0, 1)],
        [(2, -1), (3, -1), (2.5, -2)],
    ])
    orig.datatypes=[0, 1]
    p = gdspy.copy(orig, 0, 5)
    p.layers = [1, 1]
    p.fillet(0.3, max_points=0)
    cell.add(p)
    p = gdspy.copy(orig, 5, 5)
    p.layers = [2, 2]
    p.fillet([0.3, 0.2, 0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.1, 0.2, 0], max_points=0)
    cell.add(p)
    p = gdspy.copy(orig, 5, 0)
    p.layers = [3, 3]
    p.fillet([[0.1, 0.1, 0.4, 0, 0.4, 0.1, 0.1, 0.4, 0.4, 0.1], [0.2, 0.2, 0.5]], max_points=0)
    cell.add(p)
    p = gdspy.PolygonSet([
        [(0, 0), (0, 0), (-1, 0), (0, -1), (0.5, -0.5), (1, 0), (1, 0), (1, 1), (4, -1), (1, 3), (1, 2), (0, 1)],
        [(2, -1), (3, -1), (2.5, -2), (2, -1)],
    ], layer=4)
    p.datatypes=[0, 1]
    p.fillet([0.8, [10.0, 10.0, 20.0, 20.0]], max_points=199, precision=1e-6)
    cell.add(p)
    assertsame(cell, target['PolygonSet_fillet'])
