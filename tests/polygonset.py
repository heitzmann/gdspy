######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import numpy
import pytest
import gdspy


@pytest.fixture
def target():
    return gdspy.GdsLibrary(infile='tests/test.gds').cell_dict

def assertsame(c1, c2):
    d1 = c1.get_polygons(by_spec=True)
    d2 = c2.get_polygons(by_spec=True)
    for key in d1:
        assert key in d2
        result = gdspy.boolean(d1[key], d2[key], 'xor', precision=1e-6, layer=key[0] + 1)
        if result is not None:
            c1.add(result)
            c2.add(result)
            result = gdspy.offset(result, -1e-6, precision=1e-7)
            if result is not None:
                gdspy.LayoutViewer(cells=[c1, c2])
        assert result is None

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
    assertsame(cell, target['PolygonSet0'])

def test_rotation(target):
    ps = gdspy.PolygonSet([
        [(10, 0), (11, 0), (10, 1)],
        numpy.array([(11.0, 0), (10, 1), (11, 1)]),
        [numpy.array((11, 1)), numpy.array((12.0, 1.0)), (11, 2)],
    ], 2, 3)
    ps.rotate(numpy.pi / 3, center=(10, 1))
    cell = gdspy.Cell('test', True).add(ps)
    assertsame(cell, target['PolygonSet1'])

def test_scale(target):
    ps = gdspy.PolygonSet([
        [(10, 0), (11, 0), (10, 1)],
        numpy.array([(11.0, 0), (10, 1), (11, 1)]),
        [numpy.array((11, 1)), numpy.array((12.0, 1.0)), (11, 2)],
    ], 3, 4)
    ps.scale(0.5)
    ps.scale(1, 2, center=(5, 1))
    cell = gdspy.Cell('test', True).add(ps)
    assertsame(cell, target['PolygonSet2'])

def test_polygonset3(tmpdir):
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

def test_fracture(target):
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
