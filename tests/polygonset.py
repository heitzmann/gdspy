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
            gdspy.LayoutViewer(cells=[c1, c2])
        assert result is None

def test_polygonset0(target):
    ps = gdspy.PolygonSet([
        [(10, 0), (11, 0), (10, 1)],
        numpy.array([(11.0, 0), (10, 1), (11, 1)]),
        [numpy.array((11, 1)), numpy.array((12.0, 1.0)), (11, 2)],
    ], 1, 2)
    assert str(ps) == "PolygonSet (3 polygons, 9 vertices, layers [1], datatypes [2])"
    bb = ps.get_bounding_box()
    assert bb.shape == (2, 2)
    assert numpy.max(numpy.abs(bb - numpy.array(((10, 0), (12, 2))))) == 0
    assert gdspy.PolygonSet([]).get_bounding_box() == None
    cell = gdspy.Cell('test', True).add(ps)
    assertsame(cell, target['PolygonSet0'])

def test_polygonset1(target):
    ps = gdspy.PolygonSet([
        [(10, 0), (11, 0), (10, 1)],
        numpy.array([(11.0, 0), (10, 1), (11, 1)]),
        [numpy.array((11, 1)), numpy.array((12.0, 1.0)), (11, 2)],
    ], 2, 3)
    ps.rotate(numpy.pi / 3, center=(10, 1))
    cell = gdspy.Cell('test', True).add(ps)
    assertsame(cell, target['PolygonSet1'])

def test_polygonset2(target):
    ps = gdspy.PolygonSet([
        [(10, 0), (11, 0), (10, 1)],
        numpy.array([(11.0, 0), (10, 1), (11, 1)]),
        [numpy.array((11, 1)), numpy.array((12.0, 1.0)), (11, 2)],
    ], 3, 4)
    ps.scale(0.5)
    ps.scale(1, 2, center=(5, 1))
    cell = gdspy.Cell('test', True).add(ps)
    assertsame(cell, target['PolygonSet2'])
