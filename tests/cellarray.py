######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import pytest
import gdspy
import numpy


def test_noreference():
    name = 'ca_noreference'
    with pytest.warns(UserWarning):
        ref = gdspy.CellArray(name, 2, 3, (3, 2), (1, -1), 90, 2, True)
    ref.translate(-1, 1)
    assert ref.ref_cell == name
    assert ref.area() == 0
    assert ref.area(True) == dict()
    assert ref.get_bounding_box() is None
    assert ref.get_polygons() == []
    assert ref.get_polygons(True) == dict()
    assert ref.origin[0] == ref.origin[1] == 0


def test_empty():
    name = 'ca_empty'
    c = gdspy.Cell(name)
    ref = gdspy.CellArray(name, 2, 3, (3, 2), (1, -1), 90, 2, True)
    ref.translate(-1, 1)
    assert ref.area() == 0
    assert ref.area(True) == dict()
    assert ref.get_bounding_box() is None
    assert ref.get_polygons() == []
    assert ref.get_polygons(True) == dict()
    assert ref.origin[0] == ref.origin[1] == 0


def test_notempty():
    name = 'ca_notempty'
    c = gdspy.Cell(name)
    ref = gdspy.CellArray(name, 2, 3, (3, 2), (1, -1), 90, 2, True)
    ref.translate(-1, 1)
    c.add(gdspy.Rectangle((0, 0), (1, 2), 2, 3))
    assert ref.area() == 48
    assert ref.area(True) == {(2, 3): 48}
    err = numpy.array(((0, 0), (8, 5))) - ref.get_bounding_box()
    assert numpy.max(numpy.abs(err)) < 1e-15
    assert ref.origin[0] == ref.origin[1] == 0
    r = gdspy.boolean([gdspy.Rectangle((0, 0), (8, 2)), gdspy.Rectangle((0, 3), (8, 5))], ref.get_polygons(), 'xor', 1e-6, 0)
    assert r is None
    d = ref.get_polygons(True)
    assert len(d.keys()) == 1
    r = gdspy.boolean([gdspy.Rectangle((0, 0), (8, 2)), gdspy.Rectangle((0, 3), (8, 5))], d[(2, 3)], 'xor', 1e-6, 0)
    assert r is None
