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
    name = 'cr_noreference'
    with pytest.warns(UserWarning):
        ref = gdspy.CellReference(name, (1, -1), 90, 2, True)
    ref.translate(-1, 1)
    assert ref.ref_cell == name
    assert ref.area() == 0
    assert ref.area(True) == dict()
    assert ref.get_bounding_box() is None
    assert ref.get_polygons() == []
    assert ref.get_polygons(True) == dict()
    assert ref.origin[0] == ref.origin[1] == 0


def test_empty():
    name = 'cr_empty'
    c = gdspy.Cell(name)
    ref = gdspy.CellReference(name, (1, -1), 90, 2, True)
    ref.translate(-1, 1)
    assert ref.area() == 0
    assert ref.area(True) == dict()
    assert ref.get_bounding_box() is None
    assert ref.get_polygons() == []
    assert ref.get_polygons(True) == dict()
    assert ref.origin[0] == ref.origin[1] == 0


def test_notempty():
    name = 'cr_notempty'
    c = gdspy.Cell(name)
    ref = gdspy.CellReference(name, (1, -1), 90, 2, True)
    ref.translate(-1, 1)
    c.add(gdspy.Rectangle((0, 0), (1, 2), 2, 3))
    assert ref.area() == 8
    assert ref.area(True) == {(2, 3): 8}
    err = numpy.array(((0, 0), (4, 2))) - ref.get_bounding_box()
    assert numpy.max(numpy.abs(err)) < 1e-15
    assert ref.origin[0] == ref.origin[1] == 0
    r = gdspy.boolean(ref.get_polygons(), gdspy.Rectangle((0, 0), (4, 2)), 'xor', 1e-6, 0)
    assert r is None
    d = ref.get_polygons(True)
    assert len(d.keys()) == 1
    r = gdspy.boolean(d[(2, 3)], gdspy.Rectangle((0, 0), (4, 2)), 'xor', 1e-6, 0)
    assert r is None
