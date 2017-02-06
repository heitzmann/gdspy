######################################################################
#                                                                    #
#  Copyright 2009-2017 Lucas Heitzmann Gabrielli                     #
#                                                                    #
#  This file is part of gdspy.                                       #
#                                                                    #
#  gdspy is free software: you can redistribute it and/or modify it  #
#  under the terms of the GNU General Public License as published    #
#  by the Free Software Foundation, either version 3 of the          #
#  License, or any later version.                                    #
#                                                                    #
#  gdspy is distributed in the hope that it will be useful, but      #
#  WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     #
#  GNU General Public License for more details.                      #
#                                                                    #
#  You should have received a copy of the GNU General Public         #
#  License along with gdspy.  If not, see                            #
#  <http://www.gnu.org/licenses/>.                                   #
#                                                                    #
######################################################################

import gdspy
import numpy


def test_noreference():
    name = 'cr_noreference'
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
    r = gdspy.fast_boolean(ref.get_polygons(),
                           gdspy.Rectangle((0, 0), (4, 2)), 'xor', 1e-6, 0)
    assert r is None
    d = ref.get_polygons(True)
    assert len(d.keys()) == 1
    r = gdspy.fast_boolean(d[(2, 3)],
                           gdspy.Rectangle((0, 0), (4, 2)), 'xor', 1e-6, 0)
    assert r is None
