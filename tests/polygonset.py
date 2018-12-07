######################################################################
#                                                                    #
#  Copyright 2009-2018 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import pytest
import gdspy


def test_fillet():
    r = gdspy.Rectangle((0, 0), (1, 1))
    r.fillet(0.25, points_per_2pi=64)
    #gdspy.Cell('1').add(r)
    assert len(r.polygons) == 1
    assert r.polygons[0].shape[0] == 64
    r = gdspy.Rectangle((0, 0), (1, 1))
    r.fillet([[0.1, 0.2, 0.3, 0.4]], points_per_2pi=64)
    #gdspy.Cell('2').add(r)
    assert len(r.polygons) == 1
    assert r.polygons[0].shape[0] == 64
    r = gdspy.Rectangle((0, 0), (1, 1))
    r.fillet([0.4, 0.3, 0.2, 0.1], points_per_2pi=64)
    #gdspy.Cell('3').add(r)
    assert len(r.polygons) == 1
    assert r.polygons[0].shape[0] == 64
    r = gdspy.PolygonSet([[(0, 0), (0, 1), (0.5, 1), (0.5, 0)], [(0.5, 0), (0.5, 1), (1, 1), (1, 0)]])
    r.fillet([0.1, [0.1, 0.3, 0.1, 0.3]], points_per_2pi=64)
    #gdspy.Cell('4').add(r)
    assert len(r.polygons) == 2
    assert r.polygons[0].shape[0] == 64
    assert r.polygons[1].shape[0] == 64
    r = gdspy.Rectangle((0, 0), (1, 1))
    r = gdspy.PolygonSet([[(0, 0), (0, 1), (0.5, 1), (0.5, 0)], [(0.5, 0), (0.5, 1), (1, 1), (1, 0)]])
    r.fillet([0.1, 0.3, 0.1, 0.3, 0.1, 0.3, 0.1, 0.3], points_per_2pi=64)
    #gdspy.Cell('5').add(r)
    assert len(r.polygons) == 2
    assert r.polygons[0].shape[0] == 64
    assert r.polygons[1].shape[0] == 64
    #gdspy.LayoutViewer()
