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


def test_hobby1(target):
    cell = gdspy.Cell('test', True)
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)])
    cell.add(gdspy.Polygon(c.get_points(), layer=1))
    c = gdspy.Curve(2, 0, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[numpy.pi / 3, None, None, None])
    cell.add(gdspy.Polygon(c.get_points(), layer=3))
    c = gdspy.Curve(4, 0, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[None, None, None, 2 / 3.0 * numpy.pi])
    cell.add(gdspy.Polygon(c.get_points(), layer=5))
    c = gdspy.Curve(0, 2, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[numpy.pi / 3, None, None, 3 / 4.0 * numpy.pi])
    cell.add(gdspy.Polygon(c.get_points(), layer=7))
    c = gdspy.Curve(2, 2, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[None, None, numpy.pi / 2, None])
    cell.add(gdspy.Polygon(c.get_points(), layer=9))
    c = gdspy.Curve(4, 2, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[None, 0, None, None])
    cell.add(gdspy.Polygon(c.get_points(), layer=11))
    c = gdspy.Curve(0, 4, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[None, 0, None, -numpy.pi / 2])
    cell.add(gdspy.Polygon(c.get_points(), layer=13))
    c = gdspy.Curve(2, 4, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[None, 0, -numpy.pi, -numpy.pi / 2])
    cell.add(gdspy.Polygon(c.get_points(), layer=15))
    c = gdspy.Curve(4, 4, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[-numpy.pi / 4, 0, numpy.pi / 2, -numpy.pi])
    cell.add(gdspy.Polygon(c.get_points(), layer=17))

    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=2))
    c = gdspy.Curve(2, 0, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[numpy.pi / 3, None, None, None], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=4))
    c = gdspy.Curve(4, 0, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[None, None, None, 2 / 3.0 * numpy.pi], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=6))
    c = gdspy.Curve(0, 2, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[numpy.pi / 3, None, None, 3 / 4.0 * numpy.pi], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=8))
    c = gdspy.Curve(2, 2, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[None, None, numpy.pi / 2, None], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=10))
    c = gdspy.Curve(4, 2, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[None, 0, None, None], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=12))
    c = gdspy.Curve(0, 4, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[None, 0, None, -numpy.pi / 2], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=14))
    c = gdspy.Curve(2, 4, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[None, 0, -numpy.pi, -numpy.pi / 2], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=16))
    c = gdspy.Curve(4, 4, tolerance=1e-3)
    c.i([(1, 0), (1, 1), (0, 1)], angles=[-numpy.pi / 4, 0, numpy.pi / 2, -numpy.pi], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=18))
    assertsame(target['Hobby1'], cell)

def test_hobby2(target):
    cell = gdspy.Cell('test', True)
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)])
    cell.add(gdspy.Polygon(c.get_points(), layer=1))
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], curl_start=0)
    cell.add(gdspy.Polygon(c.get_points(), layer=2))
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], curl_end=0)
    cell.add(gdspy.Polygon(c.get_points(), layer=3))
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], curl_start=0, curl_end=0)
    cell.add(gdspy.Polygon(c.get_points(), layer=4))
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], angles=[numpy.pi / 2, None, None, None, -numpy.pi / 2],
        curl_start=0, curl_end=0)
    cell.add(gdspy.Polygon(c.get_points(), layer=5))
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], angles=[None, 0, None, 0, None],
        curl_start=0, curl_end=1)
    cell.add(gdspy.Polygon(c.get_points(), layer=6))
    assertsame(target['Hobby2'], cell)

def test_hobby3(target):
    cell = gdspy.Cell('test', True)
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)])
    cell.add(gdspy.Polygon(c.get_points(), layer=1))
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], t_in=2)
    cell.add(gdspy.Polygon(c.get_points(), layer=2))
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], t_out=2)
    cell.add(gdspy.Polygon(c.get_points(), layer=3))
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], t_in=2, t_out=2)
    cell.add(gdspy.Polygon(c.get_points(), layer=4))
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], t_in=[2, 1, 1, 1, 1], t_out=[1, 1, 1, 1, 2])
    cell.add(gdspy.Polygon(c.get_points(), layer=5))
    c = gdspy.Curve(0, 0, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], t_in=[1, 1, 2, 1, 1], t_out=[1, 2, 1, 1, 1])
    cell.add(gdspy.Polygon(c.get_points(), layer=6))
    assertsame(target['Hobby3'], cell)

def test_hobby4(target):
    cell = gdspy.Cell('test', True)
    c = gdspy.Curve(0, 3, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=10))
    c = gdspy.Curve(0, 3, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], t_in=[2, 1, 1, 1, 1],
        t_out=[1, 1, 1, 1, 2], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=11))
    c = gdspy.Curve(0, 3, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)], t_in=[1, 1, 2, 1, 1],
        t_out=[1, 2, 1, 1, 1], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=12))
    c = gdspy.Curve(0, 3, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)],
        angles=[numpy.pi * 3 / 4.0, None, None, None, -numpy.pi * 3 / 4.0],
        t_in=[2, 1, 1, 1, 1], t_out=[1, 1, 1, 1, 2], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=13))
    c = gdspy.Curve(0, 3, tolerance=1e-3)
    c.i([(1, 2), (2, 1), (3, 2), (4, 0)],
        angles=[numpy.pi * 3 / 4.0, None, None, None, -numpy.pi * 3 / 4.0],
        t_in=[1, 1, 1, 1, 1], t_out=[1, 1, 1, 1, 1], cycle=True)
    cell.add(gdspy.Polygon(c.get_points(), layer=14))
    assertsame(target['Hobby4'], cell)

