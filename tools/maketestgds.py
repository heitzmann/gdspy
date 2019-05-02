######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import gdspy
import numpy


### PolygonSet

cell = gdspy.Cell('PolygonSet')

p = gdspy.PolygonSet([
    [(10, 0), (11, 0), (10, 1)],
    [(11, 0), (10, 1), (11, 1)],
    [(11, 1), (12, 1), (11, 2)],
], 1, 2)
cell.add(p)

cell = gdspy.Cell('PolygonSet_fillet')

orig = gdspy.PolygonSet([
    [(0, 0), (-1, 0), (0, -1), (0.5, -0.5), (1, 0), (1, 1), (4, -1),
     (1, 3), (1, 2), (0, 1)],
    [(2, -1), (3, -1), (2.5, -2)],
])
orig.datatypes=[0, 1]
p = gdspy.copy(orig, 0, 5)
p.layers = [1, 1]
p.fillet(0.3, max_points=0)
cell.add(p)
p = gdspy.copy(orig, 5, 5)
p.layers = [2, 2]
p.fillet([0.3, 0.2, 0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.1, 0.2, 0],
         max_points=0)
cell.add(p)
p = gdspy.copy(orig, 5, 0)
p.layers = [3, 3]
p.fillet([[0.1, 0.1, 0.4, 0, 0.4, 0.1, 0.1, 0.4, 0.4, 0.1], [0.2, 0.2, 0.5]],
         max_points=0)
cell.add(p)
p = gdspy.copy(orig, 0, 0)
p.layers = [4, 4]
p.fillet([0.8, [10.0, 10.0, 20.0]], max_points=199, precision=1e-6)
cell.add(p)


### FlexPath

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


cell = gdspy.Cell('FlexPath1')

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

cell = gdspy.Cell('FlexPath2')

fp = gdspy.FlexPath([(0, 0), (0.5, 0), (1, 0), (1, 1), (0, 1), (-1, -2), (-2, 0)], 0.05,
                    [0, -0.1, 0, 0.1], corners=['natural', 'circular bend',
                                                'circular bend', 'circular bend'],
                    ends=['flush', 'extended', (0.1, 0.2), 'round'], tolerance=1e-4,
                    layer=[0, 1, 1, 2], bend_radius=[0, 0.3, 0.3, 0.2], max_points=10)
cell.add(fp)

cell = gdspy.Cell('FlexPath3')

pts = numpy.array([(0, 0), (0.5, 0), (1, 0), (1, 2), (3, 0), (2, -1), (2, -2),
                   (0, -1), (1, -2), (1, -3)])
fp = gdspy.FlexPath(pts + numpy.array((0, 5)), [0.1, 0.1, 0.1], 0.15,
                    layer=[1, 2, 3], corners=['natural', 'miter', 'bevel'],
                    ends=(0.5, 0))
cell.add(fp)
fp = gdspy.FlexPath(pts + numpy.array((5, 0)), [0.1, 0.1, 0.1], 0.15,
                    layer=[4, 5, 6], corners=['round', 'smooth', broken],
                    ends=[pointy, 'smooth', (0, 0.5)])
cell.add(fp)

cell = gdspy.Cell('FlexPath4')

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


### RobustPath

cell = gdspy.Cell('RobustPath1')

rp = gdspy.RobustPath((0, 0), 0.1, layer=[1], gdsii_path=True)
rp.segment((1, 1))
cell.add(rp)
rp = gdspy.RobustPath((1, 0), 0.1, [-0.1, 0.1], tolerance=1e-5,
                      ends=['round', 'extended'], layer=[2, 3], max_points=6)
rp.segment((2, 1))
cell.add(rp)
rp = gdspy.RobustPath((2, 0), [0.1, 0.2], 0.2, ends=(0.2, 0.1),
                      layer=4, datatype=[1, 1])
rp.segment((3, 1))
cell.add(rp)
rp = gdspy.RobustPath((3, 0), [0.1, 0.2, 0.1], [-0.2, 0, 0.2],
                      ends=[(0.2, 0.1), 'smooth', 'flush'], datatype=5)
rp.segment((4, 1))
cell.add(rp)

cell = gdspy.Cell('RobustPath2')

rp = gdspy.RobustPath((0, 0), [0.1, 0.2, 0.1], 0.15, layer=[1, 2, 3])
rp.segment((1, 0))
rp.segment((1, 1), 0.1, 0.05)
rp.segment((1, 1), [0.2, 0.1, 0.1], -0.05, True)
rp.segment((-1, 1), 0.2, [-0.2, 0, 0.3], True)
rp.arc(2, 0, 0.5 * numpy.pi)
rp.arc(3, 0.7 * numpy.pi, numpy.pi, 0.1, 0)
rp.arc(2, 0.4 * numpy.pi, -0.4 * numpy.pi, [0.1, 0.2, 0.1], [0.2, 0, -0.2])
rp.turn(1, -0.3 * numpy.pi)
rp.turn(1, 'rr', 0.15)
rp.turn(0.5, 'l', [0.05, 0.1, 0.05], [0.15, 0, -0.15])
cell.add(rp)
rp = gdspy.RobustPath((-5, 6), 0.8, layer=20, ends='round', tolerance=1e-4)
rp.segment((1, 1), 0.1, relative=True)
cell.add(rp)
rp = gdspy.RobustPath((-5, 6), 0.8, layer=21, ends='extended', tolerance=1e-4)
rp.segment((1, 1), 0.1, relative=True)
cell.add(rp)
rp = gdspy.RobustPath((-5, 6), 0.8, layer=22, ends=(0.1, 0.2), tolerance=1e-4)
rp.segment((1, 1), 0.1, relative=True)
cell.add(rp)
rp = gdspy.RobustPath((-5, 6), 0.8, layer=23, ends='smooth', tolerance=1e-4)
rp.segment((1, 1), 0.1, relative=True)
cell.add(rp)
rp = gdspy.RobustPath((-3, 6), 0.8, layer=10, ends='round', tolerance=1e-5)
rp.segment((1, 0), 0.1, relative=True)
rp.segment((0, 1), 0.8, relative=True)
cell.add(rp)
rp = gdspy.RobustPath((-3, 6), 0.8, layer=11, ends='extended', tolerance=1e-5)
rp.segment((1, 0), 0.1, relative=True)
rp.segment((0, 1), 0.8, relative=True)
cell.add(rp)
rp = gdspy.RobustPath((-3, 6), 0.8, layer=12, ends='smooth', tolerance=1e-5)
rp.segment((1, 0), 0.1, relative=True)
rp.segment((0, 1), 0.8, relative=True)
cell.add(rp)
rp = gdspy.RobustPath((-3, 8), 0.1, layer=13, ends='round', tolerance=1e-5)
rp.segment((1, 0), 0.8, relative=True)
rp.segment((0, 1), 0.1, relative=True)
cell.add(rp)
rp = gdspy.RobustPath((-3, 8), 0.1, layer=14, ends=(0.2, 0.2), tolerance=1e-5)
rp.segment((1, 0), 0.8, relative=True)
rp.segment((0, 1), 0.1, relative=True)
cell.add(rp)
rp = gdspy.RobustPath((-3, 8), 0.1, layer=15, ends='smooth', tolerance=1e-5)
rp.segment((1, 0), 0.8, relative=True)
rp.segment((0, 1), 0.1, relative=True)
cell.add(rp)
rp = gdspy.RobustPath((5, 2), [0.05, 0.1, 0.2], [-0.2, 0, 0.4], layer=[4, 5, 6])
rp.parametric(lambda u: numpy.array((5.5 + 3 * u, 2 + 3 * u**2)), relative=False)
rp.segment((0, 1), relative=True)
rp.parametric(lambda u: numpy.array((2 * numpy.cos(0.5 * numpy.pi * u) - 2,
                                     3 * numpy.sin(0.5 * numpy.pi * u))),
              width=[0.2, 0.1, 0.05], offset=[-0.3, 0, 0.3])
rp.parametric(lambda u: numpy.array((-2*u, 0)), width=0.1, offset=0.2)
rp.bezier([(-3, 0), (-2, -3), (0, -4), (0, -5)], offset=[-0.2, 0, 0.2])
rp.bezier([(4.5, 0), (1, -1),  (1, 5), (3, 2), (5, 2)], width=[0.05, 0.1, 0.2],
          offset=[-0.2, 0, 0.4], relative=False)
cell.add(rp)
rp = gdspy.RobustPath((2, -1), 0.1, layer=7, tolerance=1e-4, max_points=0)
rp.smooth([(1, 0), (1, -1), (0, -1)], angles=[numpy.pi / 3, None, -2 / 3.0 * numpy.pi, None],
          cycle=True)
cell.add(rp)
rp = gdspy.RobustPath((2.5, -1.5), 0.1, layer=8)
rp.smooth([(3, -1.5), (4, -2), (5, -1), (6, -2), (7, -1.5), (7.5, -1.5)], relative=False,
          width=0.2)
cell.add(rp)

cell = gdspy.Cell('RobustPath3')

rp = gdspy.RobustPath((0, 0), 0.1)
rp.parametric(lambda u: numpy.array((3 * numpy.sin(numpy.pi * u),
                                     -3 * numpy.cos(numpy.pi * u))),
              relative=False)
rp.parametric(lambda u: numpy.array((3.5 - 3 * numpy.cos(numpy.pi * u),
                                     -0.5 + 3 * numpy.sin(numpy.pi * u))),
              lambda u: numpy.array((numpy.sin(numpy.pi * u),
                                     numpy.cos(numpy.pi * u))),
              relative=True)
cell.add(rp)


### Curve

cell = gdspy.Cell('Hobby1')

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

cell = gdspy.Cell('Hobby2')

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

cell = gdspy.Cell('Hobby3')

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

cell = gdspy.Cell('Hobby4')

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


### END

gdspy.write_gds('tests/test.gds', unit=1, precision=1e-7)
gdspy.LayoutViewer(cells=[cell])
