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


def test_robustpath_warnings():
    with pytest.warns(UserWarning):
        gdspy.RobustPath((0, 0), 1, ends='smooth', gdsii_path=True)
    with pytest.warns(UserWarning):
        gdspy.RobustPath((0, 0), [1, 1], ends=['flush', 'smooth'], gdsii_path=True)

def test_robustpath_len():
    rp = gdspy.RobustPath((0, 0), [0.1, 0.2, 0.1], 0.15, layer=[1, 2, 3])
    assert(len(rp) == 0)
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
    assert(len(rp) == 10)

def test_robustpath_call():
    rp = gdspy.RobustPath((0, 0), [0.1, 0.2, 0.1], 0.15, layer=[1, 2, 3])
    rp.segment((1, 0))
    rp.turn(5, 'll')
    rp.segment((-1, 0), relative=True)
    assert(numpy.sum((rp(0) - numpy.array([(0, 0.15), (0, 0), (0, -0.15)]))**2) < 1e-12)
    assert(numpy.sum((rp(0.5) - numpy.array([(0.5, 0.15), (0.5, 0), (0.5, -0.15)]))**2) < 1e-12)
    assert(numpy.sum((rp(1) - numpy.array([(1, 0.15), (1, 0), (1, -0.15)]))**2) < 1e-12)
    assert(numpy.sum((rp(1.5, arm=1) - numpy.array([(5.9, 5), (6.1, 5), (6.2, 5)]))**2) < 1e-12)
    assert(numpy.sum((rp(1.5, arm=-1) - numpy.array([(5.8, 5), (5.9, 5), (6.1, 5)]))**2) < 1e-12)
    assert(numpy.sum((rp(2) - numpy.array([(1, 9.85), (1, 10), (1, 10.15)]))**2) < 1e-12)
    assert(numpy.sum((rp(3) - numpy.array([(0, 9.85), (0, 10), (0, 10.15)]))**2) < 1e-12)

def test_robustpath_grad():
    rp = gdspy.RobustPath((0, 0), 0.1)
    rp.segment((1, 0), 0.3)
    rp.segment((1, 1), 0.1)
    assert(numpy.sum((numpy.array((1, 0)) - rp.grad(0))**2) < 1e-12)
    assert(numpy.sum((numpy.array((1, 0)) - rp.grad(1))**2) < 1e-12)
    assert(numpy.sum((numpy.array((0, 1)) - rp.grad(1, side='+'))**2) < 1e-12)
    assert(numpy.sum((numpy.array((0, 1)) - rp.grad(2))**2) < 1e-12)
    assert(numpy.sum((numpy.array((0.1, 1)) - rp.grad(2, arm=-1))**2) < 1e-12)
    assert(numpy.sum((numpy.array((-0.1, 1)) - rp.grad(2, arm=1))**2) < 1e-12)

def test_robustpath_width():
    rp = gdspy.RobustPath((0, 0), [0.1, 0.3])
    rp.segment((1, 0), 0.3)
    rp.segment((1, 1), [0.1, 0.2])
    assert(numpy.sum((numpy.array((0.1, 0.3)) - rp.width(0))**2) < 1e-12)
    assert(numpy.sum((numpy.array((0.2, 0.3)) - rp.width(0.5))**2) < 1e-12)
    assert(numpy.sum((numpy.array((0.3, 0.3)) - rp.width(1))**2) < 1e-12)
    assert(numpy.sum((numpy.array((0.2, 0.25)) - rp.width(1.5))**2) < 1e-12)
    assert(numpy.sum((numpy.array((0.1, 0.2)) - rp.width(2))**2) < 1e-12)

def test_robustpath1(target):
    cell = gdspy.Cell('test', True)
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
    assertsame(target['RobustPath1'], cell)

def test_robustpath2(target):
    cell = gdspy.Cell('test', True)
    rp = gdspy.RobustPath((0, 0), [0.1, 0.2, 0.1], 0.15, layer=[1, 2, 3])
    assert(len(rp) == 0)
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
    assert(len(rp) == 10)
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
    assertsame(target['RobustPath2'], cell)

def test_robustpath3(target):
    cell = gdspy.Cell('test', True)
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
    assertsame(target['RobustPath3'], cell)

def test_robustpath_gdsiipath():
    cells = []
    for gdsii_path in [True, False]:
        cells.append(gdspy.Cell(str(gdsii_path), True))
        rp = gdspy.RobustPath((0, 0), 0.05, [-0.1, 0.1], ends=['extended', (0.1, 0.2)],
                              layer=[0, 1], gdsii_path=gdsii_path)
        rp.segment((1, 1))
        rp.parametric(lambda u: numpy.array((u, u - u**2)))
        cells[-1].add(rp)
    assertsame(*cells)

def test_robustpath_getpolygons():
    rp = gdspy.RobustPath((0, 0),
                          0.05, [-0.1, 0.1], ends=['extended', (0.1, 0.2)],
                          layer=[0, 1], datatype=[1, 0])
    rp.segment((1, 1))
    rp.parametric(lambda u: numpy.array((u, u - u**2)))
    d = rp.get_polygons(True)
    l = rp.get_polygons()
    assert len(d) == 2
    assert (1, 0) in d
    assert (0, 1) in d
    assert sum(len(p) for p in d.values()) == len(l)
    assert sum(sum(len(x) for x in p) for p in d.values()) == sum(len(x) for x in l)
    ps = rp.to_polygonset()
    assert len(ps.layers) == len(ps.datatypes) == len(ps.polygons)
    assert gdspy.RobustPath((0, 0), 1).to_polygonset() == None

def test_robustpath_togds(tmpdir):
    cell = gdspy.Cell('robustpath')
    rp = gdspy.RobustPath((0, 0), 0.1, layer=[1])
    rp.segment((1, 1))
    rp.segment((2, 3), 0)
    cell.add(rp)
    rp = gdspy.RobustPath((2, 0), [0.1, 0.2], 0.2, ends=['round', (0.2, 0.1)],
                          layer=4, datatype=[1, 1], gdsii_path=True)
    rp.segment((3, 1))
    cell.add(rp)
    rp = gdspy.RobustPath((0, 0), 0.1, layer=5, tolerance=1e-5, max_points=0,
                          max_evals=1e6, gdsii_path=True)
    rp.segment((10, 0))
    rp.turn(20, 'll')
    rp.turn(20, 'rr')
    rp.turn(20, 'll')
    cell.add(rp)
    fname = str(tmpdir.join('test.gds'))
    with pytest.warns(UserWarning):
        gdspy.write_gds(fname, unit=1, precision=1e-7)
    lib = gdspy.GdsLibrary(infile=fname, rename={'robustpath': 'file'})
    assertsame(lib.cell_dict['file'], cell, tolerance=1e-3)
    gdspy.current_library = gdspy.GdsLibrary()
