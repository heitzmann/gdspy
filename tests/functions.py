######################################################################
#                                                                    #
#  Copyright 2009-2018 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import numpy
import gdspy


def test_8b_f():
    f = gdspy._eight_byte_real_to_float
    assert f(b'\x00\x00\x00\x00\x00\x00\x00\x00') == 0
    assert f(b'\x41\x10\x00\x00\x00\x00\x00\x00') == 1
    assert f(b'\x41\x20\x00\x00\x00\x00\x00\x00') == 2
    assert f(b'\xC1\x30\x00\x00\x00\x00\x00\x00') == -3


def test_f_8b():
    g = gdspy._eight_byte_real
    assert b'\x00\x00\x00\x00\x00\x00\x00\x00' == g(0)
    assert b'\x41\x10\x00\x00\x00\x00\x00\x00' == g(1)
    assert b'\x41\x20\x00\x00\x00\x00\x00\x00' == g(2)
    assert b'\xC1\x30\x00\x00\x00\x00\x00\x00' == g(-3)


def test_twoway():
    f = gdspy._eight_byte_real_to_float
    g = gdspy._eight_byte_real
    for x in [0, 1.5, -numpy.pi, 1 / 3.0e12, -1.0e12 / 7, 1.1e75, -0.9e-78]:
        assert x == f(g(x))
    for _ in range(10000):
        x = 10**(numpy.random.random() * 150 - 75)
        assert x == f(g(x))


def test_inside():
    polygons = [
        gdspy.Round((0, 0), 10, inner_radius=5, number_of_points=180),
        gdspy.Rectangle((20, -10), (40, 10)),
        gdspy.Rectangle((-10, 0), (10, 20))
    ]
    assert gdspy.inside([(0, 0)], polygons) == (True, )
    assert gdspy.inside([(0, 0), (0, 30), (30, 0), (0, -1)], polygons) == \
        (True, False, True, False)
    assert gdspy.inside([[(0, 0), (0, 30), (30, 0),
                          (0, -1)], [(0, -1), (0, 30)], [(0, 0), (30, 0)]],
                        polygons, 'any') == (True, False, True)
    assert gdspy.inside([[(0, 0), (0, 30), (30, 0),
                          (0, -1)], [(0, -1), (0, 30)], [(0, 0), (30, 0)]],
                        polygons, 'all') == (False, False, True)


def test_copy():
    p = gdspy.Rectangle((0, 0), (1, 1))
    q = gdspy.copy(p, 1, -1)
    assert set(p.polygons[0][:, 0]) == {0, 1}
    assert set(p.polygons[0][:, 1]) == {0, 1}
    assert set(q.polygons[0][:, 0]) == {1, 2}
    assert set(q.polygons[0][:, 1]) == {-1, 0}
    p = gdspy.PolygonSet([[(0, 0), (1, 0), (0, 1)], [(2, 2), (3, 2), (2, 3)]])
    q = gdspy.copy(p, 1, -1)
    assert set(p.polygons[0][:, 0]) == {0, 1}
    assert set(p.polygons[0][:, 1]) == {0, 1}
    assert set(q.polygons[0][:, 0]) == {1, 2}
    assert set(q.polygons[0][:, 1]) == {-1, 0}
    assert set(p.polygons[1][:, 0]) == {2, 3}
    assert set(p.polygons[1][:, 1]) == {2, 3}
    assert set(q.polygons[1][:, 0]) == {3, 4}
    assert set(q.polygons[1][:, 1]) == {1, 2}
    l = gdspy.Label('text', (0, 1))
    m = gdspy.copy(l, -1, 1)
    assert l.position[0] == 0 and l.position[1] == 1
    assert m.position[0] == -1 and m.position[1] == 2
    c = gdspy.CellReference('empty', (0, 1), ignore_missing=True)
    d = gdspy.copy(c, -1, 1)
    assert c.origin == (0, 1)
    assert d.origin == (-1, 2)
    c = gdspy.CellArray('empty', 2, 3, (1, 0), (0, 1), ignore_missing=True)
    d = gdspy.copy(c, -1, 1)
    assert c.origin == (0, 1)
    assert d.origin == (-1, 2)


def test_write_gds(tmpdir):
    gdspy.current_library = gdspy.GdsLibrary()
    c1 = gdspy.Cell('fu_rw_gds_1')
    c1.add(gdspy.Rectangle((0, -1), (1, 2), 2, 4))
    c1.add(gdspy.Label('label', (1, -1), 'w', 45, 1.5, True, 5, 6))
    c2 = gdspy.Cell('fu_rw_gds_2')
    c2.add(gdspy.Round((0, 0), 1, number_of_points=32, max_points=20))
    c3 = gdspy.Cell('fu_rw_gds_3')
    c3.add(gdspy.CellReference(c1, (0, 1), -90, 2, True))
    c4 = gdspy.Cell('fu_rw_gds_4')
    c4.add(gdspy.CellArray(c2, 2, 3, (1, 4), (-1, -2), 180, 0.5, True))

    fname1 = str(tmpdir.join('test1.gds'))
    gdspy.write_gds(fname1, name='lib', unit=2e-6, precision=1e-8)
    lib1 = gdspy.GdsLibrary(
        infile=fname1,
        units='convert',
        rename={'fu_rw_gds_1': '1'},
        layers={2: 4},
        datatypes={4: 2},
        texttypes={6: 7})
    assert lib1.name == 'lib'
    assert len(lib1.cell_dict) == 4
    assert set(lib1.cell_dict.keys()) == {
        '1', 'fu_rw_gds_2', 'fu_rw_gds_3', 'fu_rw_gds_4'
    }
    c = lib1.cell_dict['1']
    assert len(c.elements) == len(c.labels) == 1
    assert c.elements[0].area() == 12.0
    assert c.elements[0].layers == [4]
    assert c.elements[0].datatypes == [2]
    assert c.labels[0].text == 'label'
    assert c.labels[0].position[0] == 2 and c.labels[0].position[1] == -2
    assert c.labels[0].anchor == 4
    assert c.labels[0].rotation == 45
    assert c.labels[0].magnification == 1.5
    assert c.labels[0].x_reflection == True
    assert c.labels[0].layer == 5
    assert c.labels[0].texttype == 7

    c = lib1.cell_dict['fu_rw_gds_2']
    assert len(c.elements) == 2
    assert isinstance(c.elements[0], gdspy.Polygon) \
           and isinstance(c.elements[1], gdspy.Polygon)

    c = lib1.cell_dict['fu_rw_gds_3']
    assert len(c.elements) == 1
    assert isinstance(c.elements[0], gdspy.CellReference)
    assert c.elements[0].ref_cell == lib1.cell_dict['1']
    assert c.elements[0].origin[0] == 0 and c.elements[0].origin[1] == 2
    assert c.elements[0].rotation == -90
    assert c.elements[0].magnification == 2
    assert c.elements[0].x_reflection == True

    c = lib1.cell_dict['fu_rw_gds_4']
    assert len(c.elements) == 1
    assert isinstance(c.elements[0], gdspy.CellArray)
    assert c.elements[0].ref_cell == lib1.cell_dict['fu_rw_gds_2']
    assert c.elements[0].origin[0] == -2 and c.elements[0].origin[1] == -4
    assert c.elements[0].rotation == 180
    assert c.elements[0].magnification == 0.5
    assert c.elements[0].x_reflection == True
    assert c.elements[0].spacing[0] == 2 and c.elements[0].spacing[1] == 8
    assert c.elements[0].columns == 2
    assert c.elements[0].rows == 3

    fname2 = str(tmpdir.join('test2.gds'))
    with open(fname2, 'wb') as fout:
        gdspy.write_gds(fout, name='lib2', unit=2e-3, precision=1e-5)
    with open(fname2, 'rb') as fin:
        lib2 = gdspy.GdsLibrary()
        lib2.read_gds(fin)
    assert lib2.name == 'lib2'
    assert len(lib2.cell_dict) == 4
