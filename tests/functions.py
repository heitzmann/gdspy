######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import pytest
import datetime
import numpy
import gdspy


def equals(x, y):
    return gdspy.boolean(x, y, 'xor') is None


@pytest.fixture()
def library():
    lib = gdspy.GdsLibrary()
    c1 = gdspy.Cell('cell1', True)
    c1.add(gdspy.Rectangle((0, -1), (1, 2), 2, 4))
    c1.add(gdspy.Label('label', (1, -1), 'w', 45, 1.5, True, 5, 6))
    c2 = gdspy.Cell('cell2', True)
    c2.add(gdspy.Round((0, 0), 1, number_of_points=32, max_points=20))
    c3 = gdspy.Cell('cell3', True)
    c3.add(gdspy.CellReference(c1, (0, 1), -90, 2, True))
    c4 = gdspy.Cell('cell04', True)
    c4.add(gdspy.CellArray(c2, 2, 3, (1, 4), (-1, -2), 180, 0.5, True))
    lib.add([c1, c2, c3, c4])
    return lib


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


def test_gather():
    def same_points(x, y):
        for px, py in zip(x, y):
            for ptx, pty in zip(px, py):
                for cx, cy in zip(ptx, pty):
                    if cx != cy: return False
        return True

    gdspy.current_library = gdspy.GdsLibrary()
    pts = [(0, 0), (1, 1), (1, 0)]
    ps1 = gdspy.Round((10, 10), 1, inner_radius=0.2)
    ps2 = gdspy.Path(0.1, (-1, -1), 2, 1).segment(2, '-x')
    c = gdspy.Cell('C1').add(gdspy.Rectangle((-4, 3), (-5, 4)))
    cr = gdspy.CellReference(c, (10, -10))
    ca = gdspy.CellArray(c, 2, 1, (2, 0))
    assert gdspy._gather_polys(None) == []
    assert same_points(gdspy._gather_polys([pts]), [pts])
    assert same_points(gdspy._gather_polys(ps1), ps1.polygons)
    assert same_points(gdspy._gather_polys(ps2), ps2.polygons)
    assert same_points(gdspy._gather_polys(cr), cr.get_polygons())
    assert same_points(gdspy._gather_polys(ca), ca.get_polygons())
    result = [pts]
    result.extend(ps2.polygons)
    result.extend(cr.get_polygons())
    assert same_points(gdspy._gather_polys([pts, ps2, cr]), result)


def test_slice():
    poly = gdspy.Path(1, (1, 0), 2, 3).segment(2, '-x')
    left = gdspy.Path(1, (0, 0), 2, 3).segment(1, '-x')
    right = gdspy.Path(1, (1, 0), 2, 3).segment(1, '-x')
    result = gdspy.slice(poly, 0, 0)
    assert equals(result[0], left)
    assert equals(result[1], right)
    bot = gdspy.Path(1, (1, -1.5)).segment(2, '-x')
    top = gdspy.Path(1, (1, 1.5)).segment(2, '-x')
    result = gdspy.slice(poly, [0.1, -0.1], 1)
    assert equals(result[0], bot)
    assert equals(result[2], top)
    assert result[1] is None


def test_offset():
    gdspy.current_library = gdspy.GdsLibrary()
    r = gdspy.Rectangle((0, 0), (1, 2))
    result = gdspy.Rectangle((-1, -1), (2, 3))
    assert equals(gdspy.offset(r, 1), result)
    c = gdspy.Cell('OFFSET').add(r)
    ca = gdspy.CellArray(c, 2, 1, (1, 0))
    result = gdspy.Rectangle((0.2, 0.2), (1.8, 1.8))
    assert equals(gdspy.offset([ca], -0.2, join_first=True), result)
    v = [gdspy.Rectangle((-1, -1), (1, 1)), [(0, 0), (1, 0), (1, 1), (0, 1)]]
    x = 1 + 0.1 * numpy.tan(numpy.pi / 8)
    result = gdspy.Polygon([(-1.1, -x), (-1.1, x), (-x, 1.1), (x, 1.1), (1.1, x), (1.1, -x),
                            (x, -1.1), (-x, -1.1)], layer=8)
    assert equals(gdspy.offset(v, 0.1, join='bevel', layer=12), result)


def test_boolean():
    op1 = gdspy.Rectangle((0, 0), (3, 3))
    op2 = gdspy.Rectangle((1, 1), (2, 2))
    result = [[(0, 0), (3, 0), (3, 3), (0, 3), (0, 0), (1, 1), (1, 2), (2, 2), (2, 1), (1, 1)]]
    assert equals(gdspy.boolean(op1, op2, 'not'), result)
    op3 = gdspy.Rectangle((0, 0), (2, 2))
    assert equals(gdspy.boolean([op2, op3], None, 'or'), op3)


def test_inside():
    gdspy.current_library = gdspy.GdsLibrary()
    polygons = [
        gdspy.Round((0, 0), 10, inner_radius=5, number_of_points=180),
        gdspy.Rectangle((20, -10), (40, 10)).polygons[0],
        gdspy.CellReference(gdspy.Cell('X').add(gdspy.Rectangle((-10, 0), (10, 20)))),
    ]
    assert gdspy.inside([(0, 0)], polygons[0]) == (False, )
    assert gdspy.inside([(0, 0)], polygons[2]) == (True, )
    assert gdspy.inside([(0, 0)], polygons) == (True, )
    assert gdspy.inside([(0, 0), (0, 30), (30, 0), (0, -1)], polygons) == \
        (True, False, True, False)
    assert gdspy.inside([[(0, 0), (0, 30), (30, 0), (0, -1)], [(0, -1), (0, 30)], [(0, 0), (30, 0)]],
                        polygons, 'any') == (True, False, True)
    assert gdspy.inside([[(0, 0), (0, 30), (30, 0), (0, -1)], [(0, -1), (0, 30)], [(0, 0), (30, 0)]],
                        polygons, 'all') == (False, False, True)


def test_copy():
    gdspy.current_library = gdspy.GdsLibrary()
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


def test_write_gds(library, tmpdir):
    gdspy.current_library = library
    fname1 = str(tmpdir.join('test1.gds'))
    gdspy.write_gds(fname1, name='lib', unit=2e-6, precision=1e-8)
    lib1 = gdspy.GdsLibrary(infile=fname1, units='convert', rename={'cell1': '1'}, layers={2: 4},
                            datatypes={4: 2}, texttypes={6: 7})
    assert lib1.name == 'lib'
    assert len(lib1.cell_dict) == 4
    assert set(lib1.cell_dict.keys()) == {'1', 'cell2', 'cell3', 'cell04'}
    c = lib1.cell_dict['1']
    assert len(c.polygons) == len(c.labels) == 1
    assert c.polygons[0].area() == 12.0
    assert c.polygons[0].layers == [4]
    assert c.polygons[0].datatypes == [2]
    assert c.labels[0].text == 'label'
    assert c.labels[0].position[0] == 2 and c.labels[0].position[1] == -2
    assert c.labels[0].anchor == 4
    assert c.labels[0].rotation == 45
    assert c.labels[0].magnification == 1.5
    assert c.labels[0].x_reflection == True
    assert c.labels[0].layer == 5
    assert c.labels[0].texttype == 7

    c = lib1.cell_dict['cell2']
    assert len(c.polygons) == 2
    assert isinstance(c.polygons[0], gdspy.Polygon) and isinstance(c.polygons[1], gdspy.Polygon)

    c = lib1.cell_dict['cell3']
    assert len(c.references) == 1
    assert c.references[0].ref_cell == lib1.cell_dict['1']
    assert c.references[0].origin[0] == 0 and c.references[0].origin[1] == 2
    assert c.references[0].rotation == -90
    assert c.references[0].magnification == 2
    assert c.references[0].x_reflection == True

    c = lib1.cell_dict['cell04']
    assert len(c.references) == 1
    assert c.references[0].ref_cell == lib1.cell_dict['cell2']
    assert c.references[0].origin[0] == -2 and c.references[0].origin[1] == -4
    assert c.references[0].rotation == 180
    assert c.references[0].magnification == 0.5
    assert c.references[0].x_reflection == True
    assert c.references[0].spacing[0] == 2 and c.references[0].spacing[1] == 8
    assert c.references[0].columns == 2
    assert c.references[0].rows == 3

    fname2 = str(tmpdir.join('test2.gds'))
    with open(fname2, 'wb') as fout:
        gdspy.write_gds(fout, name='lib2', unit=2e-3, precision=1e-5)
    with open(fname2, 'rb') as fin:
        lib2 = gdspy.GdsLibrary()
        lib2.read_gds(fin)
    assert lib2.name == 'lib2'
    assert len(lib2.cell_dict) == 4


def test_gdsii_hash(library, tmpdir):
    out1 = str(tmpdir.join('test1.gds'))
    out2 = str(tmpdir.join('test2.gds'))
    library.write_gds(out1)
    library.write_gds(out2, timestamp=datetime.datetime.today() + datetime.timedelta(1))
    assert gdspy.gdsii_hash(out1) == gdspy.gdsii_hash(out2)


def test_get_gds_units(tmpdir):
    out = str(tmpdir.join('test1.gds'))
    lib = gdspy.GdsLibrary(unit=10.0, precision=0.1)
    lib.write_gds(out)
    assert (10.0, 0.1) == gdspy.get_gds_units(out)
    lib.unit = 0.2
    lib.precision = 5e-5
    out = str(tmpdir.join('test2.gds'))
    lib.write_gds(out)
    with open(out, 'rb') as fin:
        assert (0.2, 5e-5) == gdspy.get_gds_units(fin)


def test_get_binary_cells(library, tmpdir):
    out = str(tmpdir.join('test.gds'))
    now = datetime.datetime.today()
    library.write_gds(out, timestamp=now)
    bincells = gdspy.get_binary_cells(out)
    for name, cell in library.cell_dict.items():
        bindata = cell.to_gds(library.unit / library.precision, timestamp=now)
        assert bindata == bincells[name]
    with open(out, 'rb') as fin:
        bincells = gdspy.get_binary_cells(fin)
    for name, cell in library.cell_dict.items():
        bindata = cell.to_gds(library.unit / library.precision, timestamp=now)
        assert bindata == bincells[name]
