######################################################################
#                                                                    #
#  Copyright 2009 Lucas Heitzmann Gabrielli.                         #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################
import datetime
import hashlib
import pathlib

import pytest
import gdspy
import numpy
import uuid
import os

gdspy.library.use_current_library = False

tmpdir = pathlib.Path(__file__).parent / 'tmp'
if not tmpdir.is_dir():
    tmpdir.mkdir()


def unique():
    return str(uuid.uuid4())


@pytest.fixture
def tree():
    p1 = gdspy.Polygon(((0, 0), (0, 1), (1, 0)), 0, 0)
    p2 = gdspy.Polygon(((2, 0), (2, 1), (1, 0)), 1, 1)
    l1 = gdspy.Label("label1", (0, 0), layer=11)
    l2 = gdspy.Label("label2", (2, 1), layer=12)
    c1 = gdspy.Cell("tree_" + unique())
    c1.add(p1)
    c1.add(l1)
    c2 = gdspy.Cell("tree_" + unique())
    c2.add(l2)
    c2.add(p2)
    c2.add(gdspy.CellReference(c1))
    c3 = gdspy.Cell("tree_" + unique())
    c3.add(gdspy.CellArray(c2, 3, 2, (3, 3)))
    return c3, c2, c1


def test_duplicate():
    lib = gdspy.GdsLibrary()
    name = "c_duplicate"
    c = gdspy.Cell(name)
    lib.add(c)
    with pytest.raises(ValueError) as e:
        lib.add(gdspy.Cell(name), overwrite_duplicate=False)
    assert name in str(e.value)


def test_ignore_duplicate():
    lib = gdspy.GdsLibrary()
    c1 = gdspy.Cell("c_ignore_duplicate")
    lib.add(c1)
    lib.add(c1, overwrite_duplicate=False)
    lib.add(gdspy.Cell(c1.name), overwrite_duplicate=True)
    assert lib.cells[c1.name] is not c1


def test_str():
    c = gdspy.Cell("c_str")
    assert str(c) == 'Cell ("c_str", 0 polygons, 0 paths, 0 labels, 0 references)'


def test_add_element():
    p = gdspy.Polygon(((0, 0), (1, 0), (0, 1)))
    c = gdspy.Cell("c_add_element")
    assert c.add(p) is c
    assert c.add([p, p]) is c
    assert c.polygons == [p, p, p]


def test_add_label():
    lbl = gdspy.Label("label", (0, 0))
    c = gdspy.Cell("c_add_label")
    assert c.add(lbl) is c
    assert c.add([lbl, lbl]) is c
    assert c.labels == [lbl, lbl, lbl]


def test_copy():
    name = "c_copy"
    p = gdspy.Polygon(((0, 0), (1, 0), (0, 1)))
    lbl = gdspy.Label("label", (0, 0))
    c1 = gdspy.Cell(name)
    c1.add(p)
    c1.add(lbl)
    c3 = c1.copy(name, False)
    assert c3.polygons == c1.polygons and c3.polygons is not c1.polygons
    assert c3.labels == c1.labels and c3.labels is not c1.labels
    cref = gdspy.Cell("c_ref").add(gdspy.Rectangle((-1, -1), (-2, -2)))
    c1.add(gdspy.CellReference(cref))
    c1.get_bounding_box()
    c4 = c1.copy("c_copy_1", True)
    assert c4.polygons != c1.polygons
    assert c4.labels != c1.labels
    assert c1._bb_valid
    assert cref._bb_valid
    assert not c4._bb_valid


def test_remove(tree):
    c3, c2, c1 = tree
    c1.remove_polygons(lambda p, layer, d: layer == 1)
    assert len(c1.polygons) == 1
    c1.remove_polygons(lambda p, layer, d: layer == 0)
    assert len(c1.polygons) == 0
    c1.remove_labels(lambda lbl: lbl.layer == 12)
    assert len(c1.labels) == 1
    c1.remove_labels(lambda lbl: lbl.layer == 11)
    assert len(c1.labels) == 0


def test_area():
    c = gdspy.Cell("c_area")
    c.add(gdspy.Rectangle((0, 0), (1, 1), layer=0))
    c.add(gdspy.Rectangle((0, 0), (1, 1), layer=1))
    c.add(gdspy.Rectangle((1, 1), (2, 2), layer=1))
    c.add(gdspy.Rectangle((1, 1), (2, 2), datatype=2))
    assert c.area() == 4.0
    assert c.area(True) == {(0, 0): 1.0, (1, 0): 2.0, (0, 2): 1}


def test_flatten_00(tree):
    c3, c2, c1 = tree
    c3.flatten()
    assert len(c3.polygons) == 12
    for i in range(12):
        assert c3.polygons[i].layers == [0] or c3.polygons[i].layers == [1]
        assert c3.polygons[i].layers == c3.polygons[i].datatypes
    assert len(c3.labels) == 12


def test_flatten_01(tree):
    c3, c2, c1 = tree
    c3.flatten(None, 2, 3)
    assert len(c3.polygons) == 12
    for i in range(12):
        assert c3.polygons[i].layers == [0] or c3.polygons[i].layers == [1]
        assert c3.polygons[i].datatypes == [2]
    assert len(c3.labels) == 12
    assert all(lbl.texttype == 3 for lbl in c3.labels)


def test_flatten_10(tree):
    c3, c2, c1 = tree
    c3.flatten(2)
    assert len(c3.polygons) == 12
    for i in range(12):
        assert c3.polygons[i].datatypes == [0] or c3.polygons[i].datatypes == [1]
        assert c3.polygons[i].layers == [2]
    assert len(c3.labels) == 12
    assert all(lbl.layer == 2 for lbl in c3.labels)


def test_flatten_11(tree):
    c3, c2, c1 = tree
    c3.flatten(2, 3, 4)
    assert len(c3.polygons) == 12
    assert all(p.layers == [2] for p in c3.polygons)
    assert all(p.datatypes == [3] for p in c3.polygons)
    assert len(c3.labels) == 12
    assert all(lbl.layer == 2 for lbl in c3.labels)
    assert all(lbl.texttype == 4 for lbl in c3.labels)


def test_bb(tree):
    c3, c2, c1 = tree
    err = numpy.array(((0, 0), (8, 4))) - c3.get_bounding_box()
    assert numpy.max(numpy.abs(err)) == 0

    p2 = gdspy.Polygon(((-1, 2), (-1, 1), (0, 2)), 2, 2)
    c2.add(p2)
    err = numpy.array(((-1, 0), (8, 5))) - c3.get_bounding_box()
    assert numpy.max(numpy.abs(err)) == 0

    p1 = gdspy.Polygon(((0, 3), (0, 2), (1, 3)), 3, 3)
    c1.add(p1)
    err = numpy.array(((-1, 0), (8, 6))) - c3.get_bounding_box()
    assert numpy.max(numpy.abs(err)) == 0


def test_layers(tree):
    assert tree[0].get_layers() == {0, 1, 11, 12}


def test_datatypes(tree):
    assert tree[0].get_datatypes() == {0, 1}


def test_get_polygons1(tree):
    c3, c2, c1 = tree
    p1 = gdspy.Polygon(((0, 3), (0, 2), (1, 3)), 3, 3)
    c1.add(p1)
    assert len(c3.get_polygons()) == 18
    assert len(c3.get_polygons(False, 0)) == 6
    assert len(c3.get_polygons(False, 1)) == 12
    assert set(c3.get_polygons(True).keys()) == {(0, 0), (1, 1), (3, 3)}
    assert set(c3.get_polygons(True, 0).keys()) == {c2.name}
    assert set(c3.get_polygons(True, 1).keys()) == {c1.name, (1, 1)}


def test_get_polygons2(tree):
    c3, c2, c1 = tree
    c1.add(gdspy.Rectangle((0, 0), (1, 1), 0, 0))
    assert len(c1.get_polygons()) == 2
    d = c1.get_polygons(True)
    assert len(d) == 1
    assert (0, 0) in d
    assert len(d[(0, 0)]) == 2
    c3.add(gdspy.CellReference(c1))
    d = c3.get_polygons(True)
    assert len(d) == 2
    assert (0, 0) in d and (1, 1) in d
    assert len(d[(0, 0)]) == 14
    assert len(d[(1, 1)]) == 6


def test_get_polygons3():
    c0 = gdspy.Cell("empty")
    assert len(c0.get_polygons()) == 0
    assert len(c0.get_polygons(True)) == 0
    assert len(c0.get_polygons(False, -1)) == 0


def test_get_polygons4(tree):
    c3, c2, c1 = tree
    c3.add(gdspy.Rectangle((0, 0), (1, 1), 0, 0))
    assert len(c3.get_polygons((0, 0))) == 7
    assert len(c3.get_polygons((0, 0), 0)) == 1
    assert len(c3.get_polygons((1, 1), 0)) == 0
    assert len(c3.get_polygons((0, 0), 1)) == 1
    assert len(c3.get_polygons((1, 1), 1)) == 6


def test_write_svg(tree, tmpdir):
    _, _, c1 = tree
    fname = str(tmpdir.join("c1.svg"))
    c1.write_svg(fname)
    assert os.path.isfile(fname)
    assert not os.stat(fname).st_size == 0


def test_write_svg_with_style(tree, tmpdir):
    _, _, c1 = tree
    fname = str(tmpdir.join("c1.svg"))
    style = {(0, 0): {"fill": "CC00FF"}}
    c1.write_svg(fname, style=style)
    assert os.path.isfile(fname)
    assert not os.stat(fname).st_size == 0


def assert_bb(bb1, bb2):
    for p1, p2 in zip(bb1, bb2):
        for c1, c2 in zip(p1, p2):
            assert abs(c1 - c2) < 1e-12


def test_bounding_box():
    cell1 = gdspy.Cell("1")
    cell1.add(gdspy.Rectangle((0, 0), (1, 1)))

    # Cell1: plain 1x1 um square, with bottom-left corner at 0,0
    assert_bb(cell1.get_bounding_box(), ((0, 0), (1, 1)))

    # Cell2: place and rotate Cell1 45deg about origin
    cell2 = gdspy.Cell("2")
    cell2.add(gdspy.CellReference(cell1, rotation=45))

    assert_bb(cell2.get_bounding_box(), ((-(0.5 ** 0.5), 0), (0.5 ** 0.5, 2 ** 0.5)))

    # Cell3: place and rotate Cell2 an additional 45deg about origin (Cell1 is now rotated total 90deg)
    cell3 = gdspy.Cell("3")
    cell3.add(gdspy.CellReference(cell2, rotation=45))

    assert_bb(cell3.get_bounding_box(), ((-1, 0), (0, 1)))

    # Cell4: nest Cell2 one level deeper with no transform (geometric equivalent to Cell2)
    cell4 = gdspy.Cell("4")
    cell4.add(gdspy.CellReference(cell2))

    assert_bb(cell4.get_bounding_box(), ((-(0.5 ** 0.5), 0), (0.5 ** 0.5, 2 ** 0.5)))

    # Cell5: rotate Cell4 an addition 45 degrees (geometric equivalent to Cell3)
    cell5 = gdspy.Cell("5")
    cell5.add(gdspy.CellReference(cell4, rotation=45))

    assert_bb(cell5.get_bounding_box(), ((-1, 0), (0, 1)))

    # Cell6: translate Cell1 by 2um east
    cell6 = gdspy.Cell("6")
    cell6.add(gdspy.CellReference(cell1, origin=(2, 0)))

    assert_bb(cell6.get_bounding_box(), ((2, 0), (3, 1)))

    cell7a = gdspy.Cell("7a")

    assert cell7a.get_bounding_box() is None

    cell7b = gdspy.Cell("7b")
    cell7b.add(gdspy.CellReference(cell7a, rotation=90))

    assert cell7b.get_bounding_box() is None

    cell7c = gdspy.Cell("7c")
    cell7c.add(gdspy.CellReference(cell7b))

    assert cell7c.get_bounding_box() is None

    cell7 = gdspy.Cell("7")
    cell7.add(cell1)
    cell7.add(cell7c)

    assert_bb(cell7.get_bounding_box(), ((0, 0), (1, 1)))


def test_bounding_box2():
    cell0 = gdspy.Cell("0")
    cell0.add(gdspy.Rectangle((0, 0), (0.5, 0.25)))

    cell1 = gdspy.Cell("1")
    cell1.add(gdspy.CellArray(cell0, 2, 4, (0.5, 0.25)))

    # Cell1: plain 1x1 um square, with bottom-left corner at 0,0
    assert_bb(cell1.get_bounding_box(), ((0, 0), (1, 1)))

    # Cell2: place and rotate Cell1 45deg about origin
    cell2 = gdspy.Cell("2")
    cell2.add(gdspy.CellReference(cell1, rotation=45))

    assert_bb(cell2.get_bounding_box(), ((-(0.5 ** 0.5), 0), (0.5 ** 0.5, 2 ** 0.5)))

    # Cell3: place and rotate Cell2 an additional 45deg about origin (Cell1 is now rotated total 90deg)
    cell3 = gdspy.Cell("3")
    cell3.add(gdspy.CellReference(cell2, rotation=45))

    assert_bb(cell3.get_bounding_box(), ((-1, 0), (0, 1)))

    # Cell4: nest Cell2 one level deeper with no transform (geometric equivalent to Cell2)
    cell4 = gdspy.Cell("4")
    cell4.add(gdspy.CellArray(cell2, 1, 1, (1, 1)))

    assert_bb(cell4.get_bounding_box(), ((-(0.5 ** 0.5), 0), (0.5 ** 0.5, 2 ** 0.5)))

    # Cell5: rotate Cell4 an addition 45 degrees (geometric equivalent to Cell3)
    cell5 = gdspy.Cell("5")
    cell5.add(gdspy.CellArray(cell4, 1, 1, (1, 1), rotation=45))

    assert_bb(cell5.get_bounding_box(), ((-1, 0), (0, 1)))

    # Cell6: translate Cell1 by 2um east
    cell6 = gdspy.Cell("6")
    cell6.add(gdspy.CellArray(cell1, 1, 1, (1, 1), origin=(2, 0)))

    assert_bb(cell6.get_bounding_box(), ((2, 0), (3, 1)))


def hash_file(filepath):
    filepath = pathlib.Path(filepath)
    md5 = hashlib.md5()
    md5.update(filepath.read_bytes())
    return md5.hexdigest()


def test_time_changes_gds_hash():
    fn1 = tmpdir / 'nofreeze1.gds'
    fn2 = tmpdir / 'nofreeze2.gds'
    date1 = datetime.datetime(1988, 8, 28)
    date2 = datetime.datetime(2020, 12, 25)
    lib = gdspy.GdsLibrary(name='speedy')
    lib.write_gds(fn1, timestamp=date1)
    hash1 = hash_file(fn1)
    lib.write_gds(fn2, timestamp=date2)
    hash2 = hash_file(fn2)

    assert hash1 != hash2


def test_frozen_gds_has_constant_hash():
    fn1 = tmpdir / 'freeze1.gds'
    fn2 = tmpdir / 'freeze2.gds'
    frozen_date = datetime.datetime(1988, 8, 28)
    lib = gdspy.GdsLibrary(name='Elsa')
    lib.write_gds(fn1, timestamp=frozen_date)
    hash1 = hash_file(fn1)
    lib.write_gds(fn2, timestamp=frozen_date)
    hash2 = hash_file(fn2)

    assert hash1 == hash2


def test_frozen_gds_with_cell_has_constant_hash():
    fn1 = tmpdir / 'freezec1.gds'
    fn2 = tmpdir / 'freezec2.gds'
    frozen_date = datetime.datetime(1988, 8, 28)
    lib = gdspy.GdsLibrary(name='Elsa')
    cell = gdspy.Cell(name='Anna')
    cell.add(gdspy.Rectangle((0, 0), (100, 1000)))
    lib.add(cell)
    lib.write_gds(fn1, timestamp=frozen_date)
    hash1 = hash_file(fn1)
    lib.write_gds(fn2, timestamp=frozen_date)
    hash2 = hash_file(fn2)

    assert hash1 == hash2


def test_frozen_gds_with_cell_array_has_constant_hash():
    fn1 = tmpdir / 'freezea1.gds'
    fn2 = tmpdir / 'freezea2.gds'
    frozen_date = datetime.datetime(1988, 8, 28)
    lib = gdspy.GdsLibrary(name='Elsa')
    cell = gdspy.Cell(name='Anna')
    cell.add(gdspy.Rectangle((0, 0), (100, 1000)))
    cell2 = gdspy.Cell(name='Olaf')
    cell2.add(gdspy.Rectangle((0, 0), (50, 100)))

    cell_array = gdspy.CellArray(ref_cell=cell2, columns=5, rows=2, spacing=(60, 120), origin=(1000, 0))
    cell.add(cell_array)
    lib.add(cell)
    lib.write_gds(fn1, timestamp=frozen_date)
    hash1 = hash_file(fn1)
    lib.write_gds(fn2, timestamp=frozen_date)
    hash2 = hash_file(fn2)

    assert hash1 == hash2
