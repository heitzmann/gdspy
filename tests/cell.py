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
import numpy
import uuid


def unique():
    return str(uuid.uuid4())


@pytest.fixture
def tree():
    p1 = gdspy.Polygon(((0, 0), (0, 1), (1, 0)), 0, 0)
    p2 = gdspy.Polygon(((2, 0), (2, 1), (1, 0)), 1, 1)
    l1 = gdspy.Label('label1', (0, 0), layer=11)
    l2 = gdspy.Label('label2', (2, 1), layer=12)
    c1 = gdspy.Cell('tree_' + unique())
    c1.add(p1)
    c1.add(l1)
    c2 = gdspy.Cell('tree_' + unique())
    c2.add(l2)
    c2.add(p2)
    c2.add(gdspy.CellReference(c1))
    c3 = gdspy.Cell('tree_' + unique())
    c3.add(gdspy.CellArray(c2, 3, 2, (3, 3)))
    return c3, c2, c1


def test_duplicate():
    gdspy.current_library = gdspy.GdsLibrary()
    name = 'c_duplicate'
    c1 = gdspy.Cell(name)
    with pytest.raises(ValueError) as e:
        c2 = gdspy.Cell(name)
    assert name in str(e.value)


def test_ignore_duplicate():
    gdspy.current_library = gdspy.GdsLibrary()
    c1 = gdspy.Cell('c_ignore_duplicate')
    c2 = gdspy.Cell(c1.name, True)


def test_str():
    gdspy.current_library = gdspy.GdsLibrary()
    c = gdspy.Cell('c_str')
    assert str(c) == 'Cell ("c_str", 0 elements, 0 labels)'


def test_add_element():
    gdspy.current_library = gdspy.GdsLibrary()
    p = gdspy.Polygon(((0, 0), (1, 0), (0, 1)))
    c = gdspy.Cell('c_add_element')
    assert c.add(p) is c
    assert c.add([p, p]) is c
    assert c.elements == [p, p, p]


def test_add_label():
    gdspy.current_library = gdspy.GdsLibrary()
    l = gdspy.Label('label', (0, 0))
    c = gdspy.Cell('c_add_label')
    assert c.add(l) is c
    assert c.add([l, l]) is c
    assert c.labels == [l, l, l]


def test_copy():
    gdspy.current_library = gdspy.GdsLibrary()
    name = 'c_copy'
    p = gdspy.Polygon(((0, 0), (1, 0), (0, 1)))
    l = gdspy.Label('label', (0, 0))
    c1 = gdspy.Cell(name)
    c1.add(p)
    c1.add(l)
    with pytest.raises(ValueError) as e:
        c2 = c1.copy(name)
    assert name in str(e.value)
    c3 = c1.copy(name, True)
    assert c3.elements == c1.elements and c3.elements is not c1.elements
    assert c3.labels == c1.labels and c3.labels is not c1.labels
    cref = gdspy.Cell('c_ref').add(gdspy.Rectangle((-1, -1), (-2, -2)))
    c1.add(gdspy.CellReference(cref))
    bb = c1.get_bounding_box()
    c4 = c1.copy('c_copy_1', False, True)
    assert c4.elements != c1.elements
    assert c4.labels != c1.labels
    assert c1._bb_valid == True
    assert cref._bb_valid == True
    assert c4._bb_valid == False


def test_remove(tree):
    c3, c2, c1 = tree
    c1.remove_polygons(lambda p, l, d: l == 1)
    assert len(c1.elements) == 1
    c1.remove_polygons(lambda p, l, d: l == 0)
    assert len(c1.elements) == 0
    c1.remove_labels(lambda lbl: lbl.layer == 12)
    assert len(c1.labels) == 1
    c1.remove_labels(lambda lbl: lbl.layer == 11)
    assert len(c1.labels) == 0


def test_area():
    gdspy.current_library = gdspy.GdsLibrary()
    c = gdspy.Cell('c_area')
    c.add(gdspy.Rectangle((0, 0), (1, 1), layer=0))
    c.add(gdspy.Rectangle((0, 0), (1, 1), layer=1))
    c.add(gdspy.Rectangle((1, 1), (2, 2), layer=1))
    c.add(gdspy.Rectangle((1, 1), (2, 2), datatype=2))
    assert c.area() == 4.0
    assert c.area(True) == {(0, 0): 1.0, (1, 0): 2.0, (0, 2): 1}


def test_flatten_00(tree):
    c3, c2, c1 = tree
    c3.flatten()
    assert len(c3.elements) == 2
    for i in range(2):
        assert (c3.elements[i].layers == [0] * 6 or
                c3.elements[i].layers == [1] * 6)
        assert c3.elements[i].layers == c3.elements[i].datatypes
    assert len(c3.labels) == 12


def test_flatten_01(tree):
    c3, c2, c1 = tree
    c3.flatten(None, 2, 3)
    assert len(c3.elements) == 2
    for i in range(2):
        assert (c3.elements[i].layers == [0] * 6 or
                c3.elements[i].layers == [1] * 6)
        assert c3.elements[i].datatypes == [2] * 6
    assert len(c3.labels) == 12
    assert all(lbl.texttype == 3 for lbl in c3.labels)


def test_flatten_10(tree):
    c3, c2, c1 = tree
    c3.flatten(2)
    assert len(c3.elements) == 2
    for i in range(2):
        assert (c3.elements[i].datatypes == [0] * 6 or
                c3.elements[i].datatypes == [1] * 6)
        assert c3.elements[i].layers == [2] * 6
    assert len(c3.labels) == 12
    assert all(lbl.layer == 2 for lbl in c3.labels)


def test_flatten_11(tree):
    c3, c2, c1 = tree
    c3.flatten(2, 3, 4)
    assert len(c3.elements) == 1
    assert c3.elements[0].layers == [2] * 12
    assert c3.elements[0].datatypes == [3] * 12
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

def test_get_polygons2():
    gdspy.current_library = gdspy.GdsLibrary()
    c0 = gdspy.Cell('empty', False)
    assert len(c0.get_polygons()) == 0
    assert len(c0.get_polygons(True)) == 0
    assert len(c0.get_polygons(False, -1)) == 0


