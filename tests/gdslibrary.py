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
import uuid


def unique():
    return str(uuid.uuid4())


def test_add():
    lib = gdspy.GdsLibrary()
    c1 = gdspy.Cell('gl_add_1', exclude_from_current=True)
    c2 = gdspy.Cell('gl_add_2', exclude_from_current=True)
    c3 = gdspy.Cell('gl_add_3', exclude_from_current=True)
    lib.add(c1)
    lib.add((c2, c3))
    assert lib.cell_dict == {'gl_add_1': c1, 'gl_add_2': c2, 'gl_add_3': c3}


def test_duplicate():
    name = 'gl_duplicate'
    lib = gdspy.GdsLibrary()
    lib.add(gdspy.Cell(name, exclude_from_current=True))
    c = gdspy.Cell(name, exclude_from_current=True)
    with pytest.raises(ValueError) as e:
        lib.add(c)
    assert name in str(e.value)

    lib.add(c, True)
    assert lib.cell_dict == {name: c}

    cl = [gdspy.Cell(name, exclude_from_current=True), gdspy.Cell(name + '1', exclude_from_current=True)]
    with pytest.raises(ValueError) as e:
        lib.add(cl)
    assert name in str(e.value)

    lib.add(cl, True)
    assert lib.cell_dict == {name: cl[0], name + '1': cl[1]}


@pytest.fixture
def tree():
    c = [gdspy.Cell('tree_' + unique(), True) for _ in range(8)]
    lib = gdspy.GdsLibrary()
    lib.add(c)
    c[0].add(gdspy.CellReference(c[1]))
    c[0].add(gdspy.CellReference(c[3]))
    c[1].add(gdspy.CellReference(c[2]))
    c[1].add(gdspy.CellArray(c[2], 2, 1, (0, 0)))
    c[1].add(gdspy.CellArray(c[3], 2, 1, (0, 0)))
    c[4].add(gdspy.CellReference(c[3]))
    c[6].add(gdspy.CellArray(c[5], 2, 1, (0, 0)))
    return lib, c


def test_top_level_1(tree):
    lib, c = tree
    tl = lib.top_level()
    assert len(tl) == 4 and c[0] in tl and c[4] in tl and c[6] in tl and c[7] in tl


def test_top_level_2(tree):
    lib, c = tree
    c[7].add(gdspy.CellReference(c[0]))
    c[7].add(gdspy.CellReference(c[4]))
    c[7].add(gdspy.CellReference(c[6]))
    assert lib.top_level() == [c[7]]


def test_top_level_3(tree):
    lib, c = tree
    c[7].add(gdspy.CellReference(c[0]))
    c[3].add(gdspy.CellReference(c[4]))
    c[2].add(gdspy.CellReference(c[6]))
    c[1].add(gdspy.CellReference(c[7]))
    assert lib.top_level() == []


def test_extract():
    c = [gdspy.Cell('tree_' + unique(), True) for _ in range(8)]
    gdspy.current_library = gdspy.GdsLibrary()
    lib = gdspy.GdsLibrary()
    lib.add(c)
    c[0].add(gdspy.CellReference(c[1]))
    c[0].add(gdspy.CellReference(c[3]))
    c[1].add(gdspy.CellReference(c[2]))
    c[1].add(gdspy.CellArray(c[2], 2, 1, (0, 0)))
    c[1].add(gdspy.CellArray(c[3], 2, 1, (0, 0)))
    c[4].add(gdspy.CellReference(c[3]))
    c[6].add(gdspy.CellArray(c[5], 2, 1, (0, 0)))
    assert len(gdspy.current_library.cell_dict) == 0

    lib.extract(c[7])
    assert gdspy.current_library.cell_dict == {c[7].name: c[7]}

    lib.extract(c[1])
    assert gdspy.current_library.cell_dict == {c[7].name: c[7], c[1].name: c[1], c[2].name: c[2], c[3].name: c[3]}

    lib.extract(c[0])
    assert gdspy.current_library.cell_dict == {
        c[7].name: c[7],
        c[0].name: c[0],
        c[1].name: c[1],
        c[2].name: c[2],
        c[3].name: c[3]
    }


def test_rw_gds(tmpdir):
    lib = gdspy.GdsLibrary('lib', unit=2e-3, precision=1e-5)
    c1 = gdspy.Cell('gl_rw_gds_1', True)
    c1.add(gdspy.Rectangle((0, -1), (1, 2), 2, 4))
    c1.add(gdspy.Label('label', (1, -1), 'w', 45, 1.5, True, 5, 6))
    c2 = gdspy.Cell('gl_rw_gds_2', True)
    c2.add(gdspy.Round((0, 0), 1, number_of_points=32, max_points=20))
    c3 = gdspy.Cell('gl_rw_gds_3', True)
    c3.add(gdspy.CellReference(c1, (0, 1), -90, 2, True))
    c4 = gdspy.Cell('gl_rw_gds_4', True)
    c4.add(gdspy.CellArray(c2, 2, 3, (1, 4), (-1, -2), 180, 0.5, True))
    lib.add((c1, c2, c3, c4))

    fname1 = str(tmpdir.join('test1.gds'))
    lib.write_gds(fname1)
    lib1 = gdspy.GdsLibrary(infile=fname1, unit=1e-3, precision=1e-6, units='convert',
                            rename={'gl_rw_gds_1': '1'}, layers={2: 4}, datatypes={4: 2},
                            texttypes={6: 7})
    assert lib1.name == 'lib'
    assert len(lib1.cell_dict) == 4
    assert set(lib1.cell_dict.keys()) == {'1', 'gl_rw_gds_2', 'gl_rw_gds_3', 'gl_rw_gds_4'}
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

    c = lib1.cell_dict['gl_rw_gds_2']
    assert len(c.polygons) == 2
    assert isinstance(c.polygons[0], gdspy.Polygon) and isinstance(c.polygons[1], gdspy.Polygon)

    c = lib1.cell_dict['gl_rw_gds_3']
    assert len(c.references) == 1
    assert isinstance(c.references[0], gdspy.CellReference)
    assert c.references[0].ref_cell == lib1.cell_dict['1']
    assert c.references[0].origin[0] == 0 and c.references[0].origin[1] == 2
    assert c.references[0].rotation == -90
    assert c.references[0].magnification == 2
    assert c.references[0].x_reflection == True

    c = lib1.cell_dict['gl_rw_gds_4']
    assert len(c.references) == 1
    assert isinstance(c.references[0], gdspy.CellArray)
    assert c.references[0].ref_cell == lib1.cell_dict['gl_rw_gds_2']
    assert c.references[0].origin[0] == -2 and c.references[0].origin[1] == -4
    assert c.references[0].rotation == 180
    assert c.references[0].magnification == 0.5
    assert c.references[0].x_reflection == True
    assert c.references[0].spacing[0] == 2 and c.references[0].spacing[1] == 8
    assert c.references[0].columns == 2
    assert c.references[0].rows == 3

    fname2 = str(tmpdir.join('test2.gds'))
    lib.name = 'lib2'
    with open(fname2, 'wb') as fout:
        lib.write_gds(fout)
    with open(fname2, 'rb') as fin:
        lib2 = gdspy.GdsLibrary()
        lib2.read_gds(fin)
    assert lib2.name == 'lib2'
    assert len(lib2.cell_dict) == 4
