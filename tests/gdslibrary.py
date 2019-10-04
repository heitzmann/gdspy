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

gdspy.library.use_current_library = False


def unique():
    return str(uuid.uuid4())


def test_add():
    lib = gdspy.GdsLibrary()
    c1 = gdspy.Cell("gl_add_1")
    c2 = gdspy.Cell("gl_add_2")
    c3 = gdspy.Cell("gl_add_3")
    r1 = gdspy.CellReference(c1)
    with pytest.warns(UserWarning):
        r2 = gdspy.CellReference("gl_add_2")
    c3.add([r1, r2])
    lib.add(c1)
    lib.add((c2, c3))
    assert lib.cells == {"gl_add_1": c1, "gl_add_2": c2, "gl_add_3": c3}
    lib = gdspy.GdsLibrary()
    lib.add(c3)
    assert lib.cells == {"gl_add_1": c1, "gl_add_3": c3}


def test_duplicate():
    name = "gl_duplicate"
    lib = gdspy.GdsLibrary()
    lib.add(gdspy.Cell(name))
    c = gdspy.Cell(name)
    with pytest.raises(ValueError) as e:
        lib.add(c)
    assert name in str(e.value)

    lib.add(c, overwrite_duplicate=True)
    assert lib.cells == {name: c}

    cl = [gdspy.Cell(name), gdspy.Cell(name + "1")]
    with pytest.raises(ValueError) as e:
        lib.add(cl)
    assert name in str(e.value)

    lib.add(cl, overwrite_duplicate=True)
    assert lib.cells == {name: cl[0], name + "1": cl[1]}


def test_add_update():
    lib = gdspy.GdsLibrary()
    main = gdspy.Cell("MAIN")
    c1 = gdspy.Cell("C1")
    c2 = gdspy.Cell("C2")
    c3 = gdspy.Cell("C1")
    r1 = gdspy.CellReference(c1)
    main.add(r1)
    with pytest.warns(UserWarning):
        r2 = gdspy.CellArray("C1", 1, 1, (1, 1))
    main.add(r2)
    r3 = gdspy.CellReference(c2)
    main.add(r3)
    r4 = gdspy.CellReference(c3)
    c2.add(r4)
    with pytest.warns(UserWarning):
        r5 = gdspy.CellReference("C3")
    c1.add(r5)
    with pytest.warns(UserWarning):
        r6 = gdspy.CellReference("C2")
    main.add(r6)
    lib.add([main, c1, c2], include_dependencies=False)
    lib.add(c3, include_dependencies=False, overwrite_duplicate=True)
    assert r1.ref_cell is c3
    assert r2.ref_cell is c3
    assert r3.ref_cell is c2
    assert r4.ref_cell is c3
    assert r5.ref_cell == "C3"
    assert r6.ref_cell == "C2"


def test_add_update2():
    lib = gdspy.GdsLibrary()
    main = gdspy.Cell("MAIN")
    c1 = gdspy.Cell("C1")
    c2 = gdspy.Cell("C2")
    c3 = gdspy.Cell("C1")
    r1 = gdspy.CellReference(c1)
    main.add(r1)
    with pytest.warns(UserWarning):
        r2 = gdspy.CellArray("C1", 1, 1, (1, 1))
    main.add(r2)
    r3 = gdspy.CellReference(c2)
    main.add(r3)
    r4 = gdspy.CellReference(c3)
    c2.add(r4)
    with pytest.warns(UserWarning):
        r5 = gdspy.CellReference("C3")
    c1.add(r5)
    with pytest.warns(UserWarning):
        r6 = gdspy.CellReference("C2")
    main.add(r6)
    lib.add([main, c1, c2], include_dependencies=False)
    lib.add(
        c3,
        include_dependencies=False,
        overwrite_duplicate=True,
        update_references=False,
    )
    assert r1.ref_cell is c1
    assert r2.ref_cell == "C1"
    assert r3.ref_cell is c2
    assert r4.ref_cell is c3
    assert r5.ref_cell == "C3"
    assert r6.ref_cell == "C2"


def test_remove():
    lib = gdspy.GdsLibrary()
    main = gdspy.Cell("MAIN")
    c1 = gdspy.Cell("C1")
    c2 = gdspy.Cell("C2")
    c3 = gdspy.Cell("C1")
    r1 = gdspy.CellReference(c1)
    main.add(r1)
    with pytest.warns(UserWarning):
        r2 = gdspy.CellArray("C1", 1, 1, (1, 1))
    main.add(r2)
    r3 = gdspy.CellReference(c2)
    main.add(r3)
    r4 = gdspy.CellReference(c3)
    c2.add(r4)
    with pytest.warns(UserWarning):
        r5 = gdspy.CellReference("C3")
    c1.add(r5)
    with pytest.warns(UserWarning):
        r6 = gdspy.CellReference("C2")
    main.add(r6)
    lib.add([main, c1, c2], include_dependencies=False)
    assert lib.remove("C3") == 1
    assert len(c1.references) == 0
    assert len(c2.references) == 1
    assert c2.references[0] is r4
    assert lib.remove(c1) == 3
    assert "C1" not in lib.cells
    assert len(main.references) == 2
    assert main.references[0] is r3
    assert main.references[1] is r6
    assert len(c2.references) == 0


def test_replace():
    lib = gdspy.GdsLibrary()
    main = gdspy.Cell("MAIN")
    c1 = gdspy.Cell("C1")
    c2 = gdspy.Cell("C2")
    c3 = gdspy.Cell("C1")
    r1 = gdspy.CellReference(c1)
    main.add(r1)
    with pytest.warns(UserWarning):
        r2 = gdspy.CellArray("C1", 1, 1, (1, 1))
    main.add(r2)
    r3 = gdspy.CellReference(c2)
    main.add(r3)
    r4 = gdspy.CellReference(c3)
    c2.add(r4)
    with pytest.warns(UserWarning):
        r5 = gdspy.CellReference("C3")
    c1.add(r5)
    with pytest.warns(UserWarning):
        r6 = gdspy.CellReference("C2")
    main.add(r6)
    lib.add([main, c1, c2], include_dependencies=False)
    assert lib.replace_references(c2, c3) == 2
    assert r3.ref_cell is c3
    assert r6.ref_cell is c3
    assert lib.replace_references("C3", "C1") == 1
    assert r5.ref_cell is c1
    assert lib.replace_references("C1", c2) == 6
    assert r1.ref_cell is c2
    assert r2.ref_cell is c2
    assert r3.ref_cell is c2
    assert r4.ref_cell is c2
    assert r5.ref_cell is c2
    assert r6.ref_cell is c2


def test_rename():
    lib = gdspy.GdsLibrary()
    main = gdspy.Cell("MAIN")
    c1 = gdspy.Cell("C1")
    c2 = gdspy.Cell("C2")
    c3 = gdspy.Cell("C1")
    r1 = gdspy.CellReference(c1)
    main.add(r1)
    with pytest.warns(UserWarning):
        r2 = gdspy.CellArray("C1", 1, 1, (1, 1))
    main.add(r2)
    r3 = gdspy.CellReference(c2)
    main.add(r3)
    r4 = gdspy.CellReference(c3)
    c2.add(r4)
    with pytest.warns(UserWarning):
        r5 = gdspy.CellReference("C3")
    c1.add(r5)
    with pytest.warns(UserWarning):
        r6 = gdspy.CellReference("C2")
    main.add(r6)
    lib.add([main, c1, c2], include_dependencies=False)
    with pytest.raises(ValueError):
        lib.rename_cell(c3, "C3")
    assert c3.name == "C1"
    with pytest.raises(ValueError):
        lib.rename_cell(c2, "C1")
    assert c2.name == "C2"
    assert lib.rename_cell(c1, "C3") == 2
    assert c1.name == "C3"
    assert lib.cells["C3"] is c1
    assert lib.rename_cell(c2, "X2", False) == 0
    assert c2.name == "X2"
    assert lib.cells["X2"] is c2
    assert r1.ref_cell is c1
    assert r2.ref_cell is c1
    assert r3.ref_cell is c2
    assert r4.ref_cell is c1
    assert r5.ref_cell == "C3"
    assert r6.ref_cell == "C2"


@pytest.fixture
def tree():
    c = [gdspy.Cell("tree_" + unique()) for _ in range(8)]
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


def test_rw_gds(tmpdir):
    lib = gdspy.GdsLibrary("lib", unit=2e-3, precision=1e-5)
    c1 = gdspy.Cell("gl_rw_gds_1")
    c1.add(gdspy.Rectangle((0, -1), (1, 2), 2, 4))
    c1.add(gdspy.Label("label", (1, -1), "w", 45, 1.5, True, 5, 6))
    c2 = gdspy.Cell("gl_rw_gds_2")
    c2.add(gdspy.Round((0, 0), 1, number_of_points=32, max_points=20))
    c3 = gdspy.Cell("gl_rw_gds_3")
    c3.add(gdspy.CellReference(c1, (0, 1), -90, 2, True))
    c4 = gdspy.Cell("gl_rw_gds_4")
    c4.add(gdspy.CellArray(c2, 2, 3, (1, 4), (-1, -2), 180, 0.5, True))
    lib.add((c1, c2, c3, c4))

    fname1 = str(tmpdir.join("test1.gds"))
    lib.write_gds(fname1)
    lib1 = gdspy.GdsLibrary(
        infile=fname1,
        unit=1e-3,
        precision=1e-6,
        units="convert",
        rename={"gl_rw_gds_1": "1"},
        layers={2: 4},
        datatypes={4: 2},
        texttypes={6: 7},
    )
    assert lib1.name == "lib"
    assert len(lib1.cells) == 4
    assert set(lib1.cells.keys()) == {"1", "gl_rw_gds_2", "gl_rw_gds_3", "gl_rw_gds_4"}
    c = lib1.cells["1"]
    assert len(c.polygons) == len(c.labels) == 1
    assert c.polygons[0].area() == 12.0
    assert c.polygons[0].layers == [4]
    assert c.polygons[0].datatypes == [2]
    assert c.labels[0].text == "label"
    assert c.labels[0].position[0] == 2 and c.labels[0].position[1] == -2
    assert c.labels[0].anchor == 4
    assert c.labels[0].rotation == 45
    assert c.labels[0].magnification == 1.5
    assert c.labels[0].x_reflection == True
    assert c.labels[0].layer == 5
    assert c.labels[0].texttype == 7

    c = lib1.cells["gl_rw_gds_2"]
    assert len(c.polygons) == 2
    assert isinstance(c.polygons[0], gdspy.Polygon) and isinstance(
        c.polygons[1], gdspy.Polygon
    )

    c = lib1.cells["gl_rw_gds_3"]
    assert len(c.references) == 1
    assert isinstance(c.references[0], gdspy.CellReference)
    assert c.references[0].ref_cell == lib1.cells["1"]
    assert c.references[0].origin[0] == 0 and c.references[0].origin[1] == 2
    assert c.references[0].rotation == -90
    assert c.references[0].magnification == 2
    assert c.references[0].x_reflection == True

    c = lib1.cells["gl_rw_gds_4"]
    assert len(c.references) == 1
    assert isinstance(c.references[0], gdspy.CellArray)
    assert c.references[0].ref_cell == lib1.cells["gl_rw_gds_2"]
    assert c.references[0].origin[0] == -2 and c.references[0].origin[1] == -4
    assert c.references[0].rotation == 180
    assert c.references[0].magnification == 0.5
    assert c.references[0].x_reflection == True
    assert c.references[0].spacing[0] == 2 and c.references[0].spacing[1] == 8
    assert c.references[0].columns == 2
    assert c.references[0].rows == 3

    fname2 = str(tmpdir.join("test2.gds"))
    lib.name = "lib2"
    with open(fname2, "wb") as fout:
        lib.write_gds(fout)
    with open(fname2, "rb") as fin:
        lib2 = gdspy.GdsLibrary()
        lib2.read_gds(fin)
    assert lib2.name == "lib2"
    assert len(lib2.cells) == 4
