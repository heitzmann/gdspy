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
import gdspy

gdspy.library.use_current_library = False


def test_writer_gds(tmpdir):
    lib = gdspy.GdsLibrary()
    c1 = gdspy.Cell("gw_rw_gds_1")
    c1.add(gdspy.Rectangle((0, -1), (1, 2), 2, 4))
    c1.add(gdspy.Label("label", (1, -1), "w", 45, 1.5, True, 5, 6))
    c2 = gdspy.Cell("gw_rw_gds_2")
    c2.add(gdspy.Round((0, 0), 1, number_of_points=32, max_points=20))
    c3 = gdspy.Cell("gw_rw_gds_3")
    c3.add(gdspy.CellReference(c1, (0, 1), -90, 2, True))
    c4 = gdspy.Cell("gw_rw_gds_4")
    c4.add(gdspy.CellArray(c2, 2, 3, (1, 4), (-1, -2), 180, 0.5, True))
    lib.add((c1, c2, c3, c4))

    fname1 = str(tmpdir.join("test1.gds"))
    writer1 = gdspy.GdsWriter(fname1, name="lib", unit=2e-3, precision=1e-5)
    for c in lib.cells.values():
        writer1.write_cell(c)
    writer1.close()
    lib1 = gdspy.GdsLibrary(unit=1e-3)
    lib1.read_gds(
        fname1,
        units="convert",
        rename={"gw_rw_gds_1": "1"},
        layers={2: 4},
        datatypes={4: 2},
        texttypes={6: 7},
    )
    assert lib1.name == "lib"
    assert len(lib1.cells) == 4
    assert set(lib1.cells.keys()) == {"1", "gw_rw_gds_2", "gw_rw_gds_3", "gw_rw_gds_4"}
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

    c = lib1.cells["gw_rw_gds_2"]
    assert len(c.polygons) == 2
    assert isinstance(c.polygons[0], gdspy.Polygon) and isinstance(
        c.polygons[1], gdspy.Polygon
    )

    c = lib1.cells["gw_rw_gds_3"]
    assert len(c.references) == 1
    assert isinstance(c.references[0], gdspy.CellReference)
    assert c.references[0].ref_cell == lib1.cells["1"]
    assert c.references[0].origin[0] == 0 and c.references[0].origin[1] == 2
    assert c.references[0].rotation == -90
    assert c.references[0].magnification == 2
    assert c.references[0].x_reflection == True

    c = lib1.cells["gw_rw_gds_4"]
    assert len(c.references) == 1
    assert isinstance(c.references[0], gdspy.CellArray)
    assert c.references[0].ref_cell == lib1.cells["gw_rw_gds_2"]
    assert c.references[0].origin[0] == -2 and c.references[0].origin[1] == -4
    assert c.references[0].rotation == 180
    assert c.references[0].magnification == 0.5
    assert c.references[0].x_reflection == True
    assert c.references[0].spacing[0] == 2 and c.references[0].spacing[1] == 8
    assert c.references[0].columns == 2
    assert c.references[0].rows == 3

    fname2 = str(tmpdir.join("test2.gds"))
    with open(fname2, "wb") as fout:
        writer2 = gdspy.GdsWriter(fout, name="lib2", unit=2e-3, precision=1e-5)
        for c in lib.cells.values():
            writer2.write_cell(c)
        writer2.close()
    with open(fname2, "rb") as fin:
        lib2 = gdspy.GdsLibrary()
        lib2.read_gds(fin)
    assert lib2.name == "lib2"
    assert len(lib2.cells) == 4


def hash_file(filepath):
    md5 = hashlib.md5()
    with open(filepath, mode='rb') as f:
        content = f.read()
        md5.update(content)
    return md5.hexdigest()


def test_time_changes_gds_hash(tmpdir):
    fn1 = str(tmpdir.join('nofreeze1.gds'))
    fn2 = str(tmpdir.join('nofreeze2.gds'))
    date1 = datetime.datetime(1988, 8, 28)
    date2 = datetime.datetime(2020, 12, 25)
    lib = gdspy.GdsLibrary(name='speedy')
    lib.write_gds(fn1, timestamp=date1)
    hash1 = hash_file(fn1)
    lib.write_gds(fn2, timestamp=date2)
    hash2 = hash_file(fn2)

    assert hash1 != hash2


def test_frozen_gds_has_constant_hash(tmpdir):
    fn1 = str(tmpdir.join('freeze1.gds'))
    fn2 = str(tmpdir.join('freeze2.gds'))
    frozen_date = datetime.datetime(1988, 8, 28)
    lib = gdspy.GdsLibrary(name='Elsa')
    lib.write_gds(fn1, timestamp=frozen_date)
    hash1 = hash_file(fn1)
    lib.write_gds(fn2, timestamp=frozen_date)
    hash2 = hash_file(fn2)

    assert hash1 == hash2


def test_frozen_gds_with_cell_has_constant_hash(tmpdir):
    fn1 = str(tmpdir.join('freezec1.gds'))
    fn2 = str(tmpdir.join('freezec2.gds'))
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


def test_frozen_gds_with_cell_array_has_constant_hash(tmpdir):
    fn1 = str(tmpdir.join('freezea1.gds'))
    fn2 = str(tmpdir.join('freezea2.gds'))
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