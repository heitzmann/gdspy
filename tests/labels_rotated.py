######################################################################
#                                                                    #
#  Copyright 2009 Lucas Heitzmann Gabrielli.                         #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import gdspy

gdspy.library.use_current_library = False


def test_add_label_rotation_x_reflection():
    # add a child cell with a label
    c1 = gdspy.Cell("child")
    label = gdspy.Label("label", (0, 0))
    p = gdspy.Polygon(((0, 0), (1, 0), (0, 1)))
    c1.add(p)
    c1.add(label)

    # add parent with rotated cell
    rotation = 90
    dx = 5
    dy = 10
    c2 = gdspy.Cell("parent")
    c1ref = gdspy.CellReference(ref_cell=c1, rotation=rotation, x_reflection=True)
    c1ref.translate(dx, dy)
    c2.add(c1ref)

    label2 = c2.get_labels(set_transform=True)[0]
    label_rotation = label2.rotation
    label_x_reflection = label2.x_reflection

    assert label_rotation == rotation, [label_rotation, rotation]
    assert label_x_reflection == True, label_x_reflection
    assert label2.position[0] == dx
    assert label2.position[1] == dy


def test_add_label_rotation():
    # add a child cell with a label
    c1 = gdspy.Cell("child")
    label = gdspy.Label("label", (0, 0))
    p = gdspy.Polygon(((0, 0), (1, 0), (0, 1)))
    c1.add(p)
    c1.add(label)

    # add parent with rotated cell
    rotation = 90
    c2 = gdspy.Cell("parent")
    c1ref = gdspy.CellReference(ref_cell=c1, rotation=rotation)
    c2.add(c1ref)

    label_rotation = c2.get_labels(set_transform=True)[0].rotation
    assert label_rotation == rotation, [label_rotation, rotation]


if __name__ == "__main__":
    # add a child cell with a label
    c1 = gdspy.Cell("child")
    label = gdspy.Label("label", (0, 0))
    p = gdspy.Polygon(((0, 0), (1, 0), (0, 1)))
    c1.add(p)
    c1.add(label)

    # add parent with rotated cell
    rotation = 90
    dx = 5
    dy = 10
    c2 = gdspy.Cell("parent")
    c1ref = gdspy.CellReference(ref_cell=c1, rotation=rotation, x_reflection=True)
    c1ref.translate(dx, dy)
    c2.add(c1ref)

    label2 = c2.get_labels(set_transform=True)[0]
    label_rotation = label2.rotation
    label_x_reflection = label2.x_reflection

    assert label_rotation == rotation, [label_rotation, rotation]
    assert label_x_reflection == True, label_x_reflection
    assert label2.position[0] == dx
    assert label2.position[1] == dy

    gdspy.write_gds("a.gds", cells=[c1, c2])

    # import gdsfactory as gf
    # gf.show("a.gds")
