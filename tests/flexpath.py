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


def test_flexpath_warnings():
    f = lambda *args: None
    for ends in ['smooth', f]:
        with pytest.warns(UserWarning):
            gdspy.FlexPath([(0, 0)], 1, ends=ends, gdsii_path=True)
    for corners in ['miter', 'bevel', 'round', 'smooth', f]:
        with pytest.warns(UserWarning):
            gdspy.FlexPath([(0, 0)], 1, corners=corners, gdsii_path=True)

def test_flexpath1(target):
    cell = gdspy.Cell('test', True)
    fp = gdspy.FlexPath([(0, 0), (1, 1)], 0.1, layer=[1])
    cell.add(fp)
    fp = gdspy.FlexPath([(1, 0), (2, 1)], 0.1, [-0.1, 0.1], ends=['round', 'extended'], tolerance=1e-3, layer=[2, 3])
    cell.add(fp)
    fp = gdspy.FlexPath([(2, 0), (3, 1)], [0.1, 0.2], 0.2, ends=(0.2, 0.1), layer=4, datatype=[1, 1])
    cell.add(fp)
    fp = gdspy.FlexPath([(3, 0), (4, 1)], [0.1, 0.2], [-0.1, 0.1], ends=[(0.2, 0.1), 'smooth'], tolerance=1e-3, datatype=5)
    cell.add(fp)
    assertsame(target['FlexPath1'], cell)
