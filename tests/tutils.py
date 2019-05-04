######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import os
import pytest
import gdspy


@pytest.fixture
def target():
    return gdspy.GdsLibrary(infile='tests' + os.sep + 'test.gds').cell_dict

def assertsame(c1, c2, tolerance=1e-6):
    d1 = c1.get_polygons(by_spec=True)
    d2 = c2.get_polygons(by_spec=True)
    for key in d1:
        assert key in d2
        result = gdspy.boolean(d1[key], d2[key], 'xor', precision=1e-7,
                               layer=key[0], datatype=100)
        if result is not None:
            r1 = gdspy.boolean(d1[key], gdspy.offset(d2[key], tolerance, precision=1e-7),
                               'not', precision=1e-7, layer=key[0], datatype=99)
            r2 = gdspy.boolean(d2[key], gdspy.offset(d1[key], tolerance, precision=1e-7),
                               'not', precision=1e-7, layer=key[0], datatype=99)
            #if not (r1 is None and r2 is None):
            #    c1.add(result)
            #    c2.add(result)
            #    if r1 is not None:
            #        c1.add(r1)
            #    if r2 is not None:
            #        c2.add(r2)
            #    gdspy.LayoutViewer(cells=[c1, c2])
            assert r1 is None
            assert r2 is None
        else:
            assert result is None
