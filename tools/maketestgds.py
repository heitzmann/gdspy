######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import numpy
import gdspy


cell = gdspy.Cell('PolygonSet')
p = gdspy.PolygonSet([
    [(10, 0), (11, 0), (10, 1)],
    [(11, 0), (10, 1), (11, 1)],
    [(11, 1), (12, 1), (11, 2)],
], 1, 2)
cell.add(p)

cell = gdspy.Cell('PolygonSet_fillet')
orig = gdspy.PolygonSet([
    [(0, 0), (-1, 0), (0, -1), (0.5, -0.5), (1, 0), (1, 1), (4, -1), (1, 3), (1, 2), (0, 1)],
    [(2, -1), (3, -1), (2.5, -2)],
])
orig.datatypes=[0, 1]
p = gdspy.copy(orig, 0, 5)
p.layers = [1, 1]
p.fillet(0.3, max_points=0)
cell.add(p)
p = gdspy.copy(orig, 5, 5)
p.layers = [2, 2]
p.fillet([0.3, 0.2, 0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.1, 0.2, 0], max_points=0)
cell.add(p)
p = gdspy.copy(orig, 5, 0)
p.layers = [3, 3]
p.fillet([[0.1, 0.1, 0.4, 0, 0.4, 0.1, 0.1, 0.4, 0.4, 0.1], [0.2, 0.2, 0.5]], max_points=0)
cell.add(p)
p = gdspy.copy(orig, 0, 0)
p.layers = [4, 4]
p.fillet([0.8, [10.0, 10.0, 20.0]], max_points=199, precision=1e-6)
cell.add(p)
gdspy.write_gds('tests/test.gds', unit=1, precision=1e-7)
gdspy.LayoutViewer()
