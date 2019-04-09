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


cell = gdspy.Cell('PolygonSet0')
p = gdspy.PolygonSet([
    [(10, 0), (11, 0), (10, 1)],
    [(11, 0), (10, 1), (11, 1)],
    [(11, 1), (12, 1), (11, 2)],
], 1, 2)
cell.add(p)

cell = gdspy.Cell('PolygonSet1')
p = gdspy.PolygonSet([
    [(10, 0), (11, 0), (10, 1)],
    [(11, 0), (10, 1), (11, 1)],
    [(11, 1), (12, 1), (11, 2)],
], 2, 3)
p.rotate(numpy.pi / 3, (10, 1))
cell.add(p)

cell = gdspy.Cell('PolygonSet2')
p = gdspy.PolygonSet([
    [(10, 0), (11, 0), (10, 1)],
    [(11, 0), (10, 1), (11, 1)],
    [(11, 1), (12, 1), (11, 2)],
], 3, 4)
p.scale(0.5)
p.scale(1, 2, center=(5, 1))
cell.add(p)

gdspy.write_gds('tests/test.gds', unit=1, precision=1e-7)
gdspy.LayoutViewer()
