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


cell = gdspy.Cell('polygonset')
p = gdspy.PolygonSet([
    [(1, 1), (1, 2), (5, 3)],
    [(0, 0.5), (0, 3), (-4, 1), (0, -1), (4, 0.5)],
], 8, 9)
cell.add(p)

gdspy.write_gds('tests/test.gds')
gdspy.LayoutViewer()
