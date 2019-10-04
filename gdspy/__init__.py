######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################
"""
gdspy is a Python module that allows the creation of GDSII stream files.

Most features of the GDSII format are implemented, including support for
polygons with any number of vertices.

GDSII format references:

- http://boolean.klaasholwerda.nl/interface/bnf/gdsformat.html
- http://www.artwork.com/gdsii/gdsii/
- http://www.buchanan1.net/stream_description.html
"""

import warnings

__version__ = "2.0.0"

from gdspy.library import (
    Cell,
    CellReference,
    CellArray,
    GdsLibrary,
    GdsWriter,
    get_gds_units,
    get_binary_cells,
)
from gdspy.curve import Curve
from gdspy.label import Label
from gdspy.path import FlexPath, RobustPath
from gdspy.polygon import PolygonSet, Polygon, Rectangle, Round, Text, Path
from gdspy.gdsiiformat import gdsii_hash
from gdspy.operation import slice, offset, boolean, inside, copy

try:
    from gdspy.viewer import LayoutViewer
except ImportError as e:
    warnings.warn(
        "[GDSPY] LayoutViewer not available: " + str(e),
        category=ImportWarning,
        stacklevel=2,
    )
