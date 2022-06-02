######################################################################
#                                                                    #
#  Copyright 2009 Lucas Heitzmann Gabrielli.                         #
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

__version__ = "1.6.12"

import warnings

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
from gdspy.polygon import (
    PolygonSet,
    Polygon,
    Rectangle,
    Round,
    Text,
    Path,
    PolyPath,
    L1Path,
)
from gdspy.gdsiiformat import gdsii_hash, set_gdsii_timestamp
from gdspy.operation import slice, offset, boolean, inside, copy

try:
    from gdspy.viewer import LayoutViewer
except ImportError as e:
    warnings.warn(
        "[GDSPY] LayoutViewer not available: " + str(e),
        category=ImportWarning,
        stacklevel=1,
    )


def fast_boolean(*args, **kwargs):
    """
    .. deprecated:: 1.5
       `fast_boolean` is deprecated in favor of `boolean` and will be
        removed in a future version of Gdspy.
    """
    warnings.warn(
        "[GDSPY] fast_boolean is deprecated.  Use boolean instead.",
        category=DeprecationWarning,
        stacklevel=1,
    )
    return boolean(*args, **kwargs)


def write_gds(
    outfile,
    cells=None,
    name="library",
    unit=1.0e-6,
    precision=1.0e-9,
    timestamp=None,
    binary_cells=None,
):
    warnings.warn(
        "[GDSPY] write_gds and the global library is deprecated.  "
        "Use GdsLibrary.write_gds instead.",
        category=DeprecationWarning,
        stacklevel=1,
    )
    current_library.name = name
    current_library.unit = unit
    current_library.precision = precision
    current_library.write_gds(outfile, cells, timestamp, binary_cells)


current_library = GdsLibrary()
