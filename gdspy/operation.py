######################################################################
#                                                                    #
#  Copyright 2009 Lucas Heitzmann Gabrielli.                         #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import

import sys

if sys.version_info.major < 3:
    from builtins import zip
    from builtins import open
    from builtins import int
    from builtins import round
    from builtins import range
    from builtins import super

    from future import standard_library

    standard_library.install_aliases()
else:
    # Python 3 doesn't have basestring, as unicode is type string
    # Python 2 doesn't equate unicode to string, but both are basestring
    # Now isinstance(s, basestring) will be True for any python version
    basestring = str

import copy as libcopy
import numpy
from gdspy import clipper
from gdspy.polygon import PolygonSet
from gdspy.path import FlexPath, RobustPath
from gdspy.library import CellReference, CellArray


def _gather_polys(args):
    """
    Gather polygons from different argument types into a list.

    Parameters
    ----------
    args : None, `PolygonSet`, `CellReference`, `CellArray` or iterable
        Polygon types.  If this is an iterable, each element must be a
        `PolygonSet`, `CellReference`, `CellArray`, or an
        array-like[N][2] of vertices of a polygon.

    Returns
    -------
    out : list of numpy array[N][2]
        List of polygons.
    """
    if args is None:
        return []
    if isinstance(args, PolygonSet):
        return [p for p in args.polygons]
    if (
        isinstance(args, RobustPath)
        or isinstance(args, FlexPath)
        or isinstance(args, CellReference)
        or isinstance(args, CellArray)
    ):
        return args.get_polygons()
    polys = []
    for p in args:
        if isinstance(p, PolygonSet):
            polys.extend(p.polygons)
        elif (
            isinstance(p, RobustPath)
            or isinstance(args, FlexPath)
            or isinstance(p, CellReference)
            or isinstance(p, CellArray)
        ):
            polys.extend(p.get_polygons())
        else:
            polys.append(p)
    return polys


def slice(polygons, position, axis, precision=1e-3, layer=0, datatype=0):
    """
    Slice polygons and polygon sets at given positions along an axis.

    Parameters
    ----------
    polygons : `PolygonSet`, `CellReference`, `CellArray` or iterable
        Operand of the slice operation.  If this is an iterable, each
        element must be a `PolygonSet`, `CellReference`, `CellArray`,
        or an array-like[N][2] of vertices of a polygon.
    position : number or list of numbers
        Positions to perform the slicing operation along the specified
        axis.
    axis : 0 or 1
        Axis along which the polygon will be sliced.
    precision : float
        Desired precision for rounding vertice coordinates.
    layer : integer, list
        The GDSII layer numbers for the elements between each division.
        If the number of layers in the list is less than the number of
        divided regions, the list is repeated.
    datatype : integer, list
        The GDSII datatype for the resulting element (between 0 and
        255).  If the number of datatypes in the list is less than the
        number of divided regions, the list is repeated.

    Returns
    -------
    out : list[N] of `PolygonSet` or None
        Result of the slicing operation, with N = len(positions) + 1.
        Each PolygonSet comprises all polygons between 2 adjacent
        slicing positions, in crescent order.

    Examples
    --------
    >>> ring = gdspy.Round((0, 0), 10, inner_radius = 5)
    >>> result = gdspy.slice(ring, [-7, 7], 0)
    >>> cell.add(result[1])
    """
    polys = _gather_polys(polygons)
    if not isinstance(layer, list):
        layer = [layer]
    if not isinstance(datatype, list):
        datatype = [datatype]
    if not isinstance(position, list):
        pos = [position]
    else:
        pos = sorted(position)
    result = [[] for _ in range(len(pos) + 1)]
    scaling = 1 / precision
    for pol in polys:
        for r, p in zip(result, clipper._chop(pol, pos, axis, scaling)):
            r.extend(p)
    for i in range(len(result)):
        if len(result[i]) == 0:
            result[i] = None
        else:
            result[i] = PolygonSet(
                result[i], layer[i % len(layer)], datatype[i % len(datatype)]
            )
    return result


def offset(
    polygons,
    distance,
    join="miter",
    tolerance=2,
    precision=0.001,
    join_first=False,
    max_points=199,
    layer=0,
    datatype=0,
):
    """
    Shrink or expand a polygon or polygon set.

    Parameters
    ----------
    polygons : `PolygonSet`, `CellReference`, `CellArray` or iterable
        Polygons to be offset.  If this is an iterable, each element
        must be a `PolygonSet`, `CellReference`, `CellArray`, or an
        array-like[N][2] of vertices of a polygon.
    distance : number
        Offset distance.  Positive to expand, negative to shrink.
    join : 'miter', 'bevel', 'round'
        Type of join used to create the offset polygon.
    tolerance : number
        For miter joints, this number must be at least 2 and it
        represents the maximal distance in multiples of offset between
        new vertices and their original position before beveling to
        avoid spikes at acute joints.  For round joints, it indicates
        the curvature resolution in number of points per full circle.
    precision : float
        Desired precision for rounding vertex coordinates.
    join_first : bool
        Join all paths before offsetting to avoid unnecessary joins in
        adjacent polygon sides.
    max_points : integer
        If greater than 4, fracture the resulting polygons to ensure
        they have at most `max_points` vertices.  This is not a
        tessellating function, so this number should be as high as
        possible.  For example, it should be set to 199 for polygons
        being drawn in GDSII files.
    layer : integer
        The GDSII layer number for the resulting element.
    datatype : integer
        The GDSII datatype for the resulting element (between 0 and
        255).

    Returns
    -------
    out : `PolygonSet` or None
        Return the offset shape as a set of polygons.
    """
    result = clipper.offset(
        _gather_polys(polygons),
        distance,
        join,
        tolerance,
        1 / precision,
        1 if join_first else 0,
    )
    if len(result) == 0:
        return None
    return PolygonSet(result, layer, datatype).fracture(max_points, precision)


def boolean(
    operand1, operand2, operation, precision=0.001, max_points=199, layer=0, datatype=0
):
    """
    Execute any boolean operation between 2 polygons or polygon sets.

    Parameters
    ----------
    operand1 : `PolygonSet`, `CellReference`, `CellArray` or iterable
        First operand.  If this is an iterable, each element must be a
        `PolygonSet`, `CellReference`, `CellArray`, or an
        array-like[N][2] of vertices of a polygon.
    operand2 : None, `PolygonSet`, `CellReference`, `CellArray` or iterable
        Second operand.  If this is an iterable, each element must be a
        `PolygonSet`, `CellReference`, `CellArray`, or an
        array-like[N][2] of vertices of a polygon.
    operation : {'or', 'and', 'xor', 'not'}
        Boolean operation to be executed.  The 'not' operation returns
        the difference ``operand1 - operand2``.
    precision : float
        Desired precision for rounding vertice coordinates.
    max_points : integer
        If greater than 4, fracture the resulting polygons to ensure
        they have at most `max_points` vertices.  This is not a
        tessellating function, so this number should be as high as
        possible.  For example, it should be set to 199 for polygons
        being drawn in GDSII files.
    layer : integer
        The GDSII layer number for the resulting element.
    datatype : integer
        The GDSII datatype for the resulting element (between 0 and
        255).

    Returns
    -------
    out : PolygonSet or None
        Result of the boolean operation.
    """
    poly1 = _gather_polys(operand1)
    poly2 = _gather_polys(operand2)
    if len(poly2) == 0:
        if operation in ["not", "xor"]:
            if len(poly1) == 0:
                return None
            return PolygonSet(poly1, layer, datatype).fracture(max_points, precision)
        poly2.append(poly1.pop())
    result = clipper.clip(poly1, poly2, operation, 1 / precision)
    if len(result) == 0:
        return None
    return PolygonSet(result, layer, datatype).fracture(max_points, precision)


def inside(points, polygons, short_circuit="any", precision=0.001):
    """
    Test whether each of the points is within the given set of polygons.

    Parameters
    ----------
    points : array-like[N][2] or sequence of array-like[N][2]
        Coordinates of the points to be tested or groups of points to be
        tested together.
    polygons : `PolygonSet`, `CellReference`, `CellArray` or iterable
        Polygons to be tested against.  If this is an iterable, each
        element must be a `PolygonSet`, `CellReference`, `CellArray`,
        or an array-like[N][2] of vertices of a polygon.
    short_circuit : {'any', 'all'}
        If `points` is a sequence of point groups, testing within each
        group will be short-circuited if any of the points in the group
        is inside ('any') or outside ('all') the polygons.  If `points`
        is simply a sequence of points, this parameter has no effect.
    precision : float
        Desired precision for rounding vertice coordinates.

    Returns
    -------
    out : tuple
        Tuple of booleans indicating if each of the points or point
        groups is inside the set of polygons.
    """
    polys = _gather_polys(polygons)
    if numpy.isscalar(points[0][0]):
        pts = (points,)
        sc = 0
    else:
        pts = points
        sc = 1 if short_circuit == "any" else -1
    return clipper.inside(pts, polys, sc, 1 / precision)


def copy(obj, dx=0, dy=0):
    """
    Create a copy of `obj` and translate it by (dx, dy).

    Parameters
    ----------
    obj : translatable object
        Object to be copied.
    dx : number
        Distance to move in the x-direction.
    dy : number
        Distance to move in the y-direction.


    Returns
    -------
    out : translatable object
        Translated copy of original `obj`

    Examples
    --------
    >>> rectangle = gdspy.Rectangle((0, 0), (10, 20))
    >>> rectangle2 = gdspy.copy(rectangle, 2,0)
    >>> myCell.add(rectangle)
    >>> myCell.add(rectangle2)
    """

    newObj = libcopy.deepcopy(obj)
    if dx != 0 or dy != 0:
        newObj.translate(dx, dy)
    return newObj
