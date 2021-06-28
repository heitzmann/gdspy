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

import struct
import itertools
import warnings
import numpy

from gdspy import clipper
from gdspy.hobby import _hobby

_directions_dict = {"+x": 0, "+y": 0.5, "-x": 1, "-y": -0.5}
_directions_list = ["+x", "+y", "-x", "-y"]
_halfpi = 0.5 * numpy.pi
_mpone = numpy.array((-1.0, 1.0))


class PolygonSet(object):
    """
    Set of polygonal objects.

    Parameters
    ----------
    polygons : iterable of array-like[N][2]
        List containing the coordinates of the vertices of each polygon.
    layer : integer
        The GDSII layer number for this element.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).

    Attributes
    ----------
    polygons : list of numpy array[N][2]
        Coordinates of the vertices of each polygon.
    layers : list of integer
        The GDSII layer number for each element.
    datatypes : list of integer
        The GDSII datatype for each element (between 0 and 255).
    properties : {integer: string} dictionary
        Properties for these elements.

    Notes
    -----
    The last point should not be equal to the first (polygons are
    automatically closed).

    The original GDSII specification supports only a maximum of 199
    vertices per polygon.
    """

    __slots__ = "layers", "datatypes", "polygons", "properties"

    def __init__(self, polygons, layer=0, datatype=0):
        self.polygons = [numpy.array(p) for p in polygons]
        self.layers = [layer] * len(self.polygons)
        self.datatypes = [datatype] * len(self.polygons)
        self.properties = {}

    def __str__(self):
        return (
            "PolygonSet ({} polygons, {} vertices, layers {}, datatypes {})"
        ).format(
            len(self.polygons),
            sum([len(p) for p in self.polygons]),
            list(set(self.layers)),
            list(set(self.datatypes)),
        )

    def get_bounding_box(self):
        """
        Calculate the bounding box of the polygons.

        Returns
        -------
        out : Numpy array[2, 2] or None
            Bounding box of this polygon in the form [[x_min, y_min],
            [x_max, y_max]], or None if the polygon is empty.
        """
        if len(self.polygons) == 0:
            return None
        return numpy.array(
            (
                (
                    min(pts[:, 0].min() for pts in self.polygons),
                    min(pts[:, 1].min() for pts in self.polygons),
                ),
                (
                    max(pts[:, 0].max() for pts in self.polygons),
                    max(pts[:, 1].max() for pts in self.polygons),
                ),
            )
        )

    def rotate(self, angle, center=(0, 0)):
        """
        Rotate this object.

        Parameters
        ----------
        angle : number
            The angle of rotation (in *radians*).
        center : array-like[2]
            Center point for the rotation.

        Returns
        -------
        out : `PolygonSet`
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle) * _mpone
        c0 = numpy.array(center)
        new_polys = []
        for points in self.polygons:
            pts = points - c0
            new_polys.append(pts * ca + pts[:, ::-1] * sa + c0)
        self.polygons = new_polys
        return self

    def scale(self, scalex, scaley=None, center=(0, 0)):
        """
        Scale this object.

        Parameters
        ----------
        scalex : number
            Scaling factor along the first axis.
        scaley : number or None
            Scaling factor along the second axis.  If None, same as
            `scalex`.
        center : array-like[2]
            Center point for the scaling operation.

        Returns
        -------
        out : `PolygonSet`
            This object.
        """
        c0 = numpy.array(center)
        s = scalex if scaley is None else numpy.array((scalex, scaley))
        self.polygons = [(points - c0) * s + c0 for points in self.polygons]
        return self

    def to_gds(self, outfile, multiplier):
        """
        Convert this object to a series of GDSII elements.

        Parameters
        ----------
        outfile : open file
            Output to write the GDSII.
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            elements.
        """
        for ii in range(len(self.polygons)):
            if len(self.polygons[ii]) > 8190:
                warnings.warn(
                    "[GDSPY] Polygons with more than 8190 are not supported by the "
                    "official GDSII specification.  This GDSII file might not be "
                    "compatible with all readers.",
                    stacklevel=4,
                )
                outfile.write(
                    struct.pack(
                        ">4Hh2Hh",
                        4,
                        0x0800,
                        6,
                        0x0D02,
                        self.layers[ii],
                        6,
                        0x0E02,
                        self.datatypes[ii],
                    )
                )
                xy = numpy.empty((self.polygons[ii].shape[0] + 1, 2), dtype=">i4")
                xy[:-1] = numpy.round(self.polygons[ii] * multiplier)
                xy[-1] = xy[0]
                i0 = 0
                while i0 < xy.shape[0]:
                    i1 = min(i0 + 8190, xy.shape[0])
                    outfile.write(struct.pack(">2H", 4 + 8 * (i1 - i0), 0x1003))
                    outfile.write(xy[i0:i1].tobytes())
                    i0 = i1
            else:
                outfile.write(
                    struct.pack(
                        ">4Hh2Hh2H",
                        4,
                        0x0800,
                        6,
                        0x0D02,
                        self.layers[ii],
                        6,
                        0x0E02,
                        self.datatypes[ii],
                        12 + 8 * len(self.polygons[ii]),
                        0x1003,
                    )
                )
                xy = numpy.round(self.polygons[ii] * multiplier).astype(">i4")
                outfile.write(xy.tobytes())
                outfile.write(xy[0].tobytes())
            if self.properties is not None and len(self.properties) > 0:
                size = 0
                for attr, value in self.properties.items():
                    if len(value) % 2 != 0:
                        value = value + "\0"
                    outfile.write(
                        struct.pack(">5H", 6, 0x2B02, attr, 4 + len(value), 0x2C06)
                    )
                    outfile.write(value.encode("ascii"))
                    size += len(value) + 2
                if size > 128:
                    warnings.warn(
                        "[GDSPY] Properties with size larger than 128 bytes are not "
                        "officially supported by the GDSII specification.  This file "
                        "might not be compatible with all readers.",
                        stacklevel=4,
                    )
            outfile.write(struct.pack(">2H", 4, 0x1100))

    def to_svg(self, outfile, scaling, precision):
        """
        Write an SVG fragment representation of this object.

        Parameters
        ----------
        outfile : open file
            Output to write the SVG representation.
        scaling : number
            Scaling factor for the geometry.
        precision : positive integer or `None`
            Maximal number of digits for coordinates after scaling.
        """
        for p, l, d in zip(self.polygons, self.layers, self.datatypes):
            outfile.write('<polygon class="l{}d{}" points="'.format(l, d))
            outfile.write(
                " ".join(
                    ",".join(
                        (
                            numpy.format_float_positional(
                                pt[0], trim="0", precision=precision
                            ),
                            numpy.format_float_positional(
                                pt[1], trim="0", precision=precision
                            ),
                        )
                    )
                    for pt in scaling * p
                )
            )
            outfile.write('"/>\n')

    def area(self, by_spec=False):
        """
        Calculate the total area of this polygon set.

        Parameters
        ----------
        by_spec : bool
            If True, the return value is a dictionary with
            ``{(layer, datatype): area}``.

        Returns
        -------
        out : number, dictionary
            Area of this object.
        """
        if by_spec:
            path_area = {}
            for poly, key in zip(self.polygons, zip(self.layers, self.datatypes)):
                poly_area = 0
                for ii in range(1, len(poly) - 1):
                    poly_area += (poly[0][0] - poly[ii + 1][0]) * (
                        poly[ii][1] - poly[0][1]
                    ) - (poly[0][1] - poly[ii + 1][1]) * (poly[ii][0] - poly[0][0])
                if key in path_area:
                    path_area[key] += 0.5 * abs(poly_area)
                else:
                    path_area[key] = 0.5 * abs(poly_area)
        else:
            path_area = 0
            for points in self.polygons:
                poly_area = 0
                for ii in range(1, len(points) - 1):
                    poly_area += (points[0][0] - points[ii + 1][0]) * (
                        points[ii][1] - points[0][1]
                    ) - (points[0][1] - points[ii + 1][1]) * (
                        points[ii][0] - points[0][0]
                    )
                path_area += 0.5 * abs(poly_area)
        return path_area

    def fracture(self, max_points=199, precision=1e-3):
        """
        Slice these polygons in the horizontal and vertical directions
        so that each resulting piece has at most `max_points`.  This
        operation occurs in place.

        Parameters
        ----------
        max_points : integer
            Maximal number of points in each resulting polygon (at least
            5 for the fracture to occur).
        precision : float
            Desired precision for rounding vertice coordinates.

        Returns
        -------
        out : `PolygonSet`
            This object.
        """
        if max_points > 4:
            ii = 0
            while ii < len(self.polygons):
                if len(self.polygons[ii]) > max_points:
                    pts0 = sorted(self.polygons[ii][:, 0])
                    pts1 = sorted(self.polygons[ii][:, 1])
                    ncuts = len(pts0) // max_points
                    if pts0[-1] - pts0[0] > pts1[-1] - pts1[0]:
                        # Vertical cuts
                        cuts = [
                            pts0[int(i * len(pts0) / (ncuts + 1.0) + 0.5)]
                            for i in range(1, ncuts + 1)
                        ]
                        chopped = clipper._chop(
                            self.polygons[ii], cuts, 0, 1 / precision
                        )
                    else:
                        # Horizontal cuts
                        cuts = [
                            pts1[int(i * len(pts1) / (ncuts + 1.0) + 0.5)]
                            for i in range(1, ncuts + 1)
                        ]
                        chopped = clipper._chop(
                            self.polygons[ii], cuts, 1, 1 / precision
                        )
                    self.polygons.pop(ii)
                    layer = self.layers.pop(ii)
                    datatype = self.datatypes.pop(ii)
                    self.polygons.extend(
                        numpy.array(x) for x in itertools.chain.from_iterable(chopped)
                    )
                    npols = sum(len(c) for c in chopped)
                    self.layers.extend(layer for _ in range(npols))
                    self.datatypes.extend(datatype for _ in range(npols))
                else:
                    ii += 1
        return self

    def fillet(self, radius, points_per_2pi=128, max_points=199, precision=1e-3):
        """
        Round the corners of these polygons and fractures them into
        polygons with less vertices if necessary.

        Parameters
        ----------
        radius : number, array-like
            Radius of the corners.  If number: all corners filleted by
            that amount.  If array: specify fillet radii on a
            per-polygon basis (length must be equal to the number of
            polygons in this `PolygonSet`).  Each element in the array
            can be a number (all corners filleted by the same amount) or
            another array of numbers, one per polygon vertex.
            Alternatively, the array can be flattened to have one radius
            per `PolygonSet` vertex.
        points_per_2pi : integer
            Number of vertices used to approximate a full circle.  The
            number of vertices in each corner of the polygon will be the
            fraction of this number corresponding to the angle
            encompassed by that corner with respect to 2 pi.
        max_points : integer
            Maximal number of points in each resulting polygon (at least
            5, otherwise the resulting polygon is not fractured).
        precision : float
            Desired precision for rounding vertice coordinates in case
            of fracturing.

        Returns
        -------
        out : `PolygonSet`
            This object.
        """
        two_pi = 2 * numpy.pi
        fracture = False
        if numpy.isscalar(radius):
            radii = [[radius] * p.shape[0] for p in self.polygons]
        else:
            if len(radius) == len(self.polygons):
                radii = []
                for r, p in zip(radius, self.polygons):
                    if numpy.isscalar(r):
                        radii.append([r] * p.shape[0])
                    else:
                        if len(r) != p.shape[0]:
                            raise ValueError(
                                "[GDSPY] Wrong length in fillet radius list.  "
                                "Found {} radii for polygon with {} vertices.".format(
                                    len(r), len(p.shape[0])
                                )
                            )
                        radii.append(r)
            else:
                total = sum(p.shape[0] for p in self.polygons)
                if len(radius) != total:
                    raise ValueError(
                        "[GDSPY] Wrong length in fillet radius list.  "
                        "Expected lengths are {} or {}; got {}.".format(
                            len(self.polygons), total, len(radius)
                        )
                    )
                radii = []
                n = 0
                for p in self.polygons:
                    radii.append(radius[n : n + p.shape[0]])
                    n += p.shape[0]

        for jj in range(len(self.polygons)):
            vec = self.polygons[jj].astype(float) - numpy.roll(self.polygons[jj], 1, 0)
            length = (vec[:, 0] ** 2 + vec[:, 1] ** 2) ** 0.5
            ii = numpy.flatnonzero(length)
            if len(ii) < len(length):
                self.polygons[jj] = numpy.array(self.polygons[jj][ii])
                radii[jj] = [radii[jj][i] for i in ii]
                vec = self.polygons[jj].astype(float) - numpy.roll(
                    self.polygons[jj], 1, 0
                )
                length = (vec[:, 0] ** 2 + vec[:, 1] ** 2) ** 0.5
            vec[:, 0] = vec[:, 0] / length
            vec[:, 1] = vec[:, 1] / length
            dvec = numpy.roll(vec, -1, 0) - vec
            norm = (dvec[:, 0] ** 2 + dvec[:, 1] ** 2) ** 0.5
            ii = numpy.flatnonzero(norm)
            dvec[ii, 0] = dvec[ii, 0] / norm[ii]
            dvec[ii, 1] = dvec[ii, 1] / norm[ii]
            dot = numpy.roll(vec, -1, 0) * vec
            theta = numpy.arccos(dot[:, 0] + dot[:, 1])
            ct = numpy.cos(theta * 0.5)
            tt = numpy.tan(theta * 0.5)

            new_points = []
            for ii in range(-1, len(self.polygons[jj]) - 1):
                if theta[ii] > 1e-6:
                    a0 = -vec[ii] * tt[ii] - dvec[ii] / ct[ii]
                    a0 = numpy.arctan2(a0[1], a0[0])
                    a1 = vec[ii + 1] * tt[ii] - dvec[ii] / ct[ii]
                    a1 = numpy.arctan2(a1[1], a1[0])
                    if a1 - a0 > numpy.pi:
                        a1 -= two_pi
                    elif a1 - a0 < -numpy.pi:
                        a1 += two_pi
                    n = max(
                        int(numpy.ceil(abs(a1 - a0) / two_pi * points_per_2pi) + 0.5), 2
                    )
                    a = numpy.linspace(a0, a1, n)
                    ll = radii[jj][ii] * tt[ii]
                    if ll > 0.49 * length[ii]:
                        r = 0.49 * length[ii] / tt[ii]
                        ll = 0.49 * length[ii]
                    else:
                        r = radii[jj][ii]
                    if ll > 0.49 * length[ii + 1]:
                        r = 0.49 * length[ii + 1] / tt[ii]
                    new_points.extend(
                        r * dvec[ii] / ct[ii]
                        + self.polygons[jj][ii]
                        + numpy.vstack((r * numpy.cos(a), r * numpy.sin(a))).transpose()
                    )
                else:
                    new_points.append(self.polygons[jj][ii])
            self.polygons[jj] = numpy.array(new_points)
            if len(new_points) > max_points:
                fracture = True

        if fracture:
            self.fracture(max_points, precision)
        return self

    def translate(self, dx, dy):
        """
        Translate this polygon.

        Parameters
        ----------
        dx : number
            Distance to move in the x-direction.
        dy : number
            Distance to move in the y-direction.

        Returns
        -------
        out : `PolygonSet`
            This object.
        """
        vec = numpy.array((dx, dy))
        self.polygons = [points + vec for points in self.polygons]
        return self

    def mirror(self, p1, p2=(0, 0)):
        """
        Mirror the polygons over a line through points 1 and 2

        Parameters
        ----------
        p1 : array-like[2]
            first point defining the reflection line
        p2 : array-like[2]
            second point defining the reflection line

        Returns
        -------
        out : `PolygonSet`
            This object.
        """
        origin = numpy.array(p1)
        vec = numpy.array(p2) - origin
        vec_r = vec * (2 / numpy.inner(vec, vec))
        self.polygons = [
            numpy.outer(numpy.inner(points - origin, vec_r), vec) - points + 2 * origin
            for points in self.polygons
        ]
        return self


class Polygon(PolygonSet):
    """
    Polygonal geometric object.

    Parameters
    ----------
    points : array-like[N][2]
        Coordinates of the vertices of the polygon.
    layer : integer
        The GDSII layer number for this element.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).

    Notes
    -----
    The last point should not be equal to the first (polygons are
    automatically closed).

    The original GDSII specification supports only a maximum of 199
    vertices per polygon.

    Examples
    --------
    >>> triangle_pts = [(0, 40), (15, 40), (10, 50)]
    >>> triangle = gdspy.Polygon(triangle_pts)
    >>> myCell.add(triangle)
    """

    __slots__ = "layers", "datatypes", "polygons", "properties"

    def __init__(self, points, layer=0, datatype=0):
        self.layers = [layer]
        self.datatypes = [datatype]
        self.polygons = [numpy.array(points)]
        self.properties = {}

    def __str__(self):
        return "Polygon ({} vertices, layer {}, datatype {})".format(
            len(self.polygons[0]), self.layers[0], self.datatypes[0]
        )


class Rectangle(PolygonSet):
    """
    Rectangular geometric object.

    Parameters
    ----------
    point1 : array-like[2]
        Coordinates of a corner of the rectangle.
    point2 : array-like[2]
        Coordinates of the corner of the rectangle opposite to `point1`.
    layer : integer
        The GDSII layer number for this element.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).

    Examples
    --------
    >>> rectangle = gdspy.Rectangle((0, 0), (10, 20))
    >>> myCell.add(rectangle)
    """

    __slots__ = "layers", "datatypes", "polygons", "properties"

    def __init__(self, point1, point2, layer=0, datatype=0):
        self.layers = [layer]
        self.datatypes = [datatype]
        self.polygons = [
            numpy.array(
                [
                    [point1[0], point1[1]],
                    [point1[0], point2[1]],
                    [point2[0], point2[1]],
                    [point2[0], point1[1]],
                ]
            )
        ]
        self.properties = {}

    def __str__(self):
        return (
            "Rectangle (({0[0]}, {0[1]}) to ({1[0]}, {1[1]}), layer {2}, datatype {3})"
        ).format(
            self.polygons[0][0], self.polygons[0][2], self.layers[0], self.datatypes[0]
        )

    def __repr__(self):
        return "Rectangle(({0[0]}, {0[1]}), ({1[0]}, {1[1]}), {2}, {3})".format(
            self.polygons[0][0], self.polygons[0][2], self.layers[0], self.datatypes[0]
        )


class Round(PolygonSet):
    """
    Circular geometric object.

    Represent a circle, ellipse, ring or their sections.

    Parameters
    ----------
    center : array-like[2]
        Coordinates of the center of the circle/ring.
    radius : number, array-like[2]
        Radius of the circle/outer radius of the ring.  To build an
        ellipse an array of 2 numbers can be used, representing the
        radii in the horizontal and vertical directions.
    inner_radius : number, array-like[2]
        Inner radius of the ring. To build an elliptical hole, an array
        of 2 numbers can be used, representing the radii in the
        horizontal and vertical directions.
    initial_angle : number
        Initial angle of the circular/ring section (in *radians*).
    final_angle : number
        Final angle of the circular/ring section (in *radians*).
    tolerance : float
        Approximate curvature resolution.  The number of points is
        automatically calculated.
    number_of_points : integer or None
        Manually define the number of vertices that form the object
        (polygonal approximation).  Overrides `tolerance`.
    max_points : integer
        If the number of points in the element is greater than
        `max_points`, it will be fractured in smaller polygons with
        at most `max_points` each.  If `max_points` is zero no fracture
        will occur.
    layer : integer
        The GDSII layer number for this element.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).

    Notes
    -----
    The original GDSII specification supports only a maximum of 199
    vertices per polygon.

    Examples
    --------
    >>> circle = gdspy.Round((30, 5), 8)
    >>> ell_ring = gdspy.Round((50, 5), (8, 7), inner_radius=(5, 4))
    >>> pie_slice = gdspy.Round((30, 25), 8, initial_angle=0,
    ...                             final_angle=-5.0*numpy.pi/6.0)
    >>> arc = gdspy.Round((50, 25), 8, inner_radius=5,
    ...                       initial_angle=-5.0*numpy.pi/6.0,
    ...                       final_angle=0)
    """

    __slots__ = "layers", "datatypes", "polygons", "properties"

    def __init__(
        self,
        center,
        radius,
        inner_radius=0,
        initial_angle=0,
        final_angle=0,
        tolerance=0.01,
        number_of_points=None,
        max_points=199,
        layer=0,
        datatype=0,
    ):
        if hasattr(radius, "__iter__"):
            orx, ory = radius
            radius = max(radius)

            def outer_transform(a):
                r = a - ((a + numpy.pi) % (2 * numpy.pi) - numpy.pi)
                t = numpy.arctan2(orx * numpy.sin(a), ory * numpy.cos(a)) + r
                t[a == numpy.pi] = numpy.pi
                return t

        else:
            orx = ory = radius

            def outer_transform(a):
                return a

        if hasattr(inner_radius, "__iter__"):
            irx, iry = inner_radius
            inner_radius = max(inner_radius)

            def inner_transform(a):
                r = a - ((a + numpy.pi) % (2 * numpy.pi) - numpy.pi)
                t = numpy.arctan2(irx * numpy.sin(a), iry * numpy.cos(a)) + r
                t[a == numpy.pi] = numpy.pi
                return t

        else:
            irx = iry = inner_radius

            def inner_transform(a):
                return a

        if isinstance(number_of_points, float):
            warnings.warn(
                "[GDSPY] Use of a floating number as number_of_points "
                "is deprecated in favor of tolerance.",
                category=DeprecationWarning,
                stacklevel=2,
            )
            tolerance = number_of_points
            number_of_points = None

        if number_of_points is None:
            full_angle = (
                2 * numpy.pi
                if final_angle == initial_angle
                else abs(final_angle - initial_angle)
            )
            number_of_points = max(
                3,
                1 + int(0.5 * full_angle / numpy.arccos(1 - tolerance / radius) + 0.5),
            )
            if inner_radius > 0:
                number_of_points *= 2

        pieces = (
            1
            if max_points == 0
            else int(numpy.ceil(number_of_points / float(max_points)))
        )
        number_of_points = number_of_points // pieces
        self.layers = [layer] * pieces
        self.datatypes = [datatype] * pieces
        self.polygons = [numpy.zeros((number_of_points, 2)) for _ in range(pieces)]
        self.properties = {}
        if final_angle == initial_angle and pieces > 1:
            final_angle += 2 * numpy.pi
        angles = numpy.linspace(initial_angle, final_angle, pieces + 1)
        oang = outer_transform(angles)
        iang = inner_transform(angles)
        for ii in range(pieces):
            if oang[ii + 1] == oang[ii]:
                if inner_radius <= 0:
                    t = (
                        numpy.arange(number_of_points)
                        * 2.0
                        * numpy.pi
                        / number_of_points
                    )
                    self.polygons[ii][:, 0] = numpy.cos(t) * orx + center[0]
                    self.polygons[ii][:, 1] = numpy.sin(t) * ory + center[1]
                else:
                    n2 = number_of_points // 2
                    n1 = number_of_points - n2
                    t = numpy.arange(n1) * 2.0 * numpy.pi / (n1 - 1)
                    self.polygons[ii][:n1, 0] = numpy.cos(t) * orx + center[0]
                    self.polygons[ii][:n1, 1] = numpy.sin(t) * ory + center[1]
                    t = numpy.arange(n2) * -2.0 * numpy.pi / (n2 - 1)
                    self.polygons[ii][n1:, 0] = numpy.cos(t) * irx + center[0]
                    self.polygons[ii][n1:, 1] = numpy.sin(t) * iry + center[1]
            else:
                if inner_radius <= 0:
                    t = numpy.linspace(oang[ii], oang[ii + 1], number_of_points - 1)
                    self.polygons[ii][1:, 0] = numpy.cos(t) * orx + center[0]
                    self.polygons[ii][1:, 1] = numpy.sin(t) * ory + center[1]
                    self.polygons[ii][0] += center
                else:
                    n2 = number_of_points // 2
                    n1 = number_of_points - n2
                    t = numpy.linspace(oang[ii], oang[ii + 1], n1)
                    self.polygons[ii][:n1, 0] = numpy.cos(t) * orx + center[0]
                    self.polygons[ii][:n1, 1] = numpy.sin(t) * ory + center[1]
                    t = numpy.linspace(iang[ii + 1], iang[ii], n2)
                    self.polygons[ii][n1:, 0] = numpy.cos(t) * irx + center[0]
                    self.polygons[ii][n1:, 1] = numpy.sin(t) * iry + center[1]

    def __str__(self):
        return ("Round ({} polygons, {} vertices, layers {}, datatypes {})").format(
            len(self.polygons),
            sum([len(p) for p in self.polygons]),
            list(set(self.layers)),
            list(set(self.datatypes)),
        )


class Text(PolygonSet):
    """
    Polygonal text object.

    Each letter is formed by a series of polygons.

    Parameters
    ----------
    text : string
        The text to be converted in geometric objects.
    size : number
        Height of the character.  The width of a character and the
        distance between characters are this value multiplied by 5 / 9
        and 8 / 9, respectively.  For vertical text, the distance is
        multiplied by 11 / 9.
    position : array-like[2]
        Text position (lower left corner).
    horizontal : bool
        If True, the text is written from left to right; if
        False, from top to bottom.
    angle : number
        The angle of rotation of the text.
    layer : integer
        The GDSII layer number for these elements.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).

    Examples
    --------
    >>> text = gdspy.Text('Sample text', 20, (-10, -100))
    >>> myCell.add(text)
    """

    # fmt: off
    _font = {
        '!': [[(2, 2), (3, 2), (3, 3), (2, 3)], [(2, 4), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (2, 7), (2, 6), (2, 5)]],
        '"': [[(1, 7), (2, 7), (2, 8), (2, 9), (1, 9), (1, 8)], [(3, 7), (4, 7), (4, 8), (4, 9), (3, 9), (3, 8)]],
        '#': [[(0, 3), (1, 3), (1, 2), (2, 2), (2, 3), (2, 4), (2, 5), (3, 5), (3, 4), (2, 4), (2, 3), (3, 3), (3, 2), (4, 2), (4, 3), (5, 3), (5, 4), (4, 4), (4, 5), (5, 5), (5, 6), (4, 6), (4, 7), (3, 7), (3, 6), (2, 6), (2, 7), (1, 7), (1, 6), (0, 6), (0, 5), (1, 5), (1, 4), (0, 4)]],
        '$': [[(0, 2), (1, 2), (2, 2), (2, 1), (3, 1), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (4, 4), (4, 5), (3, 5), (3, 6), (3, 7), (4, 7), (5, 7), (5, 8), (4, 8), (3, 8), (3, 9), (2, 9), (2, 8), (1, 8), (1, 7), (2, 7), (2, 6), (1, 6), (1, 5), (2, 5), (2, 4), (2, 3), (1, 3), (0, 3)], [(0, 6), (1, 6), (1, 7), (0, 7)], [(4, 3), (5, 3), (5, 4), (4, 4)]],
        '%': [[(0, 2), (1, 2), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 7), (1, 7), (2, 7), (2, 8), (2, 9), (1, 9), (0, 9), (0, 8)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (3, 4), (3, 3)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        '&': [[(0, 3), (1, 3), (1, 4), (1, 5), (0, 5), (0, 4)], [(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 2), (2, 2), (3, 2), (3, 3), (2, 3), (1, 3)], [(1, 5), (2, 5), (3, 5), (3, 6), (3, 7), (3, 8), (2, 8), (2, 7), (2, 6), (1, 6)], [(1, 8), (2, 8), (2, 9), (1, 9)], [(3, 3), (4, 3), (4, 4), (4, 5), (3, 5), (3, 4)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 5), (5, 5), (5, 6), (4, 6)]], "'": [[(2, 7), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8)]],
        '(': [[(1, 4), (2, 4), (2, 5), (2, 6), (2, 7), (1, 7), (1, 6), (1, 5)], [(2, 3), (3, 3), (3, 4), (2, 4)], [(2, 7), (3, 7), (3, 8), (2, 8)], [(3, 2), (4, 2), (4, 3), (3, 3)], [(3, 8), (4, 8), (4, 9), (3, 9)]],
        ')': [[(3, 4), (4, 4), (4, 5), (4, 6), (4, 7), (3, 7), (3, 6), (3, 5)], [(1, 2), (2, 2), (2, 3), (1, 3)], [(1, 8), (2, 8), (2, 9), (1, 9)], [(2, 3), (3, 3), (3, 4), (2, 4)], [(2, 7), (3, 7), (3, 8), (2, 8)]],
        '*': [[(0, 2), (1, 2), (1, 3), (0, 3)], [(0, 4), (1, 4), (1, 3), (2, 3), (2, 2), (3, 2), (3, 3), (4, 3), (4, 4), (5, 4), (5, 5), (4, 5), (4, 6), (3, 6), (3, 7), (2, 7), (2, 6), (1, 6), (1, 5), (0, 5)], [(0, 6), (1, 6), (1, 7), (0, 7)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        '+': [[(0, 4), (1, 4), (2, 4), (2, 3), (2, 2), (3, 2), (3, 3), (3, 4), (4, 4), (5, 4), (5, 5), (4, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6), (2, 5), (1, 5), (0, 5)]],
        ',': [[(1, 0), (2, 0), (2, 1), (1, 1)], [(2, 1), (3, 1), (3, 2), (3, 3), (2, 3), (2, 2)]],
        '-': [[(0, 4), (1, 4), (2, 4), (3, 4), (4, 4), (5, 4), (5, 5), (4, 5), (3, 5), (2, 5), (1, 5), (0, 5)]],
        '.': [[(2, 2), (3, 2), (3, 3), (2, 3)]],
        '/': [[(0, 2), (1, 2), (1, 3), (1, 4), (0, 4), (0, 3)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        '0': [[(0, 3), (1, 3), (1, 4), (2, 4), (2, 5), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 6), (4, 6), (4, 5), (4, 4), (4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (3, 7)]],
        '1': [[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (1, 8), (1, 7), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (1, 3)]],
        '2': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 7), (1, 7), (1, 8), (0, 8)], [(1, 4), (2, 4), (3, 4), (3, 5), (2, 5), (1, 5)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '3': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '4': [[(0, 4), (1, 4), (2, 4), (3, 4), (3, 3), (3, 2), (4, 2), (4, 3), (4, 4), (5, 4), (5, 5), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (3, 9), (2, 9), (2, 8), (3, 8), (3, 7), (3, 6), (3, 5), (2, 5), (1, 5), (1, 6), (0, 6), (0, 5)], [(1, 6), (2, 6), (2, 7), (2, 8), (1, 8), (1, 7)]],
        '5': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 5), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)]],
        '6': [[(0, 3), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)]],
        '7': [[(0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (4, 7), (4, 6), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (2, 5), (2, 4), (2, 3)], [(3, 5), (4, 5), (4, 6), (3, 6)]],
        '8': [[(0, 3), (1, 3), (1, 4), (1, 5), (0, 5), (0, 4)], [(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '9': [[(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 4), (4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (3, 6), (2, 6), (1, 6)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)]],
        ':': [[(2, 2), (3, 2), (3, 3), (2, 3)], [(2, 5), (3, 5), (3, 6), (2, 6)]],
        ';': [[(1, 0), (2, 0), (2, 1), (1, 1)], [(2, 1), (3, 1), (3, 2), (3, 3), (2, 3), (2, 2)], [(2, 4), (3, 4), (3, 5), (2, 5)]],
        '<': [[(0, 5), (1, 5), (1, 6), (0, 6)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(1, 6), (2, 6), (2, 7), (1, 7)], [(2, 3), (3, 3), (4, 3), (4, 4), (3, 4), (2, 4)], [(2, 7), (3, 7), (4, 7), (4, 8), (3, 8), (2, 8)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 8), (5, 8), (5, 9), (4, 9)]],
        '=': [[(0, 3), (1, 3), (2, 3), (3, 3), (4, 3), (5, 3), (5, 4), (4, 4), (3, 4), (2, 4), (1, 4), (0, 4)], [(0, 5), (1, 5), (2, 5), (3, 5), (4, 5), (5, 5), (5, 6), (4, 6), (3, 6), (2, 6), (1, 6), (0, 6)]],
        '>': [[(0, 2), (1, 2), (1, 3), (0, 3)], [(0, 8), (1, 8), (1, 9), (0, 9)], [(1, 3), (2, 3), (3, 3), (3, 4), (2, 4), (1, 4)], [(1, 7), (2, 7), (3, 7), (3, 8), (2, 8), (1, 8)], [(3, 4), (4, 4), (4, 5), (3, 5)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 5), (5, 5), (5, 6), (4, 6)]],
        '?': [[(0, 7), (1, 7), (1, 8), (0, 8)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 2), (3, 2), (3, 3), (2, 3)], [(2, 4), (3, 4), (3, 5), (2, 5)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '@': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 4), (3, 4), (4, 4), (4, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6), (2, 5)], [(4, 5), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6)]],
        'A': [[(0, 2), (1, 2), (1, 3), (1, 4), (2, 4), (3, 4), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (4, 5), (3, 5), (2, 5), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)]],
        'B': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        'C': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9)]],
        'D': [[(0, 2), (1, 2), (2, 2), (3, 2), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 7), (4, 7), (4, 8), (3, 8)], [(4, 4), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5)]],
        'E': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'F': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'G': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (3, 6), (2, 6), (2, 5), (3, 5), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9)]],
        'H': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'I': [[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (1, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (1, 3)]],
        'J': [[(0, 3), (1, 3), (1, 4), (0, 4)], [(0, 8), (1, 8), (2, 8), (3, 8), (3, 7), (3, 6), (3, 5), (3, 4), (3, 3), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(1, 2), (2, 2), (3, 2), (3, 3), (2, 3), (1, 3)]],
        'K': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (2, 6), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(2, 4), (3, 4), (3, 5), (2, 5)], [(2, 6), (3, 6), (3, 7), (2, 7)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 7), (4, 7), (4, 8), (3, 8)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 8), (5, 8), (5, 9), (4, 9)]],
        'L': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'M': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 7), (2, 8), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(2, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6)], [(3, 7), (4, 7), (4, 6), (4, 5), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (3, 8)]],
        'N': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 7), (2, 8), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(2, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6)], [(3, 4), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (4, 5), (3, 5)]],
        'O': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'P': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        'Q': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4), (3, 4), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 4), (3, 4), (3, 5), (2, 5)]],
        'R': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (3, 4), (4, 4), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (4, 3)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        'S': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6)], [(1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)]],
        'T': [[(0, 8), (1, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)]],
        'U': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'V': [[(0, 5), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6)], [(1, 3), (2, 3), (2, 4), (2, 5), (1, 5), (1, 4)], [(2, 2), (3, 2), (3, 3), (2, 3)], [(3, 3), (4, 3), (4, 4), (4, 5), (3, 5), (3, 4)], [(4, 5), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6)]],
        'W': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (2, 3), (1, 3)], [(2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (2, 6), (2, 5), (2, 4)], [(3, 2), (4, 2), (4, 3), (3, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'X': [[(0, 2), (1, 2), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(1, 6), (2, 6), (2, 7), (1, 7)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 4), (4, 4), (4, 5), (3, 5)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (4, 3)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        'Y': [[(0, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8)], [(1, 5), (2, 5), (2, 6), (2, 7), (1, 7), (1, 6)], [(2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (2, 5), (2, 4), (2, 3)], [(3, 5), (4, 5), (4, 6), (4, 7), (3, 7), (3, 6)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        'Z': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 6), (4, 6), (4, 7), (3, 7)]],
        '[': [[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (1, 8), (1, 7), (1, 6), (1, 5), (1, 4), (1, 3)]],
        '\\': [[(0, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8)], [(1, 6), (2, 6), (2, 7), (1, 7)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 4), (4, 4), (4, 5), (3, 5)], [(4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (4, 3)]],
        ']': [[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (1, 8), (2, 8), (3, 8), (3, 7), (3, 6), (3, 5), (3, 4), (3, 3), (2, 3), (1, 3)]],
        '^': [[(0, 6), (1, 6), (1, 7), (0, 7)], [(1, 7), (2, 7), (2, 8), (1, 8)], [(2, 8), (3, 8), (3, 9), (2, 9)], [(3, 7), (4, 7), (4, 8), (3, 8)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        '_': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)]],
        '`': [[(1, 8), (2, 8), (2, 9), (1, 9)], [(2, 7), (3, 7), (3, 8), (2, 8)]],
        'a': [[(0, 3), (1, 3), (1, 4), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (3, 5), (2, 5), (1, 5), (1, 4), (2, 4), (3, 4), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'b': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4)]],
        'c': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'd': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (3, 7), (2, 7), (1, 7), (1, 6), (2, 6), (3, 6), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)]],
        'e': [[(0, 3), (1, 3), (1, 4), (2, 4), (3, 4), (4, 4), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (3, 5), (2, 5), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'f': [[(0, 5), (1, 5), (1, 4), (1, 3), (1, 2), (2, 2), (2, 3), (2, 4), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (2, 7), (2, 8), (1, 8), (1, 7), (1, 6), (0, 6)], [(2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9)]],
        'g': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 0), (2, 0), (3, 0), (4, 0), (4, 1), (3, 1), (2, 1), (1, 1)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 1), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'h': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3)]],
        'i': [[(1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (2, 7), (1, 7)], [(2, 8), (3, 8), (3, 9), (2, 9)]],
        'j': [[(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1)], [(1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (2, 1), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (2, 7), (1, 7)], [(2, 8), (3, 8), (3, 9), (2, 9)]],
        'k': [[(0, 2), (1, 2), (1, 3), (1, 4), (2, 4), (3, 4), (3, 5), (2, 5), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        'l': [[(1, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (1, 9)]],
        'm': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3)]],
        'n': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3)]],
        'o': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4)]],
        'p': [[(0, 0), (1, 0), (1, 1), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3), (0, 2), (0, 1)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4)]],
        'q': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 1), (4, 0), (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'r': [[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (3, 6), (2, 6), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7)]],
        's': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 5), (1, 5), (1, 6), (0, 6)], [(1, 4), (2, 4), (3, 4), (4, 4), (4, 5), (3, 5), (2, 5), (1, 5)], [(1, 6), (2, 6), (3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (2, 7), (1, 7)], [(4, 3), (5, 3), (5, 4), (4, 4)]],
        't': [[(1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (2, 7), (1, 7)], [(3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3)]],
        'u': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'v': [[(0, 5), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6)], [(1, 3), (2, 3), (2, 4), (2, 5), (1, 5), (1, 4)], [(2, 2), (3, 2), (3, 3), (2, 3)], [(3, 3), (4, 3), (4, 4), (4, 5), (3, 5), (3, 4)], [(4, 5), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6)]],
        'w': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (2, 3), (1, 3)], [(2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (2, 6), (2, 5), (2, 4)], [(3, 2), (4, 2), (4, 3), (3, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'x': [[(0, 2), (1, 2), (1, 3), (0, 3)], [(0, 6), (1, 6), (1, 7), (0, 7)], [(1, 3), (2, 3), (2, 4), (1, 4)], [(1, 5), (2, 5), (2, 6), (1, 6)], [(2, 4), (3, 4), (3, 5), (2, 5)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        'y': [[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 0), (2, 0), (3, 0), (4, 0), (4, 1), (3, 1), (2, 1), (1, 1)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 1), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)]],
        'z': [[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (2, 4), (1, 4), (1, 3), (0, 3)], [(0, 6), (1, 6), (2, 6), (3, 6), (3, 5), (4, 5), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7)], [(2, 4), (3, 4), (3, 5), (2, 5)]],
        '{': [[(1, 5), (2, 5), (2, 4), (2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (2, 8), (2, 7), (2, 6), (1, 6)], [(3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3)], [(3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9)]],
        '|': [[(2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3)]],
        '}': [[(0, 2), (1, 2), (2, 2), (2, 3), (1, 3), (0, 3)], [(0, 8), (1, 8), (2, 8), (2, 9), (1, 9), (0, 9)], [(2, 3), (3, 3), (3, 4), (3, 5), (4, 5), (4, 6), (3, 6), (3, 7), (3, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4)]],
        '~': [[(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 8), (2, 8), (2, 9), (1, 9)], [(2, 7), (3, 7), (3, 8), (2, 8)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
    }
    # fmt: on

    __slots__ = "layers", "datatypes", "polygons", "properties"

    def __init__(
        self, text, size, position=(0, 0), horizontal=True, angle=0, layer=0, datatype=0
    ):
        self.polygons = []
        posX = 0
        posY = 0
        text_multiplier = size / 9.0
        if angle == 0:
            ca = 1
            sa = 0
        else:
            ca = numpy.cos(angle)
            sa = numpy.sin(angle)
        for jj in range(len(text)):
            if text[jj] == "\n":
                if horizontal:
                    posY -= 11
                    posX = 0
                else:
                    posX += 8
                    posY = 0
            elif text[jj] == "\t":
                if horizontal:
                    posX = posX + 32 - (posX + 8) % 32
                else:
                    posY = posY - 11 - (posY - 22) % 44
            else:
                if text[jj] in Text._font:
                    for p in Text._font[text[jj]]:
                        polygon = p[:]
                        for ii in range(len(polygon)):
                            xp = text_multiplier * (posX + polygon[ii][0])
                            yp = text_multiplier * (posY + polygon[ii][1])
                            polygon[ii] = (
                                position[0] + xp * ca - yp * sa,
                                position[1] + xp * sa + yp * ca,
                            )
                        self.polygons.append(numpy.array(polygon))
                if horizontal:
                    posX += 8
                else:
                    posY -= 11
        self.layers = [layer] * len(self.polygons)
        self.datatypes = [datatype] * len(self.polygons)
        self.properties = {}

    def __str__(self):
        return ("Text ({} polygons, {} vertices, layers {}, datatypes {})").format(
            len(self.polygons),
            sum([len(p) for p in self.polygons]),
            list(set(self.layers)),
            list(set(self.datatypes)),
        )


class Path(PolygonSet):
    """
    Series of geometric objects that form a path or a collection of
    parallel paths.

    Parameters
    ----------
    width : number
        The width of each path.
    initial_point : array-like[2]
        Starting position of the path.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number
        Distance between the centers of adjacent paths.

    Attributes
    ----------
    x : number
        Current position of the path in the x direction.
    y : number
        Current position of the path in the y direction.
    w : number
        *Half*-width of each path.
    n : integer
        Number of parallel paths.
    direction : '+x', '-x', '+y', '-y' or number
        Direction or angle (in *radians*) the path points to.
    distance : number
        Distance between the centers of adjacent paths.
    length : number
        Length of the central path axis.  If only one path is created,
        this is the real length of the path.
    properties : {integer: string} dictionary
        Properties for this path.
    """

    __slots__ = (
        "layers",
        "datatypes",
        "polygons",
        "x",
        "y",
        "w",
        "n",
        "direction",
        "distance",
        "length",
        "properties",
    )

    def __init__(self, width, initial_point=(0, 0), number_of_paths=1, distance=0):
        self.x = initial_point[0]
        self.y = initial_point[1]
        self.w = width * 0.5
        self.n = number_of_paths
        self.direction = "+x"
        self.distance = distance
        self.length = 0.0
        self.polygons = []
        self.layers = []
        self.datatypes = []
        self.properties = {}

    def __str__(self):
        if self.n > 1:
            return "Path (x{}, end at ({}, {}) towards {}, length {}, width {}, {} apart, {} polygons, {} vertices, layers {}, datatypes {})".format(
                self.n,
                self.x,
                self.y,
                self.direction,
                self.length,
                self.w * 2,
                self.distance,
                len(self.polygons),
                sum([len(p) for p in self.polygons]),
                list(set(self.layers)),
                list(set(self.datatypes)),
            )
        else:
            return "Path (end at ({}, {}) towards {}, length {}, width {}, {} polygons, {} vertices, layers {}, datatypes {})".format(
                self.x,
                self.y,
                self.direction,
                self.length,
                self.w * 2,
                len(self.polygons),
                sum([len(p) for p in self.polygons]),
                list(set(self.layers)),
                list(set(self.datatypes)),
            )

    def translate(self, dx, dy):
        """
        Translate this object.

        Parameters
        ----------
        dx : number
            Distance to move in the x-direction.
        dy : number
            Distance to move in the y-direction.

        Returns
        -------
        out : `Path`
            This object.
        """
        vec = numpy.array((dx, dy))
        self.polygons = [points + vec for points in self.polygons]
        self.x += dx
        self.y += dy
        return self

    def rotate(self, angle, center=(0, 0)):
        """
        Rotate this object.

        Parameters
        ----------
        angle : number
            The angle of rotation (in *radians*).
        center : array-like[2]
            Center point for the rotation.

        Returns
        -------
        out : `Path`
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle) * _mpone
        c0 = numpy.array(center)
        if isinstance(self.direction, basestring):
            self.direction = _directions_dict[self.direction] * numpy.pi
        self.direction += angle
        cur = numpy.array((self.x, self.y)) - c0
        self.x, self.y = cur * ca + cur[::-1] * sa + c0
        self.polygons = [
            (points - c0) * ca + (points - c0)[:, ::-1] * sa + c0
            for points in self.polygons
        ]
        return self

    def scale(self, scalex, scaley=None, center=(0, 0)):
        """
        Scale this object.

        Parameters
        ----------
        scalex : number
            Scaling factor along the first axis.
        scaley : number or None
            Scaling factor along the second axis.  If None, same as
            `scalex`.
        center : array-like[2]
            Center point for the scaling operation.

        Returns
        -------
        out : `Path`
            This object.

        Notes
        -----
        The direction of the path is not modified by this method and
        its width is scaled only by `scalex`.
        """
        c0 = numpy.array(center)
        s = scalex if scaley is None else numpy.array((scalex, scaley))
        self.polygons = [(points - c0) * s + c0 for points in self.polygons]
        self.x = (self.x - c0[0]) * scalex + c0[0]
        self.y = (self.y - c0[1]) * (scalex if scaley is None else scaley) + c0[1]
        self.w *= scalex
        return self

    def mirror(self, p1, p2=(0, 0)):
        """
        Mirror the polygons over a line through points 1 and 2

        Parameters
        ----------
        p1 : array-like[2]
            first point defining the reflection line
        p2 : array-like[2]
            second point defining the reflection line

        Returns
        -------
        out : `Path`
            This object.
        """
        origin = numpy.array(p1)
        vec = numpy.array(p2) - origin
        vec_r = vec * (2 / numpy.inner(vec, vec))
        self.polygons = [
            numpy.outer(numpy.inner(points - origin, vec_r), vec) - points + 2 * origin
            for points in self.polygons
        ]
        dot = (self.x - origin[0]) * vec_r[0] + (self.y - origin[1]) * vec_r[1]
        self.x = dot * vec[0] - self.x + 2 * origin[0]
        self.y = dot * vec[1] - self.y + 2 * origin[1]
        if isinstance(self.direction, basestring):
            self.direction = _directions_dict[self.direction] * numpy.pi
        self.direction = 2 * numpy.arctan2(vec[1], vec[0]) - self.direction
        return self

    def segment(
        self,
        length,
        direction=None,
        final_width=None,
        final_distance=None,
        axis_offset=0,
        layer=0,
        datatype=0,
    ):
        """
        Add a straight section to the path.

        Parameters
        ----------
        length : number
            Length of the section to add.
        direction : '+x', '-x', '+y', '-y' or number
            Direction or angle (in *radians*) of rotation of the
            segment.
        final_width : number
            If set, the paths of this segment will have their widths
            linearly changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from
            its current value to this one along this segment.
        axis_offset : number
            If set, the paths will be offset from their direction by
            this amount.
        layer : integer, list
            The GDSII layer numbers for the elements of each path.  If
            the number of layers in the list is less than the number
            of paths, the list is repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between
            0 and 255).  If the number of datatypes in the list is less
            than the number of paths, the list is repeated.

        Returns
        -------
        out : `Path`
            This object.
        """
        if direction is None:
            direction = self.direction
        else:
            self.direction = direction
        if direction == "+x":
            ca = 1
            sa = 0
        elif direction == "-x":
            ca = -1
            sa = 0
        elif direction == "+y":
            ca = 0
            sa = 1
        elif direction == "-y":
            ca = 0
            sa = -1
        else:
            ca = numpy.cos(direction)
            sa = numpy.sin(direction)
        old_x = self.x
        old_y = self.y
        self.x += length * ca + axis_offset * sa
        self.y += length * sa - axis_offset * ca
        old_w = self.w
        old_distance = self.distance
        if final_width is not None:
            self.w = final_width * 0.5
        if final_distance is not None:
            self.distance = final_distance
        if (self.w != 0) or (old_w != 0):
            for ii in range(self.n):
                d0 = ii * self.distance - (self.n - 1) * self.distance * 0.5
                old_d0 = ii * old_distance - (self.n - 1) * old_distance * 0.5
                self.polygons.append(
                    numpy.array(
                        [
                            (
                                old_x + (old_d0 - old_w) * sa,
                                old_y - (old_d0 - old_w) * ca,
                            ),
                            (
                                old_x + (old_d0 + old_w) * sa,
                                old_y - (old_d0 + old_w) * ca,
                            ),
                            (self.x + (d0 + self.w) * sa, self.y - (d0 + self.w) * ca),
                            (self.x + (d0 - self.w) * sa, self.y - (d0 - self.w) * ca),
                        ]
                    )
                )
                if self.w == 0:
                    self.polygons[-1] = self.polygons[-1][:-1]
                if old_w == 0:
                    self.polygons[-1] = self.polygons[-1][1:]
            self.length += (length ** 2 + axis_offset ** 2) ** 0.5
            if isinstance(layer, list):
                self.layers.extend((layer * (self.n // len(layer) + 1))[: self.n])
            else:
                self.layers.extend(layer for _ in range(self.n))
            if isinstance(datatype, list):
                self.datatypes.extend(
                    (datatype * (self.n // len(datatype) + 1))[: self.n]
                )
            else:
                self.datatypes.extend(datatype for _ in range(self.n))
        return self

    def arc(
        self,
        radius,
        initial_angle,
        final_angle,
        tolerance=0.01,
        number_of_points=None,
        max_points=199,
        final_width=None,
        final_distance=None,
        layer=0,
        datatype=0,
    ):
        """
        Add a curved section to the path.

        Parameters
        ----------
        radius : number
            Central radius of the section.
        initial_angle : number
            Initial angle of the curve (in *radians*).
        final_angle : number
            Final angle of the curve (in *radians*).
        tolerance : float
            Approximate curvature resolution.  The number of points is
            automatically calculated.
        number_of_points : integer or None
            Manually define the number of vertices that form the object
            (polygonal approximation).  Overrides `tolerance`.
        max_points : integer
            If the number of points in the element is greater than
            `max_points`, it will be fractured in smaller polygons with
            at most `max_points` each.  If `max_points` is zero no
            fracture will occur.
        final_width : number
            If set, the paths of this segment will have their widths
            linearly changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from
            its current value to this one along this segment.
        layer : integer, list
            The GDSII layer numbers for the elements of each path.  If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255).  If the number of datatypes in the list is less
            than the number of paths, the list is repeated.

        Returns
        -------
        out : `Path`
            This object.

        Notes
        -----
        The original GDSII specification supports only a maximum of 199
        vertices per polygon.
        """
        cx = self.x - radius * numpy.cos(initial_angle)
        cy = self.y - radius * numpy.sin(initial_angle)
        self.x = cx + radius * numpy.cos(final_angle)
        self.y = cy + radius * numpy.sin(final_angle)
        if final_angle > initial_angle:
            self.direction = final_angle + numpy.pi * 0.5
        else:
            self.direction = final_angle - numpy.pi * 0.5
        old_w = self.w
        old_distance = self.distance
        if final_width is not None:
            self.w = final_width * 0.5
        if final_distance is not None:
            self.distance = final_distance
        if isinstance(number_of_points, float):
            warnings.warn(
                "[GDSPY] Use of a floating number as number_of_points "
                "is deprecated in favor of tolerance.",
                category=DeprecationWarning,
                stacklevel=2,
            )
            tolerance = number_of_points
            number_of_points = None
        if number_of_points is None:
            r = (
                radius
                + max(old_distance, self.distance) * (self.n - 1) * 0.5
                + max(old_w, self.w)
            )
            number_of_points = max(
                6,
                2
                + 2
                * int(
                    0.5
                    * abs(final_angle - initial_angle)
                    / numpy.arccos(1 - tolerance / r)
                    + 0.5
                ),
            )
        pieces = (
            1
            if max_points == 0
            else int(numpy.ceil(number_of_points / float(max_points)))
        )
        number_of_points = number_of_points // pieces
        widths = numpy.linspace(old_w, self.w, pieces + 1)
        distances = numpy.linspace(old_distance, self.distance, pieces + 1)
        angles = numpy.linspace(initial_angle, final_angle, pieces + 1)
        if (self.w != 0) or (old_w != 0):
            for jj in range(pieces):
                for ii in range(self.n):
                    self.polygons.append(numpy.zeros((number_of_points, 2)))
                    r0 = (
                        radius
                        + ii * distances[jj + 1]
                        - (self.n - 1) * distances[jj + 1] * 0.5
                    )
                    old_r0 = (
                        radius + ii * distances[jj] - (self.n - 1) * distances[jj] * 0.5
                    )
                    pts2 = number_of_points // 2
                    pts1 = number_of_points - pts2
                    ang = numpy.linspace(angles[jj], angles[jj + 1], pts1)
                    rad = numpy.linspace(old_r0 + widths[jj], r0 + widths[jj + 1], pts1)
                    self.polygons[-1][:pts1, 0] = numpy.cos(ang) * rad + cx
                    self.polygons[-1][:pts1, 1] = numpy.sin(ang) * rad + cy
                    if widths[jj + 1] == 0:
                        pts1 -= 1
                        pts2 += 1
                    if widths[jj] == 0:
                        self.polygons[-1][: pts1 - 1] = numpy.array(
                            self.polygons[-1][1:pts1]
                        )
                        pts1 -= 1
                        pts2 += 1
                    ang = numpy.linspace(angles[jj + 1], angles[jj], pts2)
                    rad = numpy.linspace(r0 - widths[jj + 1], old_r0 - widths[jj], pts2)
                    if rad[0] <= 0 or rad[-1] <= 0:
                        warnings.warn(
                            "[GDSPY] Path arc with width larger than radius "
                            "created: possible self-intersecting polygon.",
                            stacklevel=2,
                        )
                    self.polygons[-1][pts1:, 0] = numpy.cos(ang) * rad + cx
                    self.polygons[-1][pts1:, 1] = numpy.sin(ang) * rad + cy
                self.length += abs((angles[jj + 1] - angles[jj]) * radius)
                if isinstance(layer, list):
                    self.layers.extend((layer * (self.n // len(layer) + 1))[: self.n])
                else:
                    self.layers.extend(layer for _ in range(self.n))
                if isinstance(datatype, list):
                    self.datatypes.extend(
                        (datatype * (self.n // len(datatype) + 1))[: self.n]
                    )
                else:
                    self.datatypes.extend(datatype for _ in range(self.n))
        return self

    def turn(
        self,
        radius,
        angle,
        tolerance=0.01,
        number_of_points=None,
        max_points=199,
        final_width=None,
        final_distance=None,
        layer=0,
        datatype=0,
    ):
        """
        Add a curved section to the path.

        Parameters
        ----------
        radius : number
            Central radius of the section.
        angle : 'r', 'l', 'rr', 'll' or number
            Angle (in *radians*) of rotation of the path.  The values
            'r' and 'l' represent 90-degree turns cw and ccw,
            respectively; the values 'rr' and 'll' represent analogous
            180-degree turns.
        tolerance : float
            Approximate curvature resolution.  The number of points is
            automatically calculated.
        number_of_points : integer or None
            Manually define the number of vertices that form the object
            (polygonal approximation).  Overrides `tolerance`.
        max_points : integer
            If the number of points in the element is greater than
            `max_points`, it will be fractured in smaller polygons with
            at most `max_points` each.  If `max_points` is zero no
            fracture will occur.
        final_width : number
            If set, the paths of this segment will have their widths
            linearly changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from
            its current value to this one along this segment.
        layer : integer, list
            The GDSII layer numbers for the elements of each path.  If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255).  If the number of datatypes in the list is less
            than the number of paths, the list is repeated.

        Returns
        -------
        out : `Path`
            This object.

        Notes
        -----
        The original GDSII specification supports only a maximum of 199
        vertices per polygon.
        """
        exact = True
        if angle == "r":
            delta_i = _halfpi
            delta_f = 0
        elif angle == "rr":
            delta_i = _halfpi
            delta_f = -delta_i
        elif angle == "l":
            delta_i = -_halfpi
            delta_f = 0
        elif angle == "ll":
            delta_i = -_halfpi
            delta_f = -delta_i
        elif angle < 0:
            exact = False
            delta_i = _halfpi
            delta_f = delta_i + angle
        else:
            exact = False
            delta_i = -_halfpi
            delta_f = delta_i + angle
        if self.direction == "+x":
            self.direction = 0
        elif self.direction == "-x":
            self.direction = numpy.pi
        elif self.direction == "+y":
            self.direction = _halfpi
        elif self.direction == "-y":
            self.direction = -_halfpi
        elif exact:
            exact = False
        self.arc(
            radius,
            self.direction + delta_i,
            self.direction + delta_f,
            tolerance,
            number_of_points,
            max_points,
            final_width,
            final_distance,
            layer,
            datatype,
        )
        if exact:
            self.direction = _directions_list[int(round(self.direction / _halfpi)) % 4]
        return self

    def parametric(
        self,
        curve_function,
        curve_derivative=None,
        tolerance=0.01,
        number_of_evaluations=5,
        max_points=199,
        final_width=None,
        final_distance=None,
        relative=True,
        layer=0,
        datatype=0,
    ):
        """
        Add a parametric curve to the path.

        `curve_function` will be evaluated uniformly in the interval
        [0, 1] at least `number_of_points` times.  More points will be
        added to the curve at the midpoint between evaluations if that
        points presents error larger than `tolerance`.

        Parameters
        ----------
        curve_function : callable
            Function that defines the curve.  Must be a function of one
            argument (that varies from 0 to 1) that returns a 2-element
            array with the coordinates of the curve.
        curve_derivative : callable
            If set, it should be the derivative of the curve function.
            Must be a function of one argument (that varies from 0 to 1)
            that returns a 2-element array.  If None, the derivative
            will be calculated numerically.
        tolerance : number
            Acceptable tolerance for the approximation of the curve
            function by a finite number of evaluations.
        number_of_evaluations : integer
            Initial number of points where the curve function will be
            evaluated.  According to `tolerance`, more evaluations will
            be performed.
        max_points : integer
            Elements will be fractured until each polygon has at most
            `max_points`.  If `max_points` is less than 4, no fracture
            will occur.
        final_width : number or function
            If set to a number, the paths of this segment will have
            their widths linearly changed from their current value to
            this one.  If set to a function, it must be a function of
            one argument (that varies from 0 to 1) and returns the width
            of the path.
        final_distance : number or function
            If set to a number, the distance between paths is linearly
            change from its current value to this one.  If set to a
            function, it must be a function of one argument (that varies
            from 0 to 1) and returns the width of the path.
        relative : bool
            If True, the return values of `curve_function` are used as
            offsets from the current path position, i.e., to ensure a
            continuous path, ``curve_function(0)`` must be (0, 0).
            Otherwise, they are used as absolute coordinates.
        layer : integer, list
            The GDSII layer numbers for the elements of each path.  If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255).  If the number of datatypes in the list is less
            than the number of paths, the list is repeated.

        Returns
        -------
        out : `Path`
            This object.

        Notes
        -----
        The norm of the vector returned by `curve_derivative` is not
        important.  Only the direction is used.

        The original GDSII specification supports only a maximum of 199
        vertices per polygon.

        Examples
        --------
        >>> def my_parametric_curve(t):
        ...         return (2**t, t**2)
        >>> def my_parametric_curve_derivative(t):
        ...         return (0.69315 * 2**t, 2 * t)
        >>> my_path.parametric(my_parametric_curve,
        ...                    my_parametric_curve_derivative)
        """
        err = tolerance ** 2
        points = list(numpy.linspace(0, 1, number_of_evaluations))
        values = [numpy.array(curve_function(u)) for u in points]
        delta = points[1]
        i = 1
        while i < len(points):
            midpoint = 0.5 * (points[i] + points[i - 1])
            midvalue = numpy.array(curve_function(midpoint))
            test_err = (values[i] + values[i - 1]) / 2 - midvalue
            if test_err[0] ** 2 + test_err[1] ** 2 > err:
                delta = min(delta, points[i] - midpoint)
                points.insert(i, midpoint)
                values.insert(i, midvalue)
            else:
                i += 1
        points = numpy.array(points)
        values = numpy.array(values)
        dvs = values[1:] - values[:-1]
        self.length += ((dvs[:, 0] ** 2 + dvs[:, 1] ** 2) ** 0.5).sum()

        delta *= 0.5
        if curve_derivative is None:
            derivs = numpy.vstack(
                (
                    numpy.array(curve_function(delta)) - values[0],
                    [
                        numpy.array(curve_function(u + delta))
                        - numpy.array(curve_function(u - delta))
                        for u in points[1:-1]
                    ],
                    values[-1] - numpy.array(curve_function(1 - delta)),
                )
            )
        else:
            derivs = numpy.array([curve_derivative(u) for u in points])

        if not callable(final_width):
            if final_width is None:
                width = numpy.full_like(points, self.w)
            else:
                width = self.w + (final_width * 0.5 - self.w) * points
                self.w = final_width * 0.5
        else:
            width = numpy.array([0.5 * final_width(u) for u in points])
            self.w = width[-1]

        if not callable(final_distance):
            if final_distance is None:
                dist = numpy.full_like(points, self.distance)
            else:
                dist = self.distance + (final_distance - self.distance) * points
                self.distance = final_distance
        else:
            dist = numpy.array([final_distance(u) for u in points])
            self.distance = dist[-1]

        np = points.shape[0]
        sh = (np, 1)
        if relative:
            x0 = values + numpy.array((self.x, self.y))
        else:
            x0 = values
        dx = (
            derivs[:, ::-1]
            * _mpone
            / ((derivs[:, 0] ** 2 + derivs[:, 1] ** 2) ** 0.5).reshape(sh)
        )
        width = width.reshape(sh)
        dist = dist.reshape(sh)

        self.x = x0[-1, 0]
        self.y = x0[-1, 1]
        self.direction = numpy.arctan2(-dx[-1, 0], dx[-1, 1])

        if max_points < 4:
            max_points = np
        else:
            max_points = max_points // 2

        i0 = 0
        while i0 < np - 1:
            i1 = min(i0 + max_points, np)
            for ii in range(self.n):
                p1 = x0[i0:i1] + dx[i0:i1] * (
                    dist[i0:i1] * (ii - (self.n - 1) * 0.5) + width[i0:i1]
                )
                p2 = (
                    x0[i0:i1]
                    + dx[i0:i1]
                    * (dist[i0:i1] * (ii - (self.n - 1) * 0.5) - width[i0:i1])
                )[::-1]
                if width[i1 - 1, 0] == 0:
                    p2 = p2[1:]
                if width[i0, 0] == 0:
                    p1 = p1[1:]
                self.polygons.append(numpy.concatenate((p1, p2)))
            if isinstance(layer, list):
                self.layers.extend((layer * (self.n // len(layer) + 1))[: self.n])
            else:
                self.layers.extend(layer for _ in range(self.n))
            if isinstance(datatype, list):
                self.datatypes.extend(
                    (datatype * (self.n // len(datatype) + 1))[: self.n]
                )
            else:
                self.datatypes.extend(datatype for _ in range(self.n))
            i0 = i1 - 1
        return self

    def bezier(
        self,
        points,
        tolerance=0.01,
        number_of_evaluations=5,
        max_points=199,
        final_width=None,
        final_distance=None,
        relative=True,
        layer=0,
        datatype=0,
    ):
        """
        Add a Bezier curve to the path.

        A Bezier curve is added to the path starting from its current
        position and finishing at the last point in the `points` array.

        Parameters
        ----------
        points : array-like[N][2]
            Control points defining the Bezier curve.
        tolerance : number
            Acceptable tolerance for the approximation of the curve
            function by a finite number of evaluations.
        number_of_evaluations : integer
            Initial number of points where the curve function will be
            evaluated.  According to `tolerance`, more evaluations will
            be performed.
        max_points : integer
            Elements will be fractured until each polygon has at most
            `max_points`.  If `max_points` is zero no fracture will
            occur.
        final_width : number or function
            If set to a number, the paths of this segment will have
            their widths linearly changed from their current value to
            this one.  If set to a function, it must be a function of
            one argument (that varies from 0 to 1) and returns the width
            of the path.
        final_distance : number or function
            If set to a number, the distance between paths is linearly
            change from its current value to this one.  If set to a
            function, it must be a function of one argument (that varies
            from 0 to 1) and returns the width of the path.
        relative : bool
            If True, all coordinates in the `points` array are used as
            offsets from the current path position, i.e., if the path is
            at (1, -2) and the last point in the array is (10, 25), the
            constructed Bezier will end at (1 + 10, -2 + 25) = (11, 23).
            Otherwise, the points are used as absolute coordinates.
        layer : integer, list
            The GDSII layer numbers for the elements of each path.  If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255).  If the number of datatypes in the list is less
            than the number of paths, the list is repeated.

        Returns
        -------
        out : `Path`
            This object.

        Notes
        -----
        The original GDSII specification supports only a maximum of 199
        vertices per polygon.
        """
        if relative:
            pts = numpy.vstack(([(0, 0)], points))
        else:
            pts = numpy.vstack(([(self.x, self.y)], points))
        dpts = (pts.shape[0] - 1) * (pts[1:] - pts[:-1])
        self.parametric(
            _func_bezier(pts),
            _func_bezier(dpts),
            tolerance,
            number_of_evaluations,
            max_points,
            final_width,
            final_distance,
            relative,
            layer,
            datatype,
        )
        return self

    def smooth(
        self,
        points,
        angles=None,
        curl_start=1,
        curl_end=1,
        t_in=1,
        t_out=1,
        cycle=False,
        tolerance=0.01,
        number_of_evaluations=5,
        max_points=199,
        final_widths=None,
        final_distances=None,
        relative=True,
        layer=0,
        datatype=0,
    ):
        """
        Add a smooth interpolating curve through the given points.

        Uses the Hobby algorithm [1]_ to calculate a smooth
        interpolating curve made of cubic Bezier segments between each
        pair of points.

        Parameters
        ----------
        points : array-like[N][2]
            Vertices in the interpolating curve.
        angles : array-like[N + 1] or None
            Tangent angles at each point (in *radians*).  Any angles
            defined as None are automatically calculated.
        curl_start : number
            Ratio between the mock curvatures at the first point and at
            its neighbor.  A value of 1 renders the first segment a good
            approximation for a circular arc.  A value of 0 will better
            approximate a straight segment.  It has no effect for closed
            curves or when an angle is defined for the first point.
        curl_end : number
            Ratio between the mock curvatures at the last point and at
            its neighbor.  It has no effect for closed curves or when an
            angle is defined for the first point.
        t_in : number or array-like[N + 1]
            Tension parameter when arriving at each point.  One value
            per point or a single value used for all points.
        t_out : number or array-like[N + 1]
            Tension parameter when leaving each point.  One value per
            point or a single value used for all points.
        cycle : bool
            If True, calculates control points for a closed curve,
            with an additional segment connecting the first and last
            points.
        tolerance : number
            Acceptable tolerance for the approximation of the curve
            function by a finite number of evaluations.
        number_of_evaluations : integer
            Initial number of points where the curve function will be
            evaluated.  According to `tolerance`, more evaluations will
            be performed.
        max_points : integer
            Elements will be fractured until each polygon has at most
            `max_points`.  If `max_points` is zero no fracture will
            occur.
        final_widths : array-like[M]
            Each element corresponds to the final width of a segment in
            the whole curve.  If an element is a number, the paths of
            this segment will have their widths linearly changed to this
            value.  If a function, it must be a function of one argument
            (that varies from 0 to 1) and returns the width of the path.
            The length of the array must be equal to the number of
            segments in the curve, i.e., M = N - 1 for an open curve and
            M = N for a closed one.
        final_distances : array-like[M]
            Each element corresponds to the final distance between paths
            of a segment in the whole curve.  If an element is a number,
            the distance between paths is linearly change to this value.
            If a function, it must be a function of one argument (that
            varies from 0 to 1) and returns the width of the path.  The
            length of the array must be equal to the number of segments
            in the curve, i.e., M = N - 1 for an open curve and M = N
            for a closed one.
        relative : bool
            If True, all coordinates in the `points` array are used as
            offsets from the current path position, i.e., if the path is
            at (1, -2) and the last point in the array is (10, 25), the
            constructed curve will end at (1 + 10, -2 + 25) = (11, 23).
            Otherwise, the points are used as absolute coordinates.
        layer : integer, list
            The GDSII layer numbers for the elements of each path.  If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255).  If the number of datatypes in the list is less
            than the number of paths, the list is repeated.

        Returns
        -------
        out : `Path`
            This object.

        Notes
        -----
        The original GDSII specification supports only a maximum of 199
        vertices per polygon.

        References
        ----------
        .. [1] Hobby, J.D.  *Discrete Comput. Geom.* (1986) 1: 123.
           `DOI: 10.1007/BF02187690
           <https://doi.org/10.1007/BF02187690>`_
        """
        if relative:
            points = numpy.vstack(([(0.0, 0.0)], points)) + numpy.array(
                (self.x, self.y)
            )
        else:
            points = numpy.vstack(([(self.x, self.y)], points))
        cta, ctb = _hobby(points, angles, curl_start, curl_end, t_in, t_out, cycle)
        if final_widths is None:
            final_widths = [None] * cta.shape[0]
        if final_distances is None:
            final_distances = [None] * cta.shape[0]
        for i in range(points.shape[0] - 1):
            self.bezier(
                [cta[i], ctb[i], points[i + 1]],
                tolerance,
                number_of_evaluations,
                max_points,
                final_widths[i],
                final_distances[i],
                False,
                layer,
                datatype,
            )
        if cycle:
            self.bezier(
                [cta[-1], ctb[-1], points[0]],
                tolerance,
                number_of_evaluations,
                max_points,
                final_widths[-1],
                final_distances[-1],
                False,
                layer,
                datatype,
            )
        return self


_pmone = numpy.array((1.0, -1.0))


class L1Path(PolygonSet):
    """
    Series of geometric objects that form a path or a collection of
    parallel paths with Manhattan geometry.

    .. deprecated:: 1.4
       `L1Path` is deprecated in favor of FlexPath and will be removed
       in a future version of Gdspy.

    Parameters
    ----------
    initial_point : array-like[2]
        Starting position of the path.
    direction : '+x', '+y', '-x', '-y'
        Starting direction of the path.
    width : number
        The initial width of each path.
    length : array-like
        Lengths of each section to add.
    turn : array-like
        Direction to turn before each section.  The sign indicate the
        turn direction (ccw is positive), and the modulus is a
        multiplicative factor for the path width after each turn.  Must
        have 1 element less then `length`.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number
        Distance between the centers of adjacent paths.
    max_points : integer
        The paths will be fractured in polygons with at most
        `max_points` (must be at least 6).  If `max_points` is zero no
        fracture will occur.
    layer : integer, list
        The GDSII layer numbers for the elements of each path.  If the
        number of layers in the list is less than the number of paths,
        the list is repeated.
    datatype : integer, list
        The GDSII datatype for the elements of each path (between 0 and
        255).  If the number of datatypes in the list is less than the
        number of paths, the list is repeated.

    Attributes
    ----------
    x : number
        Final position of the path in the x direction.
    y : number
        Final position of the path in the y direction.
    direction : '+x', '-x', '+y', '-y' or number
        Direction or angle (in *radians*) the path points to.  The
        numerical angle is returned only after a rotation of the object.
    properties : {integer: string} dictionary
        Properties for this path.

    Examples
    --------
    >>> length = [10, 30, 15, 15, 15, 15, 10]
    >>> turn = [1, -1, -1, 3, -1, 1]
    >>> l1path = gdspy.L1Path((0, 0), '+x', 2, length, turn)
    >>> myCell.add(l1path)
    """

    __slots__ = "layers", "datatypes", "polygons", "direction", "x", "y", "properties"

    def __init__(
        self,
        initial_point,
        direction,
        width,
        length,
        turn,
        number_of_paths=1,
        distance=0,
        max_points=199,
        layer=0,
        datatype=0,
    ):
        warnings.warn(
            "[GDSPY] L1Path is deprecated favor of FlexPath and will be "
            "removed in a future version of Gdspy.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        if not isinstance(layer, list):
            layer = [layer]
        if not isinstance(datatype, list):
            datatype = [datatype]
        layer = (layer * (number_of_paths // len(layer) + 1))[:number_of_paths]
        datatype = (datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths]
        w = width * 0.5
        points = len(turn) + 1 if max_points == 0 else max_points // 2 - 1
        paths = [[[], []] for ii in range(number_of_paths)]
        self.polygons = []
        self.layers = []
        self.datatypes = []
        self.properties = {}
        self.x = initial_point[0]
        self.y = initial_point[1]
        if direction == "+x":
            direction = 0
            for ii in range(number_of_paths):
                d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                paths[ii][0].append((initial_point[0], d0 + initial_point[1] - w))
                paths[ii][1].append((initial_point[0], d0 + initial_point[1] + w))
        elif direction == "+y":
            direction = 1
            for ii in range(number_of_paths):
                d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                paths[ii][0].append((d0 + initial_point[0] + w, initial_point[1]))
                paths[ii][1].append((d0 + initial_point[0] - w, initial_point[1]))
        elif direction == "-x":
            direction = 2
            for ii in range(number_of_paths):
                d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                paths[ii][0].append((initial_point[0], d0 + initial_point[1] + w))
                paths[ii][1].append((initial_point[0], d0 + initial_point[1] - w))
        elif direction == "-y":
            direction = 3
            for ii in range(number_of_paths):
                d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                paths[ii][0].append((d0 + initial_point[0] - w, initial_point[1]))
                paths[ii][1].append((d0 + initial_point[0] + w, initial_point[1]))
        for jj in range(len(turn)):
            points -= 1
            if direction == 0:
                for ii in range(number_of_paths):
                    d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                    paths[ii][0].append(
                        (self.x + length[jj] - (d0 - w) * turn[jj], paths[ii][0][-1][1])
                    )
                    paths[ii][1].append(
                        (self.x + length[jj] - (d0 + w) * turn[jj], paths[ii][1][-1][1])
                    )
                self.x += length[jj]
            elif direction == 1:
                for ii in range(number_of_paths):
                    d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                    paths[ii][0].append(
                        (paths[ii][0][-1][0], self.y + length[jj] - (d0 - w) * turn[jj])
                    )
                    paths[ii][1].append(
                        (paths[ii][1][-1][0], self.y + length[jj] - (d0 + w) * turn[jj])
                    )
                self.y += length[jj]
            elif direction == 2:
                for ii in range(number_of_paths):
                    d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                    paths[ii][0].append(
                        (self.x - length[jj] - (d0 + w) * turn[jj], paths[ii][0][-1][1])
                    )
                    paths[ii][1].append(
                        (self.x - length[jj] - (d0 - w) * turn[jj], paths[ii][1][-1][1])
                    )
                self.x -= length[jj]
            elif direction == 3:
                for ii in range(number_of_paths):
                    d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                    paths[ii][0].append(
                        (paths[ii][0][-1][0], self.y - length[jj] - (d0 + w) * turn[jj])
                    )
                    paths[ii][1].append(
                        (paths[ii][1][-1][0], self.y - length[jj] - (d0 - w) * turn[jj])
                    )
                self.y -= length[jj]
            if points == 0:
                for p in paths:
                    if direction % 2 == 0:
                        min_dist = 1e300
                        for x1 in [p[0][-2][0], p[1][-2][0]]:
                            for x2 in [p[0][-1][0], p[1][-1][0]]:
                                if abs(x1 - x2) < min_dist:
                                    x0 = 0.5 * (x1 + x2)
                                    min_dist = abs(x1 - x2)
                        p0 = (x0, p[0][-1][1])
                        p1 = (x0, p[1][-1][1])
                    else:
                        min_dist = 1e300
                        for y1 in [p[0][-2][1], p[1][-2][1]]:
                            for y2 in [p[0][-1][1], p[1][-1][1]]:
                                if abs(y1 - y2) < min_dist:
                                    y0 = 0.5 * (y1 + y2)
                                    min_dist = abs(y1 - y2)
                        p0 = (p[0][-1][0], y0)
                        p1 = (p[1][-1][0], y0)
                    self.polygons.append(
                        numpy.array(p[0][:-1] + [p0, p1] + p[1][-2::-1])
                    )
                    p[0] = [p0, p[0][-1]]
                    p[1] = [p1, p[1][-1]]
                self.layers.extend(layer)
                self.datatypes.extend(datatype)
                points = max_points // 2 - 2
            if turn[jj] > 0:
                direction = (direction + 1) % 4
            else:
                direction = (direction - 1) % 4
        if direction == 0:
            for ii in range(number_of_paths):
                d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                paths[ii][0].append((self.x + length[-1], paths[ii][0][-1][1]))
                paths[ii][1].append((self.x + length[-1], paths[ii][1][-1][1]))
            self.x += length[-1]
        elif direction == 1:
            for ii in range(number_of_paths):
                d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                paths[ii][0].append((paths[ii][0][-1][0], self.y + length[-1]))
                paths[ii][1].append((paths[ii][1][-1][0], self.y + length[-1]))
            self.y += length[-1]
        elif direction == 2:
            for ii in range(number_of_paths):
                d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                paths[ii][0].append((self.x - length[-1], paths[ii][0][-1][1]))
                paths[ii][1].append((self.x - length[-1], paths[ii][1][-1][1]))
            self.x -= length[-1]
        elif direction == 3:
            for ii in range(number_of_paths):
                d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                paths[ii][0].append((paths[ii][0][-1][0], self.y - length[-1]))
                paths[ii][1].append((paths[ii][1][-1][0], self.y - length[-1]))
            self.y -= length[jj]
        self.direction = ["+x", "+y", "-x", "-y"][direction]
        self.polygons.extend(numpy.array(p[0] + p[1][::-1]) for p in paths)
        self.layers.extend(layer)
        self.datatypes.extend(datatype)

    def __str__(self):
        return "L1Path (end at ({}, {}) towards {}, {} polygons, {} vertices, layers {}, datatypes {})".format(
            self.x,
            self.y,
            self.direction,
            len(self.polygons),
            sum([len(p) for p in self.polygons]),
            list(set(self.layers)),
            list(set(self.datatypes)),
        )

    def rotate(self, angle, center=(0, 0)):
        """
        Rotate this object.
        Parameters
        ----------
        angle : number
            The angle of rotation (in *radians*).
        center : array-like[2]
            Center point for the rotation.
        Returns
        -------
        out : `L1Path`
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle) * _mpone
        c0 = numpy.array(center)
        if isinstance(self.direction, basestring):
            self.direction = _directions_dict[self.direction] * numpy.pi
        self.direction += angle
        cur = numpy.array((self.x, self.y)) - c0
        self.x, self.y = cur * ca + cur[::-1] * sa + c0
        self.polygons = [
            (points - c0) * ca + (points - c0)[:, ::-1] * sa + c0
            for points in self.polygons
        ]
        return self


class PolyPath(PolygonSet):
    """
    Series of geometric objects that form a polygonal path or a
    collection of parallel polygonal paths.

    .. deprecated:: 1.4
       `PolyPath` is deprecated in favor of FlexPath and will be removed
       in a future version of Gdspy.

    Parameters
    ----------
    points : array-like[N][2]
        Points along the center of the path.
    width : number or array-like[N]
        Width of the path.  If an array is given, width at each
        endpoint.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number or array-like[N]
        Distance between the centers of adjacent paths.  If an array is
        given, distance at each endpoint.
    corners : 'miter' or 'bevel'
        Type of joins.
    ends : 'flush', 'round', 'extended'
        Type of end caps for the paths.
    max_points : integer
        The paths will be fractured in polygons with at most
        `max_points` (must be at least 4).  If `max_points` is zero no
        fracture will occur.
    layer : integer, list
        The GDSII layer numbers for the elements of each path.  If the
        number of layers in the list is less than the number of paths,
        the list is repeated.
    datatype : integer, list
        The GDSII datatype for the elements of each path (between 0 and
        255).  If the number of datatypes in the list is less than the
        number of paths, the list is repeated.
    Notes
    -----
    The bevel join will give strange results if the number of paths is
    greater than 1.
    """

    __slots__ = "layers", "datatypes", "polygons", "properties"

    def __init__(
        self,
        points,
        width,
        number_of_paths=1,
        distance=0,
        corners="miter",
        ends="flush",
        max_points=199,
        layer=0,
        datatype=0,
    ):
        warnings.warn(
            "[GDSPY] PolyPath is deprecated favor of FlexPath and will "
            "be removed in a future version of Gdspy.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        if not isinstance(layer, list):
            layer = [layer]
        if not isinstance(datatype, list):
            datatype = [datatype]
        if hasattr(width, "__iter__"):
            width = numpy.array(width) * 0.5
        else:
            width = numpy.array([width * 0.5])
        len_w = len(width)
        if hasattr(distance, "__iter__"):
            distance = numpy.array(distance)
        else:
            distance = numpy.array([distance])
        len_d = len(distance)
        points = numpy.array(points, dtype=float)
        self.polygons = []
        self.layers = []
        self.datatypes = []
        self.properties = {}
        if points.shape[0] == 2 and number_of_paths == 1:
            v = points[1] - points[0]
            v = v / (v[0] ** 2 + v[1] ** 2) ** 0.5
            w0 = width[0]
            w1 = width[1 % len_w]
            if ends == "round":
                a = numpy.arctan2(v[1], v[0]) + _halfpi
                self.polygons.append(
                    Round(
                        points[0],
                        w0,
                        initial_angle=a,
                        final_angle=a + numpy.pi,
                        number_of_points=33,
                    ).polygons[0]
                )
                self.polygons.append(
                    Round(
                        points[1],
                        w1,
                        initial_angle=a - numpy.pi,
                        final_angle=a,
                        number_of_points=33,
                    ).polygons[0]
                )
                self.layers.extend(layer[:1] * 2)
                self.datatypes.extend(datatype[:1] * 2)
            if ends == "extended":
                points[0] = points[0] - v * w0
                points[1] = points[1] + v * w1
            u = numpy.array((-v[1], v[0]))
            if w0 == 0:
                self.polygons.append(
                    numpy.array((points[0], points[1] - u * w1, points[1] + u * w1))
                )
            elif w1 == 0:
                self.polygons.append(
                    numpy.array((points[0] + u * w0, points[0] - u * w0, points[1]))
                )
            else:
                self.polygons.append(
                    numpy.array(
                        (
                            points[0] + u * w0,
                            points[0] - u * w0,
                            points[1] - u * w1,
                            points[1] + u * w1,
                        )
                    )
                )
            self.layers.append(layer[0])
            self.datatypes.append(datatype[0])
            return
        if corners not in ["miter", "bevel"]:
            if corners in [0, 1]:
                corners = ["miter", "bevel"][corners]
                warnings.warn(
                    "[GDSPY] Argument corners must be one of 'miter' or 'bevel'.",
                    category=DeprecationWarning,
                    stacklevel=2,
                )
            else:
                raise ValueError(
                    "[GDSPY] Argument corners must be one of 'miter' or 'bevel'."
                )
        bevel = corners == "bevel"
        if ends not in ["flush", "round", "extended"]:
            if ends in [0, 1, 2]:
                ends = ["flush", "round", "extended"][ends]
                warnings.warn(
                    "[GDSPY] Argument ends must be one of 'flush', "
                    "'round', or 'extended'.",
                    category=DeprecationWarning,
                    stacklevel=2,
                )
            else:
                raise ValueError(
                    "[GDSPY] Argument ends must be one of 'flush', "
                    "'round', or 'extended'."
                )
        if ends == "extended":
            v = points[0] - points[1]
            v = v / (v[0] ** 2 + v[1] ** 2) ** 0.5
            points[0] = points[0] + v * width[0]
            v = points[-1] - points[-2]
            v = v / (v[0] ** 2 + v[1] ** 2) ** 0.5
            points[-1] = points[-1] + v * width[(points.shape[0] - 1) % len_w]
        elif ends == "round":
            v0 = points[1] - points[0]
            angle0 = numpy.arctan2(v0[1], v0[0]) + _halfpi
            v0 = numpy.array((-v0[1], v0[0])) / (v0[0] ** 2 + v0[1] ** 2) ** 0.5
            d0 = 0.5 * (number_of_paths - 1) * distance[0]
            v1 = points[-1] - points[-2]
            angle1 = numpy.arctan2(v1[1], v1[0]) - _halfpi
            v1 = numpy.array((-v1[1], v1[0])) / (v1[0] ** 2 + v1[1] ** 2) ** 0.5
            j1w = (points.shape[0] - 1) % len_w
            j1d = (points.shape[0] - 1) % len_d
            d1 = 0.5 * (number_of_paths - 1) * distance[j1d]
            self.polygons.extend(
                (
                    Round(
                        points[0] + v0 * (ii * distance[0] - d0),
                        width[0],
                        initial_angle=angle0,
                        final_angle=angle0 + numpy.pi,
                        number_of_points=33,
                    ).polygons[0]
                    for ii in range(number_of_paths)
                )
            )
            self.polygons.extend(
                (
                    Round(
                        points[-1] + v1 * (ii * distance[j1d] - d1),
                        width[j1w],
                        initial_angle=angle1,
                        final_angle=angle1 + numpy.pi,
                        number_of_points=33,
                    ).polygons[0]
                )
                for ii in range(number_of_paths)
            )
            self.layers.extend(
                ((layer * (number_of_paths // len(layer) + 1))[:number_of_paths]) * 2
            )
            self.datatypes.extend(
                ((datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths])
                * 2
            )
        v = points[1] - points[0]
        v = numpy.array((-v[1], v[0])) / (v[0] ** 2 + v[1] ** 2) ** 0.5
        d0 = 0.5 * (number_of_paths - 1) * distance[0]
        d1 = 0.5 * (number_of_paths - 1) * distance[1 % len_d]
        paths = [
            [
                [points[0] + (ii * distance[0] - d0 - width[0]) * v],
                [points[0] + (ii * distance[0] - d0 + width[0]) * v],
            ]
            for ii in range(number_of_paths)
        ]
        p1 = [
            (
                points[1] + (ii * distance[1 % len_d] - d1 - width[1 % len_w]) * v,
                points[1] + (ii * distance[1 % len_d] - d1 + width[1 % len_w]) * v,
            )
            for ii in range(number_of_paths)
        ]
        for jj in range(1, points.shape[0] - 1):
            j0d = jj % len_d
            j0w = jj % len_w
            j1d = (jj + 1) % len_d
            j1w = (jj + 1) % len_w
            v = points[jj + 1] - points[jj]
            v = numpy.array((-v[1], v[0])) / (v[0] ** 2 + v[1] ** 2) ** 0.5
            d0 = d1
            d1 = 0.5 * (number_of_paths - 1) * distance[j1d]
            p0 = p1
            p1 = []
            pp = []
            for ii in range(number_of_paths):
                pp.append(
                    (
                        points[jj] + (ii * distance[j0d] - d0 - width[j0w]) * v,
                        points[jj] + (ii * distance[j0d] - d0 + width[j0w]) * v,
                    )
                )
                p1.append(
                    (
                        points[jj + 1] + (ii * distance[j1d] - d1 - width[j1w]) * v,
                        points[jj + 1] + (ii * distance[j1d] - d1 + width[j1w]) * v,
                    )
                )
                for kk in (0, 1):
                    p0m = paths[ii][kk][-1] - p0[ii][kk]
                    p1p = pp[ii][kk] - p1[ii][kk]
                    vec = p0m[0] * p1p[1] - p1p[0] * p0m[1]
                    if abs(vec) > 1e-30:
                        p = (
                            _pmone
                            * (
                                p0m * p1p[::-1] * p1[ii][kk]
                                - p1p * p0m[::-1] * p0[ii][kk]
                                + p0m * p1p * (p0[ii][kk][::-1] - p1[ii][kk][::-1])
                            )
                            / vec
                        )
                        l0 = (p - pp[ii][kk]) * p1p
                        l1 = (p - p0[ii][kk]) * p0m
                        if bevel and l0[0] + l0[1] > 0 and l1[0] + l1[1] < 0:
                            paths[ii][kk].append(p0[ii][kk])
                            paths[ii][kk].append(pp[ii][kk])
                        else:
                            paths[ii][kk].append(p)
                if (
                    max_points > 0
                    and len(paths[ii][0]) + len(paths[ii][1]) + 3 > max_points
                ):
                    diff = paths[ii][0][0] - paths[ii][1][0]
                    if diff[0] ** 2 + diff[1] ** 2 == 0:
                        paths[ii][1] = paths[ii][1][1:]
                    diff = paths[ii][0][-1] - paths[ii][1][-1]
                    if diff[0] ** 2 + diff[1] ** 2 == 0:
                        self.polygons.append(
                            numpy.array(paths[ii][0] + paths[ii][1][-2::-1])
                        )
                    else:
                        self.polygons.append(
                            numpy.array(paths[ii][0] + paths[ii][1][::-1])
                        )
                    paths[ii][0] = paths[ii][0][-1:]
                    paths[ii][1] = paths[ii][1][-1:]
                    self.layers.append(layer[ii % len(layer)])
                    self.datatypes.append(datatype[ii % len(datatype)])
        for ii in range(number_of_paths):
            diff = paths[ii][0][0] - paths[ii][1][0]
            if diff[0] ** 2 + diff[1] ** 2 == 0:
                paths[ii][1] = paths[ii][1][1:]
            diff = p1[ii][0] - p1[ii][1]
            if diff[0] ** 2 + diff[1] ** 2 != 0:
                paths[ii][0].append(p1[ii][0])
            paths[ii][1].append(p1[ii][1])
        self.polygons.extend(numpy.array(pol[0] + pol[1][::-1]) for pol in paths)
        self.layers.extend(
            (layer * (number_of_paths // len(layer) + 1))[:number_of_paths]
        )
        self.datatypes.extend(
            (datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths]
        )

    def __str__(self):
        return "PolyPath ({} polygons, {} vertices, layers {}, datatypes {})".format(
            len(self.polygons),
            sum([len(p) for p in self.polygons]),
            list(set(self.layers)),
            list(set(self.datatypes)),
        )


from gdspy.path import _func_bezier
