########################################################################
##                                                                    ##
##  Copyright 2009-2016 Lucas Heitzmann Gabrielli                     ##
##                                                                    ##
##  This file is part of gdspy.                                       ##
##                                                                    ##
##  gdspy is free software: you can redistribute it and/or modify it  ##
##  under the terms of the GNU General Public License as published    ##
##  by the Free Software Foundation, either version 3 of the          ##
##  License, or any later version.                                    ##
##                                                                    ##
##  gdspy is distributed in the hope that it will be useful, but      ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of        ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     ##
##  GNU General Public License for more details.                      ##
##                                                                    ##
##  You should have received a copy of the GNU General Public         ##
##  License along with gdspy.  If not, see                            ##
##  <http://www.gnu.org/licenses/>.                                   ##
##                                                                    ##
########################################################################

from __future__ import absolute_import

import struct
import datetime
import warnings
import numpy

from gdspy import boolext
from gdspy import clipper
from gdspy.viewer import LayoutViewer

__version__ = '0.8.1'
__doc__ = """
gdspy is a Python module that allows the creation of GDSII stream files.

Many features of the GDSII format are implemented, such as cell
references and arrays, but the support for fonts is quite limited.  Text
is only available through polygonal objects.

If the Python Imaging Library is installed, it can be used to output the
geometry created to an image file.
"""

_halfpi = 0.5 * numpy.pi

def _eight_byte_real(value):
    """
    Convert a number into the GDSII 8 byte real format.

    Parameters
    ----------
    value : number
        The number to be converted.

    Returns
    -------
    out : string
        The GDSII binary string that represents ``value``.
    """
    byte1 = 0
    byte2 = 0
    short3 = 0
    long4 = 0
    if value != 0:
        if value < 0:
            byte1 = 0x80
            value = -value
        exponent = int(numpy.floor(numpy.log2(value) * 0.25))
        mantissa = long(value * 16L**(14 - exponent))
        while mantissa >= 72057594037927936L:
            exponent += 1
            mantissa = long(value * 16L**(14 - exponent))
        byte1 += exponent + 64
        byte2 = (mantissa // 281474976710656L)
        short3 = (mantissa % 281474976710656L) // 4294967296L
        long4 = mantissa % 4294967296L
    return struct.pack(">HHL", byte1 * 256 + byte2, short3, long4)


def _eight_byte_real_to_float(value):
    """
    Convert a number from GDSII 8 byte real format to float.

    Parameters
    ----------
    value : string
        The GDSII binary string representation of the number.

    Returns
    -------
    out : float
        The number represented by ``value``.
    """
    short1, short2, long3 = struct.unpack('>HHL', value)
    exponent = (short1 & 0x7f00) // 256
    mantissa = (((short1 & 0x00ff) * 65536L + short2) * 4294967296L
                + long3) / 72057594037927936.0
    return (-1 if (short1 & 0x8000) else 1) * mantissa * 16L ** (exponent - 64)


class Polygon:
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
    verbose : bool
        If False, warnings about the number of vertices of the polygon will
        be suppressed.

    Notes
    -----
    The last point should not be equal to the first (polygons are
    automatically closed).

    The GDSII specification supports only a maximum of 199 vertices per
    polygon.

    Examples
    --------
    >>> triangle_pts = [(0, 40), (15, 40), (10, 50)]
    >>> triangle = gdspy.Polygon(triangle_pts)
    >>> myCell.add(triangle)
    """

    def __init__(self, points, layer=0, datatype=0, verbose=True) :
        if len(points) > 199 and verbose:
            warnings.warn("[GDSPY] A polygon with more than 199 points was created (not officially supported by the GDSII format).", stacklevel=2)
        self.layer = layer
        self.points = numpy.array(points)
        self.datatype = datatype

    def __str__(self):
        return "Polygon ({} vertices, layer {}, datatype {})".format(
            len(self.points), self.layer, self.datatype)

    def to_gds(self, multiplier):
        """
        Convert this object to a GDSII element.

        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            element.

        Returns
        -------
        out : string
            The GDSII binary string that represents this object.
        """
        if len(self.points) > 4094:
            raise ValueError("[GDSPY] Polygons with more than 4094 are not supported by the GDSII format.")
        return struct.pack('>10h', 4, 0x0800, 6, 0x0D02, self.layer, 6, 0x0E02, self.datatype, 12 + 8 * len(self.points), 0x1003) + b''.join(struct.pack('>2l', int(round(point[0] * multiplier)), int(round(point[1] * multiplier))) for point in self.points) + struct.pack('>2l2h', int(round(self.points[0][0] * multiplier)), int(round(self.points[0][1] * multiplier)), 4, 0x1100)

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
        out : ``Polygon``
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle)
        sa = numpy.array((-sa, sa))
        c0 = numpy.array(center)
        self.points = (self.points - c0) * ca + (self.points - c0)[:,::-1] * sa + c0
        return self

    def area(self, by_spec=False):
        """
        Calculate the total area of this object.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary ``{(layer, datatype):
            area}``.

        Returns
        -------
        out : number, dictionary
            Area of this object.
        """
        poly_area = 0
        for ii in range(1, self.points.shape[0] - 1):
            poly_area += (self.points[0][0] - self.points[ii + 1][0]) * (self.points[ii][1] - self.points[0][1]) - (self.points[0][1] - self.points[ii + 1][1]) * (self.points[ii][0] - self.points[0][0])
        if by_spec:
            return {(self.layer, self.datatype): 0.5 * abs(poly_area)}
        else:
            return 0.5 * abs(poly_area)

    def fracture(self, max_points=199):
        """
        Slice this polygon in the horizontal and vertical directions so that
        each resulting piece has at most ``max_points``.

        Parameters
        ----------
        max_points : integer
            Maximal number of points in each resulting polygon (must be greater
            than 4).

        Returns
        -------
        out : ``PolygonSet``
            Resulting polygons from the fracture operation.
        """
        out_polygons = [self.points]
        if max_points > 4:
            ii = 0
            while ii < len(out_polygons):
                if len(out_polygons[ii]) > max_points:
                    pts0 = [x[0] for x in out_polygons[ii]]
                    pts1 = [x[1] for x in out_polygons[ii]]
                    pts0.sort()
                    pts1.sort()
                    if pts0[-1] - pts0[0] > pts1[-1] - pts1[0]:
                        ## Vertical cuts
                        chopped = _chop(out_polygons[ii], (pts0[len(pts0) // 2] + pts0[len(pts0) // 2 + 1]) / 2, 0)
                    else:
                        ## Horizontal cuts
                        chopped = _chop(out_polygons[ii], (pts1[len(pts1) // 2] + pts1[len(pts1) // 2 + 1]) / 2, 1)
                    out_polygons.pop(ii)
                    out_polygons += chopped[0]
                    out_polygons += chopped[1]
                else:
                    ii += 1
        return PolygonSet(out_polygons, self.layer, self.datatype)

    def fillet(self, radius, points_per_2pi=128, max_points=199):
        """
        Round the corners of this polygon and fractures it into polygons with
        less vertices if necessary.

        Parameters
        ----------
        radius : number
            Radius of the corners.
        points_per_2pi : integer
            Number of vertices used to approximate a full circle.  The number of
            vertices in each corner of the polygon will be the fraction of this
            number corresponding to the angle encompassed by that corner with
            respect to 2 pi.
        max_points : integer
            Maximal number of points in each resulting polygon (must be greater
            than 4).

        Returns
        -------
        out : ``Polygon`` or ``PolygonSet``
            If no fracturing occurs, return this object; otherwise return a
            ``PolygonSet`` with the fractured result (this object will have more
            than ``max_points`` vertices).
        """
        two_pi = 2 * numpy.pi
        vec = self.points.astype(float) - numpy.roll(self.points, 1, 0)
        length = numpy.sqrt(numpy.sum(vec ** 2, 1))
        ii = numpy.flatnonzero(length)
        if len(ii) < len(length):
            self.points = self.points[ii]
            vec = self.points - numpy.roll(self.points, 1, 0)
            length = numpy.sqrt(numpy.sum(vec ** 2, 1))
        vec[:,0] /= length
        vec[:,1] /= length
        dvec = numpy.roll(vec, -1, 0) - vec
        norm = numpy.sqrt(numpy.sum(dvec ** 2, 1))
        ii = numpy.flatnonzero(norm)
        dvec[ii,0] /= norm[ii]
        dvec[ii,1] /= norm[ii]
        theta = numpy.arccos(numpy.sum(numpy.roll(vec, -1, 0) * vec, 1))
        ct = numpy.cos(theta * 0.5)
        tt = numpy.tan(theta * 0.5)

        new_points = []
        for ii in range(-1, len(self.points) - 1):
            if theta[ii] > 0:
                a0 = -vec[ii] * tt[ii] - dvec[ii] / ct[ii]
                a0 = numpy.arctan2(a0[1], a0[0])
                a1 = vec[ii + 1] * tt[ii] - dvec[ii] / ct[ii]
                a1 = numpy.arctan2(a1[1], a1[0])
                if a1 - a0 > numpy.pi:
                    a1 -= two_pi
                elif a1 - a0 < -numpy.pi:
                    a1 += two_pi
                n = max(int(numpy.ceil(abs(a1 - a0) / two_pi * points_per_2pi) + 0.5), 2)
                a = numpy.linspace(a0, a1, n)
                l = radius * tt[ii]
                if l > 0.49 * length[ii]:
                    r = 0.49 * length[ii] / tt[ii]
                    l = 0.49 * length[ii]
                else:
                    r = radius
                if l > 0.49 * length[ii + 1]:
                    r = 0.49 * length[ii + 1] / tt[ii]
                new_points += list(r * dvec[ii] / ct[ii] + self.points[ii] + numpy.vstack((r * numpy.cos(a), r * numpy.sin(a))).transpose())
            else:
                new_points.append(self.points[ii])

        self.points = numpy.array(new_points)
        if len(self.points) > max_points:
            return self.fracture()
        else:
            return self


class PolygonSet:
    """
    Set of polygonal objects.

    Parameters
    ----------
    polygons : list of array-like[N][2]
        List containing the coordinates of the vertices of each polygon.
    layer : integer
        The GDSII layer number for this element.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).
    verbose : bool
        If False, warnings about the number of vertices of the polygons will
        be suppressed.

    Notes
    -----
    The last point should not be equal to the first (polygons are
    automatically closed).

    The GDSII specification supports only a maximum of 199 vertices per
    polygon.
    """

    def __init__(self, polygons, layer=0, datatype=0, verbose=True):
        self.layers = [layer] * len(polygons)
        self.datatypes = [datatype] * len(polygons)
        self.polygons = [None] * len(polygons)
        for i in range(len(polygons)):
            self.polygons[i] = numpy.array(polygons[i])
            if len(polygons[i]) > 199 and verbose:
                warnings.warn("[GDSPY] A polygon with more than 199 points was created (not officially supported by the GDSII format).", stacklevel=2)

    def __str__(self):
        return "PolygonSet ({} polygons, {} vertices, layers {}, datatypes {})".format(len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))

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
        out : ``PolygonSet``
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle)
        sa = numpy.array((-sa, sa))
        c0 = numpy.array(center)
        self.polygons = [(points - c0) * ca + (points - c0)[:,::-1] * sa + c0 for points in self.polygons]
        return self

    def to_gds(self, multiplier):
        """
        Convert this object to a series of GDSII elements.

        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            elements.

        Returns
        -------
        out : string
            The GDSII binary string that represents this object.
        """
        data = []
        for ii in range(len(self.polygons)):
            if len(self.polygons[ii]) > 4094:
                raise ValueError("[GDSPY] Polygons with more than 4094 are not supported by the GDSII format.")
            data.append(struct.pack('>10h', 4, 0x0800, 6, 0x0D02, self.layers[ii], 6, 0x0E02, self.datatypes[ii], 12 + 8 * len(self.polygons[ii]), 0x1003))
            data.extend(struct.pack('>2l', int(round(point[0] * multiplier)), int(round(point[1] * multiplier))) for point in self.polygons[ii])
            data.append(struct.pack('>2l2h', int(round(self.polygons[ii][0][0] * multiplier)), int(round(self.polygons[ii][0][1] * multiplier)), 4, 0x1100))
        return b''.join(data)

    def area(self, by_spec=False):
        """
        Calculate the total area of the path(s).

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary ``{(layer, datatype):
            area}``.

        Returns
        -------
        out : number, dictionary
            Area of this object.
        """
        if by_spec:
            path_area = {}
            for jj in range(len(self.polygons)):
                poly_area = 0
                for ii in range(1, len(self.polygons[jj]) - 1):
                    poly_area += (self.polygons[jj][0][0] - self.polygons[jj][ii + 1][0]) * (self.polygons[jj][ii][1] - self.polygons[jj][0][1]) - (self.polygons[jj][0][1] - self.polygons[jj][ii + 1][1]) * (self.polygons[jj][ii][0] - self.polygons[jj][0][0])
                key = (self.layers[jj], self.datatypes[jj])
                if path_area.has_key(key):
                    path_area[key] += 0.5 * abs(poly_area)
                else:
                    path_area[key] = 0.5 * abs(poly_area)
        else:
            path_area = 0
            for points in self.polygons:
                poly_area = 0
                for ii in range(1, len(points) - 1):
                    poly_area += (points[0][0] - points[ii + 1][0]) * (points[ii][1] - points[0][1]) - (points[0][1] - points[ii + 1][1]) * (points[ii][0] - points[0][0])
                path_area += 0.5 * abs(poly_area)
        return path_area

    def fracture(self, max_points=199):
        """
        Slice these polygons in the horizontal and vertical directions so that
        each resulting piece has at most ``max_points``.  This operation occurs
        in place.

        Parameters
        ----------
        max_points : integer
            Maximal number of points in each resulting polygon (must be greater
            than 4).

        Returns
        -------
        out : ``PolygonSet``
            This object.
        """
        if max_points > 4:
            ii = 0
            while ii < len(self.polygons):
                if len(self.polygons[ii]) > max_points:
                    pts0 = [x[0] for x in self.polygons[ii]]
                    pts1 = [x[1] for x in self.polygons[ii]]
                    pts0.sort()
                    pts1.sort()
                    if pts0[-1] - pts0[0] > pts1[-1] - pts1[0]:
                        ## Vertical cuts
                        chopped = _chop(self.polygons[ii], (pts0[len(pts0) // 2] + pts0[len(pts0) // 2 + 1]) / 2, 0)
                    else:
                        ## Horizontal cuts
                        chopped = _chop(self.polygons[ii], (pts1[len(pts1) // 2] + pts1[len(pts1) // 2 + 1]) / 2, 1)
                    self.polygons.pop(ii)
                    layer = self.layers.pop(ii)
                    datatype = self.datatypes.pop(ii)
                    self.polygons += [numpy.array(x) for x in chopped[0] + chopped[1]]
                    self.layers += [layer] * (len(chopped[0]) + len(chopped[1]))
                    self.datatypes += [datatype] * (len(chopped[0]) + len(chopped[1]))
                else:
                    ii += 1
        return self

    def fillet(self, radius, points_per_2pi=128, max_points=199):
        """
        Round the corners of these polygons and fractures them into polygons
        with less vertices if necessary.

        Parameters
        ----------
        radius : number
            Radius of the corners.
        points_per_2pi : integer
            Number of vertices used to approximate a full circle.  The number of
            vertices in each corner of the polygon will be the fraction of this
            number corresponding to the angle encompassed by that corner with
            respect to 2 pi.
        max_points : integer
            Maximal number of points in each resulting polygon (must be greater
            than 4).

        Returns
        -------
        out : ``PolygonSet``
            This object.
        """
        two_pi = 2 * numpy.pi
        fracture = False

        for jj in range(len(self.polygons)):
            vec = self.polygons[jj].astype(float) - numpy.roll(self.polygons[jj], 1, 0)
            length = numpy.sqrt(numpy.sum(vec ** 2, 1))
            ii = numpy.flatnonzero(length)
            if len(ii) < len(length):
                self.polygons[jj] = self.polygons[jj][ii]
                vec = self.polygons[jj] - numpy.roll(self.polygons[jj], 1, 0)
                length = numpy.sqrt(numpy.sum(vec ** 2, 1))
            vec[:,0] /= length
            vec[:,1] /= length
            dvec = numpy.roll(vec, -1, 0) - vec
            norm = numpy.sqrt(numpy.sum(dvec ** 2, 1))
            ii = numpy.flatnonzero(norm)
            dvec[ii,0] /= norm[ii]
            dvec[ii,1] /= norm[ii]
            theta = numpy.arccos(numpy.sum(numpy.roll(vec, -1, 0) * vec, 1))
            ct = numpy.cos(theta * 0.5)
            tt = numpy.tan(theta * 0.5)

            new_points = []
            for ii in range(-1, len(self.polygons[jj]) - 1):
                if theta[ii] > 0:
                    a0 = -vec[ii] * tt[ii] - dvec[ii] / ct[ii]
                    a0 = numpy.arctan2(a0[1], a0[0])
                    a1 = vec[ii + 1] * tt[ii] - dvec[ii] / ct[ii]
                    a1 = numpy.arctan2(a1[1], a1[0])
                    if a1 - a0 > numpy.pi:
                        a1 -= two_pi
                    elif a1 - a0 < -numpy.pi:
                        a1 += two_pi
                    n = max(int(numpy.ceil(abs(a1 - a0) / two_pi * points_per_2pi) + 0.5), 2)
                    a = numpy.linspace(a0, a1, n)
                    l = radius * tt[ii]
                    if l > 0.49 * length[ii]:
                        r = 0.49 * length[ii] / tt[ii]
                        l = 0.49 * length[ii]
                    else:
                        r = radius
                    if l > 0.49 * length[ii + 1]:
                        r = 0.49 * length[ii + 1] / tt[ii]
                    new_points += list(r * dvec[ii] / ct[ii] + self.polygons[jj][ii] + numpy.vstack((r * numpy.cos(a), r * numpy.sin(a))).transpose())
                else:
                    new_points.append(self.polygons[jj][ii])
            self.polygons[jj] = numpy.array(new_points)
            if len(new_points) > max_points:
                fracture = True

        if fracture:
            self.fracture(max_points)
        return self


class Rectangle(Polygon):
    """
    Rectangular geometric object.

    Parameters
    ----------
    point1 : array-like[2]
        Coordinates of a corner of the rectangle.
    point2 : array-like[2]
        Coordinates of the corner of the rectangle opposite to ``point1``.
    layer : integer
        The GDSII layer number for this element.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).

    Examples
    --------
    >>> rectangle = gdspy.Rectangle((0, 0), (10, 20))
    >>> myCell.add(rectangle)
    """

    def __init__(self, point1, point2, layer=0, datatype=0):
        self.layer = layer
        self.points = numpy.array([[point1[0], point1[1]], [point1[0], point2[1]], [point2[0], point2[1]], [point2[0], point1[1]]])
        self.datatype = datatype

    def __str__(self):
        return "Rectangle (({0[0]}, {0[1]}) to ({1[0]}, {1[1]}), layer {2}, datatype {3})".format(self.points[0], self.points[2], self.layer, self.datatype)

    def __repr__(self):
        return "Rectangle({2}, ({0[0]}, {0[1]}), ({1[0]}, {1[1]}), {3})".format(self.points[0], self.points[2], self.layer, self.datatype)


class Round(PolygonSet):
    """
    Circular geometric object.
    Represent a circle, a circular section, a ring or a ring section.

    Parameters
    ----------
    center : array-like[2]
        Coordinates of the center of the circle/ring.
    radius : number
        Radius of the circle/outer radius of the ring.
    inner_radius : number
        Inner radius of the ring.
    initial_angle : number
        Initial angle of the circular/ring section (in *radians*).
    final_angle : number
        Final angle of the circular/ring section (in *radians*).
    number_of_points : integer or float
        If integer: number of vertices that form the object (polygonal
        approximation). If float: approximate curvature resolution. The
        actual number of points is automatically calculated.
    max_points : integer
        if ``number_of_points > max_points``, the element will be fractured
        in smaller polygons with at most ``max_points`` each.
    layer : integer
        The GDSII layer number for this element.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).

    Notes
    -----
    The GDSII specification supports only a maximum of 199 vertices per
    polygon.

    Examples
    --------
    >>> circle = gdspy.Round((30, 5), 8)
    >>> ring = gdspy.Round((50, 5), 8, inner_radius=5)
    >>> pie_slice = gdspy.Round((30, 25), 8, initial_angle=0,
    ...                             final_angle=-5.0*numpy.pi/6.0)
    >>> arc = gdspy.Round((50, 25), 8, inner_radius=5,
    ...                       initial_angle=-5.0*numpy.pi/6.0,
    ...                       final_angle=0)
    """

    def __init__(self, center, radius, inner_radius=0, initial_angle=0, final_angle=0, number_of_points=0.01, max_points=199, layer=0, datatype=0):
        if isinstance(number_of_points, float):
            if inner_radius <= 0:
                if final_angle == initial_angle:
                    number_of_points = int(2 * radius * numpy.pi / number_of_points + 0.5)
                else:
                    number_of_points = int(abs(final_angle - initial_angle) * radius / number_of_points + 0.5) + 2
            else:
                if final_angle == initial_angle:
                    number_of_points = 2 * int(2 * radius * numpy.pi / number_of_points + 0.5) + 2
                else:
                    number_of_points = 2 * int(abs(final_angle - initial_angle) * radius / number_of_points + 0.5) + 2
        number_of_points = max(number_of_points, 3)
        pieces = int(numpy.ceil(number_of_points / float(max_points)))
        number_of_points = number_of_points // pieces
        self.layers = [layer] * pieces
        self.datatypes = [datatype] * pieces
        self.polygons = [numpy.zeros((number_of_points, 2)) for _ in range(pieces)]
        if final_angle == initial_angle and pieces > 1:
            final_angle += 2 * numpy.pi
        angles = numpy.linspace(initial_angle, final_angle, pieces + 1)
        for ii in range(pieces):
            if angles[ii+1] == angles[ii]:
                if inner_radius <= 0:
                    angle = numpy.arange(number_of_points) * 2.0 * numpy.pi / number_of_points
                    self.polygons[ii][:,0] = numpy.cos(angle)
                    self.polygons[ii][:,1] = numpy.sin(angle)
                    self.polygons[ii] = self.polygons[ii] * radius + numpy.array(center)
                else:
                    n2 = number_of_points // 2
                    n1 = number_of_points - n2
                    angle = numpy.arange(n1) * 2.0 * numpy.pi / (n1 - 1.0)
                    self.polygons[ii][:n1,0] = numpy.cos(angle) * radius + center[0]
                    self.polygons[ii][:n1,1] = numpy.sin(angle) * radius + center[1]
                    angle = numpy.arange(n2) * -2.0 * numpy.pi / (n2 - 1.0)
                    self.polygons[ii][n1:,0] = numpy.cos(angle) * inner_radius + center[0]
                    self.polygons[ii][n1:,1] = numpy.sin(angle) * inner_radius + center[1]
            else:
                if inner_radius <= 0:
                    angle = numpy.linspace(angles[ii], angles[ii+1], number_of_points - 1)
                    self.polygons[ii][1:,0] = numpy.cos(angle)
                    self.polygons[ii][1:,1] = numpy.sin(angle)
                    self.polygons[ii] = self.polygons[ii] * radius + numpy.array(center)
                else:
                    n2 = number_of_points // 2
                    n1 = number_of_points - n2
                    angle = numpy.linspace(angles[ii], angles[ii+1], n1)
                    self.polygons[ii][:n1,0] = numpy.cos(angle) * radius + center[0]
                    self.polygons[ii][:n1,1] = numpy.sin(angle) * radius + center[1]
                    angle = numpy.linspace(angles[ii+1], angles[ii], n2)
                    self.polygons[ii][n1:,0] = numpy.cos(angle) * inner_radius + center[0]
                    self.polygons[ii][n1:,1] = numpy.sin(angle) * inner_radius + center[1]

    def __str__(self):
        return "Round ({} polygons, {} vertices, layers {}, datatypes {})".format(len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))


class Text(PolygonSet):
    """
    Polygonal text object.

    Each letter is formed by a series of polygons.

    Parameters
    ----------
    text : string
        The text to be converted in geometric objects.
    size : number
        Base size of each character.
    position : array-like[2]
        Text position (lower left corner).
    horizontal : bool
        If ``True``, the text is written from left to right; if ``False``,
        from top to bottom.
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
    _font = {
        '!':[[(2, 2), (3, 2), (3, 3), (2, 3)], [(2, 4), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (2, 7), (2, 6), (2, 5)]],
        '"':[[(1, 7), (2, 7), (2, 8), (2, 9), (1, 9), (1, 8)], [(3, 7), (4, 7), (4, 8), (4, 9), (3, 9), (3, 8)]],
        '#':[[(0, 3), (1, 3), (1, 2), (2, 2), (2, 3), (2, 4), (2, 5), (3, 5), (3, 4), (2, 4), (2, 3), (3, 3), (3, 2), (4, 2), (4, 3), (5, 3), (5, 4), (4, 4), (4, 5), (5, 5), (5, 6), (4, 6), (4, 7), (3, 7), (3, 6), (2, 6), (2, 7), (1, 7), (1, 6), (0, 6), (0, 5), (1, 5), (1, 4), (0, 4)]],
        '$':[[(0, 2), (1, 2), (2, 2), (2, 1), (3, 1), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (4, 4), (4, 5), (3, 5), (3, 6), (3, 7), (4, 7), (5, 7), (5, 8), (4, 8), (3, 8), (3, 9), (2, 9), (2, 8), (1, 8), (1, 7), (2, 7), (2, 6), (1, 6), (1, 5), (2, 5), (2, 4), (2, 3), (1, 3), (0, 3)], [(0, 6), (1, 6), (1, 7), (0, 7)], [(4, 3), (5, 3), (5, 4), (4, 4)]],
        '%':[[(0, 2), (1, 2), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 7), (1, 7), (2, 7), (2, 8), (2, 9), (1, 9), (0, 9), (0, 8)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (3, 4), (3, 3)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        '&':[[(0, 3), (1, 3), (1, 4), (1, 5), (0, 5), (0, 4)], [(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 2), (2, 2), (3, 2), (3, 3), (2, 3), (1, 3)], [(1, 5), (2, 5), (3, 5), (3, 6), (3, 7), (3, 8), (2, 8), (2, 7), (2, 6), (1, 6)], [(1, 8), (2, 8), (2, 9), (1, 9)], [(3, 3), (4, 3), (4, 4), (4, 5), (3, 5), (3, 4)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 5), (5, 5), (5, 6), (4, 6)]],
        '\'':[[(2, 7), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8)]],
        '(':[[(1, 4), (2, 4), (2, 5), (2, 6), (2, 7), (1, 7), (1, 6), (1, 5)], [(2, 3), (3, 3), (3, 4), (2, 4)], [(2, 7), (3, 7), (3, 8), (2, 8)], [(3, 2), (4, 2), (4, 3), (3, 3)], [(3, 8), (4, 8), (4, 9), (3, 9)]],
        ')':[[(1, 2), (2, 2), (2, 3), (1, 3)], [(1, 8), (2, 8), (2, 9), (1, 9)], [(2, 3), (3, 3), (3, 4), (2, 4)], [(2, 7), (3, 7), (3, 8), (2, 8)], [(3, 4), (4, 4), (4, 5), (4, 6), (4, 7), (3, 7), (3, 6), (3, 5)]],
        '*':[[(0, 2), (1, 2), (1, 3), (0, 3)], [(0, 4), (1, 4), (1, 3), (2, 3), (2, 2), (3, 2), (3, 3), (4, 3), (4, 4), (5, 4), (5, 5), (4, 5), (4, 6), (3, 6), (3, 7), (2, 7), (2, 6), (1, 6), (1, 5), (0, 5)], [(0, 6), (1, 6), (1, 7), (0, 7)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        '+':[[(0, 4), (1, 4), (2, 4), (2, 3), (2, 2), (3, 2), (3, 3), (3, 4), (4, 4), (5, 4), (5, 5), (4, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6), (2, 5), (1, 5), (0, 5)]],
        ',':[[(1, 0), (2, 0), (2, 1), (1, 1)], [(2, 1), (3, 1), (3, 2), (3, 3), (2, 3), (2, 2)]],
        '-':[[(0, 4), (1, 4), (2, 4), (3, 4), (4, 4), (5, 4), (5, 5), (4, 5), (3, 5), (2, 5), (1, 5), (0, 5)]],
        '.':[[(2, 2), (3, 2), (3, 3), (2, 3)]],
        '/':[[(0, 2), (1, 2), (1, 3), (1, 4), (0, 4), (0, 3)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        '0':[[(0, 3), (1, 3), (1, 4), (2, 4), (2, 5), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 6), (4, 6), (4, 5), (4, 4), (4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (3, 7)]],
        '1':[[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (1, 8), (1, 7), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (1, 3)]],
        '2':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 7), (1, 7), (1, 8), (0, 8)], [(1, 4), (2, 4), (3, 4), (3, 5), (2, 5), (1, 5)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '3':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '4':[[(0, 4), (1, 4), (2, 4), (3, 4), (3, 3), (3, 2), (4, 2), (4, 3), (4, 4), (5, 4), (5, 5), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (3, 9), (2, 9), (2, 8), (3, 8), (3, 7), (3, 6), (3, 5), (2, 5), (1, 5), (1, 6), (0, 6), (0, 5)], [(1, 6), (2, 6), (2, 7), (2, 8), (1, 8), (1, 7)]],
        '5':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 5), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)]],
        '6':[[(0, 3), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)]],
        '7':[[(0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (4, 7), (4, 6), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (2, 5), (2, 4), (2, 3)], [(3, 5), (4, 5), (4, 6), (3, 6)]],
        '8':[[(0, 3), (1, 3), (1, 4), (1, 5), (0, 5), (0, 4)], [(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '9':[[(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 4), (4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (3, 6), (2, 6), (1, 6)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)]],
        ':':[[(2, 2), (3, 2), (3, 3), (2, 3)], [(2, 5), (3, 5), (3, 6), (2, 6)]],
        ';':[[(1, 0), (2, 0), (2, 1), (1, 1)], [(2, 1), (3, 1), (3, 2), (3, 3), (2, 3), (2, 2)], [(2, 4), (3, 4), (3, 5), (2, 5)]],
        '<':[[(0, 5), (1, 5), (1, 6), (0, 6)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(1, 6), (2, 6), (2, 7), (1, 7)], [(2, 3), (3, 3), (4, 3), (4, 4), (3, 4), (2, 4)], [(2, 7), (3, 7), (4, 7), (4, 8), (3, 8), (2, 8)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 8), (5, 8), (5, 9), (4, 9)]],
        '=':[[(0, 3), (1, 3), (2, 3), (3, 3), (4, 3), (5, 3), (5, 4), (4, 4), (3, 4), (2, 4), (1, 4), (0, 4)], [(0, 5), (1, 5), (2, 5), (3, 5), (4, 5), (5, 5), (5, 6), (4, 6), (3, 6), (2, 6), (1, 6), (0, 6)]],
        '>':[[(0, 2), (1, 2), (1, 3), (0, 3)], [(0, 8), (1, 8), (1, 9), (0, 9)], [(1, 3), (2, 3), (3, 3), (3, 4), (2, 4), (1, 4)], [(1, 7), (2, 7), (3, 7), (3, 8), (2, 8), (1, 8)], [(3, 4), (4, 4), (4, 5), (3, 5)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 5), (5, 5), (5, 6), (4, 6)]],
        '?':[[(0, 7), (1, 7), (1, 8), (0, 8)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 2), (3, 2), (3, 3), (2, 3)], [(2, 4), (3, 4), (3, 5), (2, 5)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        '@':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 4), (3, 4), (4, 4), (4, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6), (2, 5)], [(4, 5), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6)]],
        'A':[[(0, 2), (1, 2), (1, 3), (1, 4), (2, 4), (3, 4), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (4, 5), (3, 5), (2, 5), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)]],
        'B':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        'C':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9)]],
        'D':[[(0, 2), (1, 2), (2, 2), (3, 2), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 7), (4, 7), (4, 8), (3, 8)], [(4, 4), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5)]],
        'E':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'F':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'G':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (3, 6), (2, 6), (2, 5), (3, 5), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9)]],
        'H':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'I':[[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (1, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (1, 3)]],
        'J':[[(0, 3), (1, 3), (1, 4), (0, 4)], [(0, 8), (1, 8), (2, 8), (3, 8), (3, 7), (3, 6), (3, 5), (3, 4), (3, 3), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(1, 2), (2, 2), (3, 2), (3, 3), (2, 3), (1, 3)]],
        'K':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (2, 6), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(2, 4), (3, 4), (3, 5), (2, 5)], [(2, 6), (3, 6), (3, 7), (2, 7)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 7), (4, 7), (4, 8), (3, 8)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 8), (5, 8), (5, 9), (4, 9)]],
        'L':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)]],
        'M':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 7), (2, 8), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(2, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6)], [(3, 7), (4, 7), (4, 6), (4, 5), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (3, 8)]],
        'N':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 7), (2, 8), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(2, 5), (3, 5), (3, 6), (3, 7), (2, 7), (2, 6)], [(3, 4), (4, 4), (4, 3), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (4, 5), (3, 5)]],
        'O':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'P':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        'Q':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4), (3, 4), (3, 3), (2, 3), (1, 3)], [(1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9)], [(2, 4), (3, 4), (3, 5), (2, 5)]],
        'R':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (3, 4), (4, 4), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6), (1, 7), (1, 8), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (4, 3)], [(4, 6), (5, 6), (5, 7), (5, 8), (4, 8), (4, 7)]],
        'S':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 5), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (1, 6)], [(1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9)], [(4, 3), (5, 3), (5, 4), (5, 5), (4, 5), (4, 4)]],
        'T':[[(0, 8), (1, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)]],
        'U':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'V':[[(0, 5), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6)], [(1, 3), (2, 3), (2, 4), (2, 5), (1, 5), (1, 4)], [(2, 2), (3, 2), (3, 3), (2, 3)], [(3, 3), (4, 3), (4, 4), (4, 5), (3, 5), (3, 4)], [(4, 5), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6)]],
        'W':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (2, 3), (1, 3)], [(2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (2, 6), (2, 5), (2, 4)], [(3, 2), (4, 2), (4, 3), (3, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'X':[[(0, 2), (1, 2), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(1, 6), (2, 6), (2, 7), (1, 7)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 4), (4, 4), (4, 5), (3, 5)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (4, 3)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        'Y':[[(0, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8)], [(1, 5), (2, 5), (2, 6), (2, 7), (1, 7), (1, 6)], [(2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (2, 5), (2, 4), (2, 3)], [(3, 5), (4, 5), (4, 6), (4, 7), (3, 7), (3, 6)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]],
        'Z':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (0, 4), (0, 3)], [(0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9), (1, 9), (0, 9)], [(1, 4), (2, 4), (2, 5), (1, 5)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 6), (4, 6), (4, 7), (3, 7)]],
        '[':[[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (3, 8), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (1, 8), (1, 7), (1, 6), (1, 5), (1, 4), (1, 3)]],
        '\\':[[(0, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8)], [(1, 6), (2, 6), (2, 7), (1, 7)], [(2, 5), (3, 5), (3, 6), (2, 6)], [(3, 4), (4, 4), (4, 5), (3, 5)], [(4, 2), (5, 2), (5, 3), (5, 4), (4, 4), (4, 3)]],
        ']':[[(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (3, 9), (2, 9), (1, 9), (1, 8), (2, 8), (3, 8), (3, 7), (3, 6), (3, 5), (3, 4), (3, 3), (2, 3), (1, 3)]],
        '^':[[(0, 6), (1, 6), (1, 7), (0, 7)], [(1, 7), (2, 7), (2, 8), (1, 8)], [(2, 8), (3, 8), (3, 9), (2, 9)], [(3, 7), (4, 7), (4, 8), (3, 8)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        '_':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)]],
        '`':[[(1, 8), (2, 8), (2, 9), (1, 9)], [(2, 7), (3, 7), (3, 8), (2, 8)]],
        'a':[[(0, 3), (1, 3), (1, 4), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (3, 5), (2, 5), (1, 5), (1, 4), (2, 4), (3, 4), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'b':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4)]],
        'c':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'd':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8), (4, 7), (3, 7), (2, 7), (1, 7), (1, 6), (2, 6), (3, 6), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)]],
        'e':[[(0, 3), (1, 3), (1, 4), (2, 4), (3, 4), (4, 4), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (3, 5), (2, 5), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'f':[[(0, 5), (1, 5), (1, 4), (1, 3), (1, 2), (2, 2), (2, 3), (2, 4), (2, 5), (3, 5), (4, 5), (4, 6), (3, 6), (2, 6), (2, 7), (2, 8), (1, 8), (1, 7), (1, 6), (0, 6)], [(2, 8), (3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9), (2, 9)]],
        'g':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 0), (2, 0), (3, 0), (4, 0), (4, 1), (3, 1), (2, 1), (1, 1)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 1), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'h':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3)]],
        'i':[[(1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (2, 7), (1, 7)], [(2, 8), (3, 8), (3, 9), (2, 9)]],
        'j':[[(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1)], [(1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (2, 1), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (2, 7), (1, 7)], [(2, 8), (3, 8), (3, 9), (2, 9)]],
        'k':[[(0, 2), (1, 2), (1, 3), (1, 4), (2, 4), (3, 4), (3, 5), (2, 5), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 9), (0, 8), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        'l':[[(1, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (1, 9)]],
        'm':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3)]],
        'n':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(4, 2), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3)]],
        'o':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4)]],
        'p':[[(0, 0), (1, 0), (1, 1), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3), (0, 2), (0, 1)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4)]],
        'q':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 1), (4, 0), (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)], [(1, 6), (2, 6), (3, 6), (4, 6), (4, 7), (3, 7), (2, 7), (1, 7)]],
        'r':[[(0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 5), (3, 6), (2, 6), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4), (0, 3)], [(3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7)]],
        's':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)], [(0, 5), (1, 5), (1, 6), (0, 6)], [(1, 4), (2, 4), (3, 4), (4, 4), (4, 5), (3, 5), (2, 5), (1, 5)], [(1, 6), (2, 6), (3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (2, 7), (1, 7)], [(4, 3), (5, 3), (5, 4), (4, 4)]],
        't':[[(1, 6), (2, 6), (2, 5), (2, 4), (2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (2, 7), (1, 7)], [(3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3)]],
        'u':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'v':[[(0, 5), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6)], [(1, 3), (2, 3), (2, 4), (2, 5), (1, 5), (1, 4)], [(2, 2), (3, 2), (3, 3), (2, 3)], [(3, 3), (4, 3), (4, 4), (4, 5), (3, 5), (3, 4)], [(4, 5), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6)]],
        'w':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 2), (2, 2), (2, 3), (1, 3)], [(2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (2, 6), (2, 5), (2, 4)], [(3, 2), (4, 2), (4, 3), (3, 3)], [(4, 3), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5), (4, 4)]],
        'x':[[(0, 2), (1, 2), (1, 3), (0, 3)], [(0, 6), (1, 6), (1, 7), (0, 7)], [(1, 3), (2, 3), (2, 4), (1, 4)], [(1, 5), (2, 5), (2, 6), (1, 6)], [(2, 4), (3, 4), (3, 5), (2, 5)], [(3, 3), (4, 3), (4, 4), (3, 4)], [(3, 5), (4, 5), (4, 6), (3, 6)], [(4, 2), (5, 2), (5, 3), (4, 3)], [(4, 6), (5, 6), (5, 7), (4, 7)]],
        'y':[[(0, 3), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (0, 7), (0, 6), (0, 5), (0, 4)], [(1, 0), (2, 0), (3, 0), (4, 0), (4, 1), (3, 1), (2, 1), (1, 1)], [(1, 2), (2, 2), (3, 2), (4, 2), (4, 1), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (4, 7), (4, 6), (4, 5), (4, 4), (4, 3), (3, 3), (2, 3), (1, 3)]],
        'z':[[(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3), (2, 3), (2, 4), (1, 4), (1, 3), (0, 3)], [(0, 6), (1, 6), (2, 6), (3, 6), (3, 5), (4, 5), (4, 6), (5, 6), (5, 7), (4, 7), (3, 7), (2, 7), (1, 7), (0, 7)], [(2, 4), (3, 4), (3, 5), (2, 5)]],
        '{':[[(1, 5), (2, 5), (2, 4), (2, 3), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (2, 8), (2, 7), (2, 6), (1, 6)], [(3, 2), (4, 2), (5, 2), (5, 3), (4, 3), (3, 3)], [(3, 8), (4, 8), (5, 8), (5, 9), (4, 9), (3, 9)]],
        '|':[[(2, 2), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (2, 9), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4), (2, 3)]],
        '}':[[(0, 2), (1, 2), (2, 2), (2, 3), (1, 3), (0, 3)], [(0, 8), (1, 8), (2, 8), (2, 9), (1, 9), (0, 9)], [(2, 3), (3, 3), (3, 4), (3, 5), (4, 5), (4, 6), (3, 6), (3, 7), (3, 8), (2, 8), (2, 7), (2, 6), (2, 5), (2, 4)]],
        '~':[[(0, 6), (1, 6), (1, 7), (1, 8), (0, 8), (0, 7)], [(1, 8), (2, 8), (2, 9), (1, 9)], [(2, 7), (3, 7), (3, 8), (2, 8)], [(3, 6), (4, 6), (4, 7), (3, 7)], [(4, 7), (5, 7), (5, 8), (5, 9), (4, 9), (4, 8)]]}

    def __init__(self, text, size, position=(0, 0), horizontal=True, angle=0, layer=0, datatype=0) :
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
            if text[jj] == '\n':
                if horizontal:
                    posY -= 11
                    posX = 0
                else:
                    posX += 8
                    posY = 0
            elif text[jj] == '\t':
                if horizontal:
                    posX = posX + 32 - (posX + 8) % 32
                else:
                    posY = posY - 11 - (posY - 22) % 44
            else:
                if Text._font.has_key(text[jj]):
                    for p in Text._font[text[jj]]:
                        polygon = p[:]
                        for ii in range(len(polygon)):
                            xp = text_multiplier * (posX + polygon[ii][0])
                            yp = text_multiplier * (posY + polygon[ii][1])
                            polygon[ii] = (position[0] + xp * ca - yp * sa, position[1] + xp * sa + yp * ca)
                        self.polygons.append(numpy.array(polygon))
                if horizontal:
                    posX += 8
                else:
                    posY -= 11
        self.layers = [layer] * len(self.polygons)
        self.datatypes = [datatype] * len(self.polygons)

    def __str__(self):
        return "Text ({} polygons, {} vertices, layers {}, datatypes {})".format(len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))


class Path(PolygonSet):
    """
    Series of geometric objects that form a path or a collection of parallel
    paths.

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
    direction : {'+x', '-x', '+y', '-y'} or number
        Direction or angle (in *radians*) the path points to.
    distance : number
        Distance between the centers of adjacent paths.
    length : number
        Length of the central path axis. If only one path is created, this
        is the real length of the path.
    """

    def __init__(self, width, initial_point=(0, 0), number_of_paths=1, distance=0):
        self.x = initial_point[0]
        self.y = initial_point[1]
        self.w = width * 0.5
        self.n = number_of_paths
        self.direction = '+x'
        self.distance = distance
        self.length = 0.0
        self.polygons = []
        self.layers = []
        self.datatypes = []

    def __str__(self):
        if self.n > 1:
            return "Path (x{}, end at ({}, {}) towards {}, length {}, width {}, {} apart, {} polygons, {} vertices, layers {}, datatypes {})".format(self.n, self.x, self.y, self.direction, self.length, self.w * 2, self.distance, len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))
        else:
            return "Path (end at ({}, {}) towards {}, length {}, width {}, {} polygons, {} vertices, layers {}, datatypes {})".format(self.x, self.y, self.direction, self.length, self.w * 2, len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))

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
        out : ``Path``
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle)
        sa = numpy.array((-sa, sa))
        c0 = numpy.array(center)
        if self.direction.__class__ == ''.__class__:
            self.direction = {'+x': 0, '+y': 0.5, '-x': 1, '-y': -0.5}[self.direction] * numpy.pi
        self.direction += angle
        cur = numpy.array((self.x, self.y)) - c0
        self.x, self.y = cur * ca + cur[::-1] * sa + c0
        self.polygons = [(points - c0) * ca + (points - c0)[:,::-1] * sa + c0 for points in self.polygons]
        return self

    def segment(self, length, direction=None, final_width=None, final_distance=None, axis_offset=0, layer=0, datatype=0) :
        """
        Add a straight section to the path.

        Parameters
        ----------
        length : number
            Length of the section to add.
        direction : {'+x', '-x', '+y', '-y'} or number
            Direction or angle (in *radians*) of rotation of the segment.
        final_width : number
            If set, the paths of this segment will have their widths linearly
            changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from its
            current value to this one along this segment.
        axis_offset : number
            If set, the paths will be offset from their direction by this
            amount.
        layer : integer, list
            The GDSII layer numbers for the elements of each path. If the number
            of layers in the list is less than the number of paths, the list is
            repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0 and
            255). If the number of datatypes in the list is less than the number
            of paths, the list is repeated.

        Returns
        -------
        out : ``Path``
            This object.
        """
        if direction == None:
            direction = self.direction
        else:
            self.direction = direction
        if direction == '+x':
            ca = 1
            sa = 0
        elif direction == '-x':
            ca = -1
            sa = 0
        elif direction == '+y':
            ca = 0
            sa = 1
        elif direction == '-y':
            ca = 0
            sa = -1
        else:
            ca = numpy.cos(direction)
            sa = numpy.sin(direction)
        old_x = self.x
        old_y = self.y
        self.x += length*ca + axis_offset*sa
        self.y += length*sa - axis_offset*ca
        old_w = self.w
        old_distance = self.distance
        if not final_width == None:
            self.w = final_width * 0.5
        if not final_distance == None:
            self.distance = final_distance
        if (self.w != 0) or (old_w != 0):
            for ii in range(self.n):
                d0 = ii * self.distance - (self.n - 1) * self.distance * 0.5
                old_d0 = ii * old_distance - (self.n - 1) * old_distance * 0.5
                self.polygons.append(numpy.array([(old_x+(old_d0-old_w)*sa, old_y-(old_d0-old_w)*ca), (old_x+(old_d0+old_w)*sa, old_y-(old_d0+old_w)*ca), (self.x+(d0+self.w)*sa, self.y-(d0+self.w)*ca), (self.x+(d0-self.w)*sa, self.y-(d0-self.w)*ca)]))
                if self.w == 0:
                    self.polygons[-1] = self.polygons[-1][:-1,:]
                if old_w == 0:
                    self.polygons[-1] = self.polygons[-1][1:,:]
            self.length += numpy.sqrt(length ** 2 + axis_offset ** 2)
            if (layer.__class__ == [].__class__):
                self.layers += (layer * (self.n // len(layer) + 1))[:self.n]
            else:
                self.layers += [layer] * self.n
            if (datatype.__class__ == [].__class__) :
                self.datatypes += (datatype * (self.n // len(datatype) + 1))[ :self.n]
            else:
                self.datatypes += [datatype] * self.n
        return self

    def arc(self, radius, initial_angle, final_angle, number_of_points=0.01, max_points=199, final_width=None, final_distance=None, layer=0, datatype=0) :
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
        number_of_points : integer or float
            If integer: number of vertices that form the object (polygonal
            approximation). If float: approximate curvature resolution. The
            actual number of points is automatically calculated.
        max_points : integer
            if ``number_of_points > max_points``, the element will be fractured
            in smaller polygons with at most ``max_points`` each.
        final_width : number
            If set, the paths of this segment will have their widths linearly
            changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from its
            current value to this one along this segment.
        layer : integer, list
            The GDSII layer numbers for the elements of each path. If the number
            of layers in the list is less than the number of paths, the list is
            repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0 and
            255). If the number of datatypes in the list is less than the number
            of paths, the list is repeated.

        Returns
        -------
        out : ``Path``
            This object.

        Notes
        -----
        The GDSII specification supports only a maximum of 199 vertices per
        polygon.
        """
        warn = True
        cx = self.x - radius*numpy.cos(initial_angle)
        cy = self.y - radius*numpy.sin(initial_angle)
        self.x = cx + radius*numpy.cos(final_angle)
        self.y = cy + radius*numpy.sin(final_angle)
        if final_angle > initial_angle:
            self.direction = final_angle + numpy.pi * 0.5
        else:
            self.direction = final_angle - numpy.pi * 0.5
        old_w = self.w
        old_distance = self.distance
        if not final_width == None:
            self.w = final_width * 0.5
        if not final_distance == None:
            self.distance = final_distance
        if isinstance(number_of_points, float):
            number_of_points = 2 * int(abs((final_angle - initial_angle) * (radius + max(old_distance, self.distance) * (self.n - 1) * 0.5 + max(old_w, self.w)) / number_of_points) + 0.5) + 2
        number_of_points = max(number_of_points, 3)
        pieces = int(numpy.ceil(number_of_points / float(max_points)))
        number_of_points = number_of_points // pieces
        widths = numpy.linspace(old_w, self.w, pieces + 1)
        distances = numpy.linspace(old_distance, self.distance, pieces + 1)
        angles = numpy.linspace(initial_angle, final_angle, pieces + 1)
        if (self.w != 0) or (old_w != 0):
            for jj in range(pieces):
                for ii in range(self.n):
                    self.polygons.append(numpy.zeros((number_of_points, 2)))
                    r0 = radius + ii * distances[jj+1] - (self.n - 1) * distances[jj+1] * 0.5
                    old_r0 = radius + ii * distances[jj] - (self.n - 1) * distances[jj] * 0.5
                    pts2 = number_of_points // 2
                    pts1 = number_of_points - pts2
                    ang = numpy.linspace(angles[jj], angles[jj+1], pts1)
                    rad = numpy.linspace(old_r0 + widths[jj], r0 + widths[jj+1], pts1)
                    self.polygons[-1][:pts1,0] = numpy.cos(ang) * rad + cx
                    self.polygons[-1][:pts1,1] = numpy.sin(ang) * rad + cy
                    if widths[jj+1] == 0:
                        pts1 -= 1
                        pts2 += 1
                    if widths[jj] == 0:
                        self.polygons[-1][:pts1-1,:] = numpy.array(self.polygons[-1][1:pts1,:])
                        pts1 -= 1
                        pts2 += 1
                    ang = numpy.linspace(angles[jj+1], angles[jj], pts2)
                    rad = numpy.linspace(r0 - widths[jj+1], old_r0 - widths[jj], pts2)
                    if (rad[0] <= 0 or rad[-1] <= 0) and warn:
                        warnings.warn("[GDSPY] Path arc with width larger than radius created: possible self-intersecting polygon.", stacklevel=2)
                        warn = False
                    self.polygons[-1][pts1:,0] = numpy.cos(ang) * rad + cx
                    self.polygons[-1][pts1:,1] = numpy.sin(ang) * rad + cy
                self.length += abs((angles[jj+1] - angles[jj])*radius)
                if (layer.__class__ == [].__class__):
                    self.layers += (layer * (self.n // len(layer) + 1))[:self.n]
                else:
                    self.layers += [layer] * self.n
                if (datatype.__class__ == [].__class__) :
                    self.datatypes += (datatype * (self.n // len(datatype) + 1))[ :self.n]
                else:
                    self.datatypes += [datatype] * self.n
        return self

    def turn(self, radius, angle, number_of_points=0.01, max_points=199, final_width=None, final_distance=None, layer=0, datatype=0):
        """
        Add a curved section to the path.

        Parameters
        ----------
        radius : number
            Central radius of the section.
        angle : {'r', 'l', 'rr', 'll'} or number
            Angle (in *radians*) of rotation of the path. The values 'r' and 'l'
            represent 90-degree turns cw and ccw, respectively; the values 'rr'
            and 'll' represent analogous 180-degree turns.
        number_of_points : integer or float
            If integer: number of vertices that form the object (polygonal
            approximation). If float: approximate curvature resolution. The
            actual number of points is automatically calculated.
        max_points : integer
            if ``number_of_points > max_points``, the element will be fractured
            in smaller polygons with at most ``max_points`` each.
        final_width : number
            If set, the paths of this segment will have their widths linearly
            changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from its
            current value to this one along this segment.
        layer : integer, list
            The GDSII layer numbers for the elements of each path. If the number
            of layers in the list is less than the number of paths, the list is
            repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0 and
            255). If the number of datatypes in the list is less than the number
            of paths, the list is repeated.

        Returns
        -------
        out : ``Path``
            This object.

        Notes
        -----
        The GDSII specification supports only a maximum of 199 vertices per
        polygon.
        """
        exact = True
        if angle == 'r':
            delta_i = _halfpi
            delta_f = 0
        elif angle == 'rr':
            delta_i = _halfpi
            delta_f = -delta_i
        elif angle == 'l':
            delta_i = -_halfpi
            delta_f = 0
        elif angle == 'll':
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
        if self.direction == '+x':
            self.direction = 0
        elif self.direction == '-x':
            self.direction = numpy.pi
        elif self.direction == '+y':
            self.direction = _halfpi
        elif self.direction == '-y':
            self.direction = -_halfpi
        elif exact:
            exact = False
        self.arc(radius, self.direction + delta_i, self.direction + delta_f, number_of_points, max_points, final_width, final_distance, layer, datatype)
        if exact:
            self.direction = ['+x', '+y', '-x', '-y'][int(round(self.direction / _halfpi)) % 4]
        return self

    def parametric(self, curve_function, curve_derivative=None, number_of_evaluations=99, max_points=199, final_width=None, final_distance=None, layer=0, datatype=0):
        """
        Add a parametric curve to the path.

        Parameters
        ----------
        curve_function : function
            Function that defines the curve. Must be a function of one argument
            (that varies from 0 to 1) that returns a 2-element list, tuple or
            array (x, y).
        curve_derivative : function
            If set, it should be the derivative of the curve function. Must be a
            function of one argument (that varies from 0 to 1) that returns a
            2-element list, tuple or array (x,y). If ``None``, the derivative
            will be calculated numerically.
        number_of_evaluations : integer
            Number of points where the curve function will be evaluated. The
            final segment will have twice this number of points.
        max_points : integer
            If ``2 * number_of_evaluations > max_points``, the element will be
            fractured in smaller polygons with at most ``max_points`` each.
        final_width : number or function
            If set to a number, the paths of this segment will have their widths
            linearly changed from their current value to this one. If set to a
            function, it must be a function of one argument (that varies from 0
            to 1) and returns the width of the path.
        final_distance : number or function
            If set to ta number, the distance between paths is linearly change
            from its current value to this one. If set to a function, it must be
            a function of one argument (that varies from 0 to 1) and returns the
            width of the path.
        layer : integer, list
            The GDSII layer numbers for the elements of each path. If the number
            of layers in the list is less than the number of paths, the list is
            repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0 and
            255). If the number of datatypes in the list is less than the number
            of paths, the list is repeated.

        Returns
        -------
        out : ``Path``
            This object.

        Notes
        -----
        The norm of the vector returned by ``curve_derivative`` is not
        important. Only the direction is used.

        The GDSII specification supports only a maximum of 199 vertices per
        polygon.

        Examples
        --------
        >>> def my_parametric_curve(t):
        ...         return (2**t, t**2)
        >>> def my_parametric_curve_derivative(t):
        ...         return (0.69315 * 2**t, 2 * t)
        >>> my_path.parametric(my_parametric_curve,
        ...                    my_parametric_curve_derivative)
        """
        pieces = int(numpy.ceil(2 * number_of_evaluations / float(max_points)))
        number_of_evaluations = number_of_evaluations // pieces
        boundaries = numpy.linspace(0, 1, pieces + 1)
        if not callable(final_width):
            old_w = self.w
            if not final_width is None:
                self.w = final_width * 0.5
            final_width = lambda u: 2 * (old_w + u * (self.w - old_w))
        if not callable(final_distance):
            old_distance = self.distance
            if not final_distance is None:
                self.distance = final_distance
            final_distance = lambda u: old_distance + u * (self.distance - old_distance)
        if curve_derivative == None:
            def numerical_diff(t):
                delta = 0.5 / (number_of_evaluations - 1.0)
                if t == 0:
                    x0,y0 = curve_function(0)
                    x1,y1 = curve_function(delta)
                elif t == 1:
                    x0,y0 = curve_function(1 - delta)
                    x1,y1 = curve_function(1)
                else:
                    x0,y0 = curve_function(t - delta)
                    x1,y1 = curve_function(t + delta)
                return (x1 - x0, y1 - y0)
            curve_derivative = numerical_diff
        orgn = numpy.array((self.x, self.y))
        refl = numpy.array((-1,1))
        for kk in range(pieces):
            uu = numpy.linspace(boundaries[kk], boundaries[kk + 1], number_of_evaluations)
            width = numpy.array([final_width(u) for u in uu]).reshape(number_of_evaluations, 1) * 0.5
            dist = numpy.array([final_distance(u) for u in uu]).reshape(number_of_evaluations, 1)
            x0 = numpy.array([curve_function(u) for u in uu]) + orgn
            dx = numpy.array([curve_derivative(u) for u in uu])
            dx = dx[:,::-1] * refl / numpy.sqrt((dx * dx).sum(1)).reshape(number_of_evaluations, 1)
            self.length += numpy.sqrt(((x0[1:,:] - x0[:-1,:])**2).sum(1)).sum()
            for ii in range(self.n):
                p1 = x0 + dx * (dist * (ii - (self.n - 1) * 0.5) + width)
                p2 = (x0 + dx * (dist * (ii - (self.n - 1) * 0.5) - width))[::-1,:]
                if width[number_of_evaluations - 1,0] == 0:
                    p2 = p2[1:]
                if width[0,0] == 0:
                    p1 = p1[1:]
                self.polygons.append(numpy.concatenate((p1, p2)))
            if (layer.__class__ == [].__class__):
                self.layers += (layer * (self.n // len(layer) + 1))[:self.n]
            else:
                self.layers += [layer] * self.n
            if (datatype.__class__ == [].__class__) :
                self.datatypes += (datatype * (self.n // len(datatype) + 1))[:self.n]
            else:
                self.datatypes += [datatype] * self.n
        self.x = x0[-1,0]
        self.y = x0[-1,1]
        self.direction = numpy.arctan2(-dx[-1,0], dx[-1,1])
        return self


class L1Path(PolygonSet):
    """
    Series of geometric objects that form a path or a collection of parallel
    paths with Manhattan geometry.

    Parameters
    ----------
    initial_point : array-like[2]
        Starting position of the path.
    direction : {'+x', '+y', '-x', '-y'}
        Starting direction of the path.
    width : number
        The initial width of each path.
    length : array-like
        Lengths of each section to add.
    turn : array-like
        Direction to turn before each section. The sign indicate the turn
        direction (ccw is positive), and the modulus is a multiplicative
        factor for the path width after each turn. Must have 1 element less
        then ``length``.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number
        Distance between the centers of adjacent paths.
    layer : integer, list
        The GDSII layer numbers for the elements of each path. If the number
        of layers in the list is less than the number of paths, the list is
        repeated.
    datatype : integer, list
        The GDSII datatype for the elements of each path (between 0 and
        255). If the number of datatypes in the list is less than the number
        of paths, the list is repeated.

    Returns
    -------
    out : ``L1Path``
        This object.

    Attributes
    ----------
    x : number
        Final position of the path in the x direction.
    y : number
        Final position of the path in the y direction.
    direction : {'+x', '-x', '+y', '-y'} or number
        Direction or angle (in *radians*) the path points to. The numerical
        angle is returned only after a rotation of the object.

    Examples
    --------
    >>> length = [10, 30, 15, 15, 15, 15, 10]
    >>> turn = [1, -1, -1, 3, -1, 1]
    >>> l1path = gdspy.L1Path((0, 0), '+x', 2, length, turn)
    >>> myCell.add(l1path)
    """
    def __init__(self, initial_point, direction, width, length, turn, number_of_paths=1, distance=0, max_points=199, layer=0, datatype=0):
        if (layer.__class__ != [].__class__):
            layer = [layer]
        if (datatype.__class__ != [].__class__) :
            datatype = [datatype]
        layer = (layer * (number_of_paths // len(layer) + 1))[:number_of_paths]
        datatype = (datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths]
        w = width * 0.5
        points = max_points // 2 - 1
        paths = [[[], []] for ii in range(number_of_paths)]
        self.polygons = []
        self.layers = []
        self.datatypes = []
        self.x = initial_point[0]
        self.y = initial_point[1]
        if direction == '+x':
            direction = 0
            for ii in range(number_of_paths):
                d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                paths[ii][0].append((initial_point[0], d0 + initial_point[1] - w))
                paths[ii][1].append((initial_point[0], d0 + initial_point[1] + w))
        elif direction == '+y':
            direction = 1
            for ii in range(number_of_paths):
                d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                paths[ii][0].append((d0 + initial_point[0] + w, initial_point[1]))
                paths[ii][1].append((d0 + initial_point[0] - w, initial_point[1]))
        elif direction == '-x':
            direction = 2
            for ii in range(number_of_paths):
                d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                paths[ii][0].append((initial_point[0], d0 + initial_point[1] + w))
                paths[ii][1].append((initial_point[0], d0 + initial_point[1] - w))
        elif direction == '-y':
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
                    paths[ii][0].append((self.x + length[jj] - (d0 - w) * turn[jj], paths[ii][0][-1][1]))
                    paths[ii][1].append((self.x + length[jj] - (d0 + w) * turn[jj], paths[ii][1][-1][1]))
                self.x += length[jj]
            elif direction == 1:
                for ii in range(number_of_paths):
                    d0 = ii * distance - (number_of_paths - 1) * distance * 0.5
                    paths[ii][0].append((paths[ii][0][-1][0], self.y + length[jj] - (d0 - w) * turn[jj]))
                    paths[ii][1].append((paths[ii][1][-1][0], self.y + length[jj] - (d0 + w) * turn[jj]))
                self.y += length[jj]
            elif direction == 2:
                for ii in range(number_of_paths):
                    d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                    paths[ii][0].append((self.x - length[jj] - (d0 + w) * turn[jj], paths[ii][0][-1][1]))
                    paths[ii][1].append((self.x - length[jj] - (d0 - w) * turn[jj], paths[ii][1][-1][1]))
                self.x -= length[jj]
            elif direction == 3:
                for ii in range(number_of_paths):
                    d0 = (number_of_paths - 1) * distance * 0.5 - ii * distance
                    paths[ii][0].append((paths[ii][0][-1][0], self.y - length[jj] - (d0 + w) * turn[jj]))
                    paths[ii][1].append((paths[ii][1][-1][0], self.y - length[jj] - (d0 - w) * turn[jj]))
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
                    self.polygons.append(numpy.array(p[0][:-1] + [p0, p1] + p[1][-2::-1]))
                    p[0] = [p0, p[0][-1]]
                    p[1] = [p1, p[1][-1]]
                self.layers += layer
                self.datatypes += datatype
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
        self.direction = ['+x', '+y', '-x', '-y'][direction]
        self.polygons += [numpy.array(p[0] + p[1][::-1]) for p in paths]
        self.layers += layer
        self.datatypes += datatype

    def __str__(self):
        return "L1Path (end at ({}, {}) towards {}, {} polygons, {} vertices, layers {}, datatypes {})".format(self.x, self.y, self.direction, len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))

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
        out : ``L1Path``
            This object.
        """
        ca = numpy.cos(angle)
        sa = numpy.sin(angle)
        sa = numpy.array((-sa, sa))
        c0 = numpy.array(center)
        self.direction = {'+x': 0, '+y': 0.5, '-x': 1, '-y': -0.5}[self.direction] * numpy.pi
        self.direction += angle
        cur = numpy.array((self.x, self.y)) - c0
        self.x, self.y = cur * ca + cur[::-1] * sa + c0
        self.polygons = [(points - c0) * ca + (points - c0)[:,::-1] * sa + c0 for points in self.polygons]
        return self


class PolyPath(PolygonSet):
    """
    Series of geometric objects that form a polygonal path or a collection
    of parallel polygonal paths.

    Parameters
    ----------
    points : array-like[N][2]
        Endpoints of each path segment.
    width : number or array-like[N]
        Width of the path. If an array is given, width at each endpoint.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number or array-like[N]
        Distance between the centers of adjacent paths. If an array is
        given, distance at each endpoint.
    corners : positive integer
        Type of joins: 0 - miter join, 1 - bevel join
    ends : positive integer
        Type of path ends: 0 - no extension, 1 - rounded, 2 - extended by
        half width
    max_points : integer
        The paths will be fractured in polygons with at most
        ``max_points``.
    layer : integer, list
        The GDSII layer numbers for the elements of each path. If the number
        of layers in the list is less than the number of paths, the list is
        repeated.
    datatype : integer, list
        The GDSII datatype for the elements of each path (between 0 and
        255). If the number of datatypes in the list is less than the number
        of paths, the list is repeated.

    Returns
    -------
    out : ``PolyPath``
        This object.

    Notes
    -----
    The bevel join will give  strange results if the number of paths is
    greater than 1.
    """
    def __init__(self, points, width, number_of_paths=1, distance=0, corners=0, ends=0, max_points=199, layer=0, datatype=0):
        if (layer.__class__ != [].__class__):
            layer = [layer]
        if (datatype.__class__ != [].__class__) :
            datatype = [datatype]
        if hasattr(width, '__iter__'):
            width = numpy.array(width) * 0.5
        else:
            width = numpy.array([width * 0.5])
        len_w = len(width)
        if not hasattr(distance, '__iter__'):
            distance = (distance,)
        len_d = len(distance)
        self.polygons = []
        self.layers = []
        self.datatypes = []
        points = numpy.array(points, dtype=float)
        if ends == 2:
            v = points[0,:] - points[1,:]
            v /= numpy.sqrt(numpy.sum(v * v))
            points[0,:] = points[0,:] + v * width[0]
            v = points[-1,:] - points[-2,:]
            v /= numpy.sqrt(numpy.sum(v * v))
            points[-1,:] = points[-1,:] + v * width[(points.shape[0] - 1) % len_w]
        elif ends == 1:
            v0 = points[1,:] - points[0,:]
            angle0 = numpy.arctan2(v0[1], v0[0]) + _halfpi
            v0 = numpy.array((-v0[1], v0[0])) / numpy.sqrt(numpy.sum(v0 * v0))
            d0 = 0.5 * (number_of_paths - 1) * distance[0]
            v1 = points[-1,:] - points[-2,:]
            angle1 = numpy.arctan2(v1[1], v1[0]) - _halfpi
            v1 = numpy.array((-v1[1], v1[0])) / numpy.sqrt(numpy.sum(v1 * v1))
            j1w = (points.shape[0] - 1) % len_w
            j1d = (points.shape[0] - 1) % len_d
            d1 = 0.5 * (number_of_paths - 1) * distance[j1d]
            self.polygons.extend((Round(points[0,:] + v0 * (ii * distance[0] - d0), width[0], initial_angle=angle0, final_angle=angle0 + numpy.pi, number_of_points=33).polygons[0] for ii in range(number_of_paths)))
            self.polygons.extend((Round(points[-1,:] + v1 * (ii * distance[j1d] - d1), width[j1w], initial_angle=angle1, final_angle=angle1 + numpy.pi, number_of_points=33).polygons[0]) for ii in range(number_of_paths))
            self.layers += ((layer * (number_of_paths // len(layer) + 1))[:number_of_paths]) * 2
            self.datatypes += ((datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths]) * 2
        v = points[1,:] - points[0,:]
        v = numpy.array((-v[1], v[0])) / numpy.sqrt(numpy.sum(v * v))
        d0 = 0.5 * (number_of_paths - 1) * distance[0]
        d1 = 0.5 * (number_of_paths - 1) * distance[1 % len_d]
        paths = [[[points[0,:] + (ii * distance[0] - d0 - width[0]) * v], [points[0,:] + (ii * distance[0] - d0 + width[0]) * v]] for ii in range(number_of_paths)]
        p1 = [(points[1,:] + (ii * distance[1 % len_d] - d1 - width[1 % len_w]) * v, points[1,:] + (ii * distance[1 % len_d] - d1 + width[1 % len_w]) * v) for ii in range(number_of_paths)]
        for jj in range(1, points.shape[0] - 1):
            j0d = jj % len_d
            j0w = jj % len_w
            j1d = (jj + 1) % len_d
            j1w = (jj + 1) % len_w
            v = points[jj + 1,:] - points[jj,:]
            v = numpy.array((-v[1], v[0])) / numpy.sqrt(numpy.sum(v * v))
            d0 = d1
            d1 = 0.5 * (number_of_paths - 1) * distance[j1d]
            p0 = p1
            p1 = []
            pp = []
            for ii in range(number_of_paths):
                pp.append((points[jj,:] + (ii * distance[j0d] - d0 - width[j0w]) * v, points[jj,:] + (ii * distance[j0d] - d0 + width[j0w]) * v))
                p1.append((points[jj + 1,:] + (ii * distance[j1d] - d1 - width[j1w]) * v, points[jj + 1,:] + (ii * distance[j1d] - d1 + width[j1w]) * v))
                for kk in (0, 1):
                    p0m = paths[ii][kk][-1] - p0[ii][kk]
                    p1p = pp[ii][kk] - p1[ii][kk]
                    vec = p0m[0] * p1p[1] - p1p[0] * p0m[1]
                    if abs(vec) > 1e-30:
                        p = numpy.array((1, -1)) * (p0m * p1p[::-1] * p1[ii][kk] - p1p * p0m[::-1] * p0[ii][kk] + p0m * p1p * (p0[ii][kk][::-1] - p1[ii][kk][::-1])) / vec
                        if corners > 0 and numpy.sum((p - pp[ii][kk]) * p1p) > 0 and numpy.sum((p - p0[ii][kk]) * p0m) < 0:
                            paths[ii][kk].append(p0[ii][kk])
                            paths[ii][kk].append(pp[ii][kk])
                        else:
                            paths[ii][kk].append(p)
                if len(paths[ii][0]) + len(paths[ii][1]) + 3 > max_points:
                    if numpy.sum((paths[ii][0][0] - paths[ii][1][0]) ** 2) == 0:
                        paths[ii][1] = paths[ii][1][1:]
                    if numpy.sum((paths[ii][0][-1] - paths[ii][1][-1]) ** 2) == 0:
                        self.polygons.append(numpy.array(paths[ii][0] + paths[ii][1][-2::-1]))
                    else:
                        self.polygons.append(numpy.array(paths[ii][0] + paths[ii][1][::-1]))
                    paths[ii][0] = paths[ii][0][-1:]
                    paths[ii][1] = paths[ii][1][-1:]
                    self.layers.append(layer[ii % len(layer)])
                    self.datatypes.append(datatype[ii % len(datatype)])
        for ii in range(number_of_paths):
            if numpy.sum((paths[ii][0][0] - paths[ii][1][0]) ** 2) == 0:
                paths[ii][1] = paths[ii][1][1:]
            if numpy.sum((p1[ii][0] - p1[ii][1]) ** 2) != 0:
                paths[ii][0].append(p1[ii][0])
            paths[ii][1].append(p1[ii][1])
        self.polygons += [numpy.array(p[0] + p[1][::-1]) for p in paths]
        self.layers += (layer * (number_of_paths // len(layer) + 1))[:number_of_paths]
        self.datatypes += (datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths]

    def __str__(self):
        return "PolyPath ({} polygons, {} vertices, layers {}, datatypes {})".format(len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))


class Label:
    """
    Text that can be used to label parts of the geometry or display
    messages. The text does not create additional geometry, it's meant for
    display and labeling purposes only.

    Parameters
    ----------
    text : string
        The text of this label.
    position : array-like[2]
        Text anchor position.
    anchor : 'n', 's', 'e', 'w', 'o', 'ne', 'nw',...
        Position of the anchor relative to the text.
    rotation : number
        Angle of rotation of the label (in *degrees*).
    magnification : number
        Magnification factor for the label.
    x_reflection : bool
        If ``True``, the label is reflected parallel to the x direction
        before being rotated (not supported by LayoutViewer).
    layer : integer
        The GDSII layer number for these elements.
    texttype : integer
        The GDSII text type for the label (between 0 and 63).

    Examples
    --------
    >>> label = gdspy.Label('Sample label', (10, 0), 'sw')
    >>> myCell.add(label)
    """

    _anchor = {'nw':0,      'top left':0,        'upper left':0,
               'n':1,       'top center':1,      'upper center':1,
               'ne':2,      'top right':2,       'upper right':2,
               'w':4,       'middle left':4,
               'o':5,       'middle center':5,
               'e':6,       'middle right':6,
               'sw':8,      'bottom left':8,     'lower left':8,
               's':9,       'bottom center':9,       'lower center':9,
               'se':10, 'bottom right':10,   'lower right':10}

    def __init__(self, text, position, anchor='o', rotation=None, magnification=None, x_reflection=False, layer=0, texttype=0):
        self.layer = layer
        self.text = text
        self.position = position
        try:
            self.anchor = Label._anchor[anchor.lower()]
        except:
            warnings.warn("[GDSPY] Label anchors must be one of: '" + "', '".join(Label._anchor.keys()) + "'.", stacklevel=2)
            self.anchor = 0
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection
        self.texttype = texttype

    def __str__(self):
        return "Label (\"{0}\", at ({1[0]}, {1[1]}), rotation {2}, magnification {3}, reflection {4}, layer {5}, texttype {6})".format(self.text, self.position, self.rotation, self.magnification, self.x_reflection, self.layer, self.texttype)

    def to_gds(self, multiplier):
        """
        Convert this label to a GDSII structure.

        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            structure.

        Returns
        -------
        out : string
            The GDSII binary string that represents this label.
        """
        text = self.text
        if len(text)%2 != 0:
            text = text + '\0'
        data = struct.pack('>11h', 4, 0x0C00, 6, 0x0D02, self.layer, 6, 0x1602, self.texttype, 6, 0x1701, self.anchor)
        if not (self.rotation is None) or not (self.magnification is None) or self.x_reflection:
            word = 0
            values = b''
            if self.x_reflection:
                word += 0x8000
            if not (self.magnification is None):
                word += 0x0004
                values += struct.pack('>2h', 12, 0x1B05) + _eight_byte_real(self.magnification)
            if not (self.rotation is None):
                word += 0x0002
                values += struct.pack('>2h', 12, 0x1C05) + _eight_byte_real(self.rotation)
            data += struct.pack('>2hH', 6, 0x1A01, word) + values
        return data + struct.pack('>2h2l2h', 12, 0x1003, int(round(self.position[0] * multiplier)), int(round(self.position[1] * multiplier)), 4 + len(text), 0x1906) + text.encode('ascii') + struct.pack('>2h', 4, 0x1100)


class Cell:
    """
    Collection of elements, both geometric objects and references to other
    cells.

    Parameters
    ----------
    name : string
        The name of the cell.
    exclude_from_global : bool
        If ``True``, the cell will not be included in the global list of
        cells maintained by ``gdspy``.
    """

    cell_dict = {}
    """
    Dictionary containing all cells created, indexed by name.  This
    dictionary is updated automatically whenever a new ``Cell`` object is
    created without the ``exclude_from_global`` flag.
    """

    _bounding_boxes = {}
    """
    Dictionary cache containing bounding box information for each cell,
    along with rotation information for cell references.
    """

    def __init__(self, name, exclude_from_global=False):
        self.name = name
        self.elements = []
        self.labels = []
        self.bb_is_valid = False
        if name in Cell.cell_dict:
            raise ValueError("[GDSPY] A cell named {0} has already been created.".format(name))
        if not exclude_from_global:
            Cell.cell_dict[name] = self

    def __str__(self):
        return "Cell (\"{}\", {} elements, {} labels)".format(self.name, len(self.elements), len(self.labels))

    def to_gds(self, multiplier):
        """
        Convert this cell to a GDSII structure.

        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            structure.

        Returns
        -------
        out : string
            The GDSII binary string that represents this cell.
        """
        now = datetime.datetime.today()
        name = self.name
        if len(name)%2 != 0:
            name = name + '\0'
        return struct.pack('>16h', 28, 0x0502, now.year, now.month, now.day, now.hour, now.minute, now.second, now.year, now.month, now.day, now.hour, now.minute, now.second, 4 + len(name), 0x0606) + name.encode('ascii') + b''.join(element.to_gds(multiplier) for element in self.elements) + b''.join(label.to_gds(multiplier) for label in self.labels) + struct.pack('>2h', 4, 0x0700)

    def copy(self, name, exclude_from_global=False):
        """
        Creates a copy of this cell.

        Parameters
        ----------
        name : string
            The name of the cell.
        exclude_from_global : bool
            If ``True``, the cell will not be included in the global list of
            cells maintained by ``gdspy``.

        Returns
        -------
        out : ``Cell``
            The new copy of this cell.
        """
        new_cell = Cell(name, exclude_from_global)
        new_cell.elements = list(self.elements)
        new_cell.labels = list(self.labels)
        new_cell.bb_is_valid = False
        return new_cell

    def add(self, element):
        """
        Add a new element or list of elements to this cell.

        Parameters
        ----------
        element : object, list
            The element or list of elements to be inserted in this cell.

        Returns
        -------
        out : ``Cell``
            This cell.
        """
        if (element.__class__ == [].__class__):
            for e in element:
                if isinstance(e, Label):
                    self.labels.append(e)
                else:
                    self.elements.append(e)
        else:
            if isinstance(element, Label):
                self.labels.append(element)
            else:
                self.elements.append(element)
        self.bb_is_valid = False
        return self

    def area(self, by_spec=False):
        """
        Calculate the total area of the elements on this cell, including cell
        references and arrays.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the areas of each
            individual pair (layer, datatype).

        Returns
        -------
        out : number, dictionary
            Area of this cell.
        """
        if by_spec:
            cell_area = {}
            for element in self.elements:
                element_area = element.area(True)
                for ll in element_area.iterkeys():
                    if cell_area.has_key(ll):
                        cell_area[ll] += element_area[ll]
                    else:
                        cell_area[ll] = element_area[ll]
        else:
            cell_area = 0
            for element in self.elements:
                cell_area += element.area()
        return cell_area

    def get_layers(self):
        """
        Returns a list of layers in this cell.

        Returns
        -------
        out : list
            List of the layers used in this cell.
        """
        layers = []
        for element in self.elements:
            if isinstance(element, Polygon):
                if element.layer not in layers:
                    layers.append(element.layer)
            elif isinstance(element, PolygonSet):
                for layer in element.layers:
                    if layer not in layers:
                        layers.append(layer)
            elif isinstance(element, CellReference) or isinstance(element, CellArray):
                for layer in element.ref_cell.get_layers():
                    if layer not in layers:
                        layers.append(layer)
        for label in self.labels:
            if label.layer not in layers:
                layers.append(label.layer)
        return layers

    def get_datatypes(self):
        """
        Returns a list of datatypes in this cell.

        Returns
        -------
        out : list
            List of the datatypes used in this cell.
        """
        datatypes = []
        for element in self.elements:
            if isinstance(element, Polygon):
                if element.datatype not in datatypes:
                    datatypes.append(element.datatype)
            elif isinstance(element, PolygonSet):
                for datatype in element.datatypes:
                    if datatype not in datatypes:
                        datatypes.append(datatype)
            elif isinstance(element, CellReference) or isinstance(element, CellArray):
                for datatype in element.ref_cell.get_datatypes():
                    if datatype not in datatypes:
                        datatypes.append(datatype)
        return datatypes

    def get_bounding_box(self):
        """
        Returns the bounding box for this cell.

        Returns
        -------
        out : Numpy array[2,2] or ``None``
            Bounding box of this cell [[x_min, y_min], [x_max, y_max]], or
            ``None`` if the cell is empty.
        """
        if len(self.elements) == 0:
            return None
        self.bb_is_valid = (self.bb_is_valid and all([ref.bb_is_valid for ref in self.get_dependencies()]))
        if not (self.bb_is_valid and self in Cell._bounding_boxes):
            bb = numpy.array(((1e300, 1e300), (-1e300, -1e300)))
            all_polygons = []
            for element in self.elements:
                if isinstance(element, Polygon):
                    all_polygons.append(element.points)
                    #bb[0,0] = min(bb[0,0], element.points[:,0].min())
                    #bb[0,1] = min(bb[0,1], element.points[:,1].min())
                    #bb[1,0] = max(bb[1,0], element.points[:,0].max())
                    #bb[1,1] = max(bb[1,1], element.points[:,1].max())
                elif isinstance(element, PolygonSet):
                    all_polygons.extend(element.polygons)
                    #for points in element.polygons:
                    #    bb[0,0] = min(bb[0,0], points[:,0].min())
                    #    bb[0,1] = min(bb[0,1], points[:,1].min())
                    #    bb[1,0] = max(bb[1,0], points[:,0].max())
                    #    bb[1,1] = max(bb[1,1], points[:,1].max())
                elif isinstance(element, CellReference) or isinstance(element, CellArray):
                    element_bb = element.get_bounding_box()
                    if not element_bb is None:
                        bb[0,0] = min(bb[0,0], element_bb[0,0])
                        bb[0,1] = min(bb[0,1], element_bb[0,1])
                        bb[1,0] = max(bb[1,0], element_bb[1,0])
                        bb[1,1] = max(bb[1,1], element_bb[1,1])
            if len(all_polygons) > 0:
                all_points = numpy.concatenate(all_polygons).transpose()
                bb[0,0] = min(bb[0,0], all_points[0].min())
                bb[0,1] = min(bb[0,1], all_points[1].min())
                bb[1,0] = max(bb[1,0], all_points[0].max())
                bb[1,1] = max(bb[1,1], all_points[1].max())
            Cell._bounding_boxes[self] = bb
            self.bb_is_valid = True
        return Cell._bounding_boxes[self]

    def get_polygons(self, by_spec=False, depth=None):
        """
        Returns a list of polygons in this cell.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the polygons of
            each individual pair (layer, datatype).
        depth : integer or ``None``
            If not ``None``, defines from how many reference levels to retrieve
            polygons. References below this level will result in a bounding box.
            If ``by_spec`` is ``True`` the key wil be the name of this cell.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each polygon, or
            dictionary with the list of polygons (if ``by_spec`` is ``True``).
        """
        if not depth is None and depth < 0:
            bb = self.get_bounding_box()
            if bb is None:
                return {} if by_spec else []
            pts = [numpy.array([(bb[0,0],bb[0,1]), (bb[0,0],bb[1,1]), (bb[1,0],bb[1,1]), (bb[1,0],bb[0,1])])]
            polygons = {self.name: pts} if by_spec else pts
        else:
            if by_spec:
                polygons = {}
                for element in self.elements:
                    if isinstance(element, Polygon):
                        key = (element.layer, element.datatype)
                        if polygons.has_key(key):
                            polygons[key].append(numpy.array(element.points))
                        else:
                            polygons[key] = [numpy.array(element.points)]
                    elif isinstance(element, PolygonSet):
                        for ii in range(len(element.polygons)):
                            key = (element.layers[ii], element.datatypes[ii])
                            if polygons.has_key(key):
                                polygons[key].append(numpy.array(element.polygons[ii]))
                            else:
                                polygons[key] = [numpy.array(element.polygons[ii])]
                    else:
                        cell_polygons = element.get_polygons(True, None if depth is None else depth - 1)
                        for kk in cell_polygons.iterkeys():
                            if polygons.has_key(kk):
                                polygons[kk] += cell_polygons[kk]
                            else:
                                polygons[kk] = cell_polygons[kk]
            else:
                polygons = []
                for element in self.elements:
                    if isinstance(element, Polygon):
                        polygons.append(numpy.array(element.points))
                    elif isinstance(element, PolygonSet):
                        for points in element.polygons:
                            polygons.append(numpy.array(points))
                    else:
                        polygons += element.get_polygons(depth=None if depth is None else depth - 1)
        return polygons

    def get_dependencies(self):
        """
        Returns a list of the cells included in this cell as references.

        Returns
        -------
        out : list of ``Cell``
            List of the cells referenced by this cell.
        """
        dependencies = []
        for element in self.elements:
            if isinstance(element, CellReference) or isinstance(element, CellArray):
                dependencies.append(element.ref_cell)
        return dependencies

    def flatten(self, single_layer=None, single_datatype=None, verbose=True):
        """
        Flatten all ``CellReference`` and ``CellArray`` elements in this cell
        into real polygons, instead of references.

        Parameters
        ----------
        single_layer : integer or None
            If not ``None``, all polygons will be transfered to the layer
            indicated by this number.
        single_datatype : integer or None
            If not ``None``, all polygons will be transfered to the datatype
            indicated by this number.
        verbose : bool
            If False, warnings about the number of vertices of the polygon will
            be suppressed.

        Returns
        -------
        out : ``Cell``
            This cell.
        """
        if single_layer is None or single_datatype is None:
            poly_dic = self.get_polygons(True)
            self.elements = []
            if single_layer is None and single_datatype is None:
                for ld in poly_dic.iterkeys():
                    self.add(PolygonSet(poly_dic[ld], *ld, verbose=verbose))
            elif single_layer is None:
                for ld in poly_dic.iterkeys():
                    self.add(PolygonSet(poly_dic[ld], ld[0], single_datatype, verbose=verbose))
            else:
                for ld in poly_dic.iterkeys():
                    self.add(PolygonSet(poly_dic[ld], single_layer, ld[1], verbose=verbose))
        else:
            polygons = self.get_polygons()
            self.elements = []
            self.add(PolygonSet(polygons, single_layer, single_datatype, verbose=verbose))
        return self


class CellReference:
    """
    Simple reference to an existing cell.

    Parameters
    ----------
    ref_cell : ``Cell`` or string
        The referenced cell or its name.
    origin : array-like[2]
        Position where the reference is inserted.
    rotation : number
        Angle of rotation of the reference (in *degrees*).
    magnification : number
        Magnification factor for the reference.
    x_reflection : bool
        If ``True``, the reference is reflected parallel to the x direction
        before being rotated.
    """

    def __init__(self, ref_cell, origin=(0, 0), rotation=None, magnification=None, x_reflection=False):
        self.origin = origin
        self.ref_cell = Cell.cell_dict.get(ref_cell, ref_cell)
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection

    def __str__(self):
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        return "CellReference (\"{0}\", at ({1[0]}, {1[1]}), rotation {2}, magnification {3}, reflection {4})".format(name, self.origin, self.rotation, self.magnification, self.x_reflection)

    def __repr__(self):
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        return "CellReference(\"{0}\", ({1[0]}, {1[1]}), {2}, {3}, {4})".format(name, self.origin, self.rotation, self.magnification, self.x_reflection)

    def to_gds(self, multiplier):
        """
        Convert this object to a GDSII element.

        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            element.

        Returns
        -------
        out : string
            The GDSII binary string that represents this object.
        """
        name = self.ref_cell.name
        if len(name)%2 != 0:
            name = name + '\0'
        data = struct.pack('>4h', 4, 0x0A00, 4 + len(name), 0x1206) + name.encode('ascii')
        if not (self.rotation is None) or not (self.magnification is None) or self.x_reflection:
            word = 0
            values = b''
            if self.x_reflection:
                word += 0x8000
            if not (self.magnification is None):
                word += 0x0004
                values += struct.pack('>2h', 12, 0x1B05) + _eight_byte_real(self.magnification)
            if not (self.rotation is None):
                word += 0x0002
                values += struct.pack('>2h', 12, 0x1C05) + _eight_byte_real(self.rotation)
            data += struct.pack('>2hH', 6, 0x1A01, word) + values
        return data + struct.pack('>2h2l2h', 12, 0x1003, int(round(self.origin[0] * multiplier)), int(round(self.origin[1] * multiplier)), 4, 0x1100)

    def area(self, by_spec=False):
        """
        Calculate the total area of the referenced cell with the magnification
        factor included.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the areas of each
            individual pair (layer, datatype).

        Returns
        -------
        out : number, dictionary
            Area of this cell.
        """
        if self.magnification is None:
            return self.ref_cell.area(by_spec)
        else:
            if by_spec:
                factor = self.magnification * self.magnification
                cell_area = self.ref_cell.area(True)
                for kk in cell_area.iterkeys():
                    cell_area[kk] *= factor
                return cell_area
            else:
                return self.ref_cell.area() * self.magnification * self.magnification

    def get_polygons(self, by_spec=False, depth=None):
        """
        Returns a list of polygons created by this reference.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the polygons of
            each individual pair (layer, datatype).
        depth : integer or ``None``
            If not ``None``, defines from how many reference levels to retrieve
            polygons. References below this level will result in a bounding box.
            If ``by_spec`` is ``True`` the key will be the name of the
            referenced cell.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each polygon, or
            dictionary with the list of polygons (if ``by_spec`` is ``True``).
        """
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0)
            st = numpy.array([-st, st])
        if self.x_reflection:
            xrefl = numpy.array([1, -1], dtype=int)
        if self.magnification is not None:
            mag = numpy.array([self.magnification, self.magnification])
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        if by_spec:
            polygons = self.ref_cell.get_polygons(True, depth)
            for kk in polygons.iterkeys():
                for ii in range(len(polygons[kk])):
                    if self.x_reflection:
                        polygons[kk][ii] *= xrefl
                    if self.magnification is not None:
                        polygons[kk][ii] *= mag
                    if self.rotation is not None:
                        polygons[kk][ii] = polygons[kk][ii] * ct + polygons[kk][ii][:,::-1] * st
                    if self.origin is not None:
                        polygons[kk][ii] = polygons[kk][ii] + orgn
        else:
            polygons = self.ref_cell.get_polygons(depth=depth)
            for ii in range(len(polygons)):
                if self.x_reflection:
                    polygons[ii] *= xrefl
                if self.magnification is not None:
                    polygons[ii] *= mag
                if self.rotation is not None:
                    polygons[ii] = polygons[ii] * ct + polygons[ii][:,::-1] * st
                if self.origin is not None:
                    polygons[ii] = polygons[ii] + orgn
        return polygons

    def get_bounding_box(self):
        """
        Returns the bounding box for this reference.

        Returns
        -------
        out : Numpy array[2,2] or ``None``
            Bounding box of this cell [[x_min, y_min], [x_max, y_max]], or
            ``None`` if the cell is empty.
        """
        self.ref_cell.bb_is_valid = (self.ref_cell.bb_is_valid and all([ref.bb_is_valid for ref in self.ref_cell.get_dependencies()]))
        if self.rotation is None and self.magnification is None and self.x_reflection is None:
            key = self
        else:
            key = (self.ref_cell, self.rotation, self.magnification, self.x_reflection)
        if not (self.ref_cell.bb_is_valid and key in Cell._bounding_boxes):
            tmp = self.origin
            self.origin = None
            polygons = self.get_polygons()
            self.origin = tmp
            if len(polygons) == 0:
                bb = None
            else:
                all_points = numpy.concatenate(polygons).transpose()
                bb = numpy.array(((all_points[0].min(), all_points[1].min()), (all_points[0].max(), all_points[1].max())))
                #bb = numpy.array(((1e300, 1e300), (-1e300, -1e300)))
                #for points in polygons:
                #    bb[0,0] = min(bb[0,0], points[:,0].min())
                #    bb[0,1] = min(bb[0,1], points[:,1].min())
                #    bb[1,0] = max(bb[1,0], points[:,0].max())
                #    bb[1,1] = max(bb[1,1], points[:,1].max())
            Cell._bounding_boxes[key] = bb
            self.ref_cell.bb_is_valid = True
            for ref in self.ref_cell.get_dependencies():
                ref.bb_is_valid = True
        bb = Cell._bounding_boxes[key]
        if self.origin is None or bb is None:
            return bb
        else:
            return bb + numpy.array(((self.origin[0], self.origin[1]), (self.origin[0], self.origin[1])))


class CellArray:
    """
    Multiple references to an existing cell in an array format.

    Parameters
    ----------
    ref_cell : ``Cell`` or string
        The referenced cell or its name.
    columns : positive integer
        Number of columns in the array.
    rows : positive integer
        Number of columns in the array.
    spacing : array-like[2]
        distances between adjacent columns and adjacent rows.
    origin : array-like[2]
        Position where the cell is inserted.
    rotation : number
        Angle of rotation of the reference (in *degrees*).
    magnification : number
        Magnification factor for the reference.
    x_reflection : bool
        If ``True``, the reference is reflected parallel to the x direction
        before being rotated.
    """

    def __init__(self, ref_cell, columns, rows, spacing, origin=(0, 0), rotation=None, magnification=None, x_reflection=False):
        self.columns = columns
        self.rows = rows
        self.spacing = spacing
        self.origin = origin
        self.ref_cell = Cell.cell_dict.get(ref_cell, ref_cell)
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection

    def __str__(self):
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        return "CellArray (\"{0}\", {1} x {2}, at ({3[0]}, {3[1]}), spacing {4[0]} x {4[1]}, rotation {5}, magnification {6}, reflection {7})".format(name, self.columns, self.rows, self.origin, self.spacing, self.rotation, self.magnification, self.x_reflection)

    def __repr__(self):
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        return "CellArray(\"{0}\", {1}, {2}, ({4[0]}, {4[1]}), ({3[0]}, {3[1]}), {5}, {6}, {7})".format(name, self.columns, self.rows, self.origin, self.spacing, self.rotation, self.magnification, self.x_reflection)

    def to_gds(self, multiplier):
        """
        Convert this object to a GDSII element.

        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            element.

        Returns
        -------
        out : string
            The GDSII binary string that represents this object.
        """
        name = self.ref_cell.name
        if len(name)%2 != 0:
            name = name + '\0'
        data = struct.pack('>4h', 4, 0x0B00, 4 + len(name), 0x1206) + name.encode('ascii')
        x2 = self.origin[0] + self.columns * self.spacing[0]
        y2 = self.origin[1]
        x3 = self.origin[0]
        y3 = self.origin[1] + self.rows * self.spacing[1]
        if not (self.rotation is None) or not (self.magnification is None) or self.x_reflection:
            word = 0
            values = b''
            if self.x_reflection:
                word += 0x8000
                y3 = 2 * self.origin[1] - y3
            if not (self.magnification is None):
                word += 0x0004
                values += struct.pack('>2h', 12, 0x1B05) + _eight_byte_real(self.magnification)
            if not (self.rotation is None):
                word += 0x0002
                sa = numpy.sin(self.rotation * numpy.pi / 180.0)
                ca = numpy.cos(self.rotation * numpy.pi / 180.0)
                tmp = (x2 - self.origin[0]) * ca - (y2 - self.origin[1]) * sa + self.origin[0]
                y2 = (x2 - self.origin[0]) * sa + (y2 - self.origin[1]) * ca + self.origin[1]
                x2 = tmp
                tmp = (x3 - self.origin[0]) * ca - (y3 - self.origin[1]) * sa + self.origin[0]
                y3 = (x3 - self.origin[0]) * sa + (y3 - self.origin[1]) * ca + self.origin[1]
                x3 = tmp
                values += struct.pack('>2h', 12, 0x1C05) + _eight_byte_real(self.rotation)
            data += struct.pack('>2hH', 6, 0x1A01, word) + values
        return data + struct.pack('>6h6l2h', 8, 0x1302, self.columns, self.rows, 28, 0x1003, int(round(self.origin[0] * multiplier)), int(round(self.origin[1] * multiplier)), int(round(x2 * multiplier)), int(round(y2 * multiplier)), int(round(x3 * multiplier)), int(round(y3 * multiplier)), 4, 0x1100)

    def area(self, by_spec=False):
        """
        Calculate the total area of the cell array with the magnification factor
        included.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the areas of each
            individual pair (layer, datatype).

        Returns
        -------
        out : number, dictionary
            Area of this cell.
        """
        if self.magnification is None:
            factor = self.columns * self.rows
        else:
            factor = self.columns * self.rows * self.magnification * self.magnification
        if by_spec:
            cell_area = self.ref_cell.area(True)
            for kk in cell_area.iterkeys():
                cell_area[kk] *= factor
            return cell_area
        else:
            return self.ref_cell.area() * factor

    def get_polygons(self, by_spec=False, depth=None):
        """
        Returns a list of polygons created by this reference.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the polygons of
            each individual pair (layer, datatype).
        depth : integer or ``None``
            If not ``None``, defines from how many reference levels to retrieve
            polygons. References below this level will result in a bounding box.
            If ``by_spec`` is ``True`` the key will be name of the referenced
            cell.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each polygon, or
            dictionary with the list of polygons (if ``by_spec`` is ``True``).
        """
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0)
            st = numpy.array([-st, st])
        if self.magnification is not None:
            mag = numpy.array([self.magnification, self.magnification])
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        if self.x_reflection:
            xrefl = numpy.array([1, -1], dtype=int)
        if by_spec:
            cell_polygons = self.ref_cell.get_polygons(True, depth)
            polygons = {}
            for kk in cell_polygons.iterkeys():
                polygons[kk] = []
                for ii in range(self.columns):
                    for jj in range(self.rows):
                        spc = numpy.array([self.spacing[0] * ii, self.spacing[1] * jj])
                        for points in cell_polygons[kk]:
                            if self.magnification:
                                polygons[kk].append(points * mag + spc)
                            else:
                                polygons[kk].append(points + spc)
                            if self.x_reflection:
                                polygons[kk][-1] *= xrefl
                            if self.rotation is not None:
                                polygons[kk][-1] = polygons[kk][-1] * ct + polygons[kk][-1][:,::-1] * st
                            if self.origin is not None:
                                polygons[kk][-1] = polygons[kk][-1] + orgn
        else:
            cell_polygons = self.ref_cell.get_polygons(depth=depth)
            polygons = []
            for ii in range(self.columns):
                for jj in range(self.rows):
                    spc = numpy.array([self.spacing[0] * ii, self.spacing[1] * jj])
                    for points in cell_polygons:
                        if self.magnification:
                            polygons.append(points * mag + spc)
                        else:
                            polygons.append(points + spc)
                        if self.x_reflection:
                            polygons[-1] *= xrefl
                        if self.rotation is not None:
                            polygons[-1] = polygons[-1] * ct + polygons[-1][:,::-1] * st
                        if self.origin is not None:
                            polygons[-1] = polygons[-1] + orgn
        return polygons

    def get_bounding_box(self):
        """
        Returns the bounding box for this reference.

        Returns
        -------
        out : Numpy array[2,2] or ``None``
            Bounding box of this cell [[x_min, y_min], [x_max, y_max]], or
            ``None`` if the cell is empty.
        """
        self.ref_cell.bb_is_valid = (self.ref_cell.bb_is_valid and all([ref.bb_is_valid for ref in self.ref_cell.get_dependencies()]))
        key = (self.ref_cell, self.rotation, self.magnification, self.x_reflection, self.columns, self.rows, self.spacing[0], self.spacing[1])
        if not (self.ref_cell.bb_is_valid and key in Cell._bounding_boxes):
            tmp = self.origin
            self.origin = None
            polygons = self.get_polygons()
            self.origin = tmp
            if len(polygons) == 0:
                bb = None
            else:
                all_points = numpy.concatenate(polygons).transpose()
                bb = numpy.array(((all_points[0].min(), all_points[1].min()), (all_points[0].max(), all_points[1].max())))
                #bb = numpy.array(((1e300, 1e300), (-1e300, -1e300)))
                #for points in polygons:
                #    bb[0,0] = min(bb[0,0], points[:,0].min())
                #    bb[0,1] = min(bb[0,1], points[:,1].min())
                #    bb[1,0] = max(bb[1,0], points[:,0].max())
                #    bb[1,1] = max(bb[1,1], points[:,1].max())
            Cell._bounding_boxes[key] = bb
            self.ref_cell.bb_is_valid = True
            for ref in self.ref_cell.get_dependencies():
                ref.bb_is_valid = True
        bb = Cell._bounding_boxes[key]
        if self.origin is None or bb is None:
            return bb
        else:
            return bb + numpy.array(((self.origin[0], self.origin[1]), (self.origin[0], self.origin[1])))


class GdsImport:
    """
    Object used to import structures from a GDSII stream file.

    Parameters
    ----------
    infile : file or string
        GDSII stream file (or path) to be imported. It must be opened for
        reading in binary format.
    unit : number
        Unit (in *meters*) to use for the imported structures. If ``None``,
        the units used to create the GDSII file will be used.
    rename : dictionary
        Dictionary used to rename the imported cells. Keys and values must
        be strings.
    layers : dictionary
        Dictionary used to convert the layers in the imported cells. Keys
        and values must be integers.
    datatypes : dictionary
        Dictionary used to convert the datatypes in the imported cells. Keys
        and values must be integers.
    texttypes : dictionary
        Dictionary used to convert the text types in the imported cells.
        Keys and values must be integers.
    verbose: bool
        If False, suppresses warnings about unsupported elements in the
        imported file.  Also supresses polygon generation warnings.

    Attributes
    ----------
    cell_dict : dictionary
        Dictionary will all imported cells, indexed by name.

    Notes
    -----
    Not all features from the GDSII specification are currently supported. A
    warning will be produced if any unsuported features are found in the
    imported file.

    Examples
    --------
    >>> gdsii = gdspy.GdsImport('gdspy-sample.gds')
    >>> for cell_name in gdsii.cell_dict:
    ...     gdsii.extract(cell_name)
    """

    _record_name = ('HEADER', 'BGNLIB', 'LIBNAME', 'UNITS', 'ENDLIB', 'BGNSTR', 'STRNAME', 'ENDSTR', 'BOUNDARY', 'PATH', 'SREF', 'AREF', 'TEXT', 'LAYER', 'DATATYPE', 'WIDTH', 'XY', 'ENDEL', 'SNAME', 'COLROW', 'TEXTNODE', 'NODE', 'TEXTTYPE', 'PRESENTATION', 'SPACING', 'STRING', 'STRANS', 'MAG', 'ANGLE', 'UINTEGER', 'USTRING', 'REFLIBS', 'FONTS', 'PATHTYPE', 'GENERATIONS', 'ATTRTABLE', 'STYPTABLE', 'STRTYPE', 'ELFLAGS', 'ELKEY', 'LINKTYPE', 'LINKKEYS', 'NODETYPE', 'PROPATTR', 'PROPVALUE', 'BOX', 'BOXTYPE', 'PLEX', 'BGNEXTN', 'ENDTEXTN', 'TAPENUM', 'TAPECODE', 'STRCLASS', 'RESERVED', 'FORMAT', 'MASK', 'ENDMASKS', 'LIBDIRSIZE', 'SRFNAME', 'LIBSECUR')
    _unused_records = (0x05, 0x00, 0x01, 0x02, 0x034, 0x38)

    def __init__(self, infile, unit=None, rename={}, layers={}, datatypes={}, texttypes={}, verbose=True):
        self.cell_dict = {}
        self._incomplete = []
        if infile.__class__ == ''.__class__:
            infile = open(infile, 'rb')
            close = True
        else:
            close = False
        emitted_warnings = []
        record = self._read_record(infile)
        kwargs = {}
        create_element = None
        while record is not None:
            ## LAYER
            if record[0] == 0x0d:
                kwargs['layer'] = layers.get(record[1][0], record[1][0])
            ## DATATYPE
            elif record[0] == 0x0e:
                kwargs['datatype'] = datatypes.get(record[1][0], record[1][0])
            ## TEXTTYPE
            elif record[0] == 0x16:
                kwargs['texttype'] = texttypes.get(record[1][0], record[1][0])
            ## XY
            elif record[0] == 0x10:
                kwargs['xy'] = factor * record[1]
            ## WIDTH
            elif record[0] == 0x0f:
                kwargs['width'] = factor * abs(record[1][0])
                if record[1][0] < 0 and record[0] not in emitted_warnings:
                    warnings.warn("[GDSPY] Paths with absolute width value are not supported. Scaling these paths will also scale their width.", stacklevel=2)
                    emitted_warnings.append(record[0])
            ## ENDEL
            elif record[0] == 0x11:
                if create_element is not None:
                    cell.add(create_element(**kwargs))
                    create_element = None
                kwargs = {}
            ## BOUNDARY
            elif record[0] == 0x08:
                kwargs['verbose'] = verbose
                create_element = self._create_polygon
            ## PATH
            elif record[0] == 0x09:
                create_element = self._create_path
            ## TEXT
            elif record[0] == 0x0c:
                create_element = self._create_label
            ## SNAME
            elif record[0] == 0x12:
                if not str is bytes:
                    if record[1][-1] == 0:
                        record[1] = record[1][:-1].decode('ascii')
                    else:
                        record[1] = record[1].decode('ascii')
                kwargs['ref_cell'] = rename.get(record[1], record[1])
            ## COLROW
            elif record[0] == 0x13:
                kwargs['columns'] = record[1][0]
                kwargs['rows'] = record[1][1]
            ## STRANS
            elif record[0] == 0x1a:
                kwargs['x_reflection'] = ((int(record[1][0]) & 0x8000) > 0)
            ## MAG
            elif record[0] == 0x1b:
                kwargs['magnification'] = record[1][0]
            ## ANGLE
            elif record[0] == 0x1c:
                kwargs['rotation'] = record[1][0]
            ## SREF
            elif record[0] == 0x0a:
                create_element = self._create_reference
            ## AREF
            elif record[0] == 0x0b:
                create_element = self._create_array
            ## STRNAME
            elif record[0] == 0x06:
                if not str is bytes:
                    if record[1][-1] == 0:
                        record[1] = record[1][:-1].decode('ascii')
                    else:
                        record[1] = record[1].decode('ascii')
                name = rename.get(record[1], record[1])
                cell = Cell(name, exclude_from_global=True)
                self.cell_dict[name] = cell
            ## STRING
            elif record[0] == 0x19:
                if not str is bytes:
                    if record[1][-1] == 0:
                        kwargs['text'] = record[1][:-1].decode('ascii')
                    else:
                        kwargs['text'] = record[1].decode('ascii')
                else:
                    kwargs['text'] = record[1]
            ## ENDSTR
            elif record[0] == 0x07:
                cell = None
            ## UNITS
            elif record[0] == 0x03:
                if unit is None:
                    factor = record[1][0]
                else:
                    factor = record[1][1] / unit
            ## PRESENTATION
            elif record[0] == 0x17:
                kwargs['anchor'] = ['nw', 'n', 'ne', None, 'w', 'o', 'e', None, 'sw', 's', 'se'][int(record[1][0]) & 0x000f]
            ## PATHTYPE
            elif record[0] == 0x21:
                if record[1][0] > 2:
                    if verbose and 0x21 not in emitted_warnings:
                        warnings.warn("[GDSPY] Path ends with custom size are not supported by gds_import.", stacklevel=2)
                        emitted_warnings.append(0x21)
                else:
                    kwargs['ends'] = record[1][0]
            ## ENDLIB
            elif record[0] == 0x04:
                for ref in self._incomplete:
                    if ref.ref_cell in self.cell_dict:
                        ref.ref_cell = self.cell_dict[ref.ref_cell]
                    else:
                        ref.ref_cell = Cell.cell_dict.get(ref.ref_cell, ref.ref_cell)
            ## Not supported
            elif verbose and record[0] not in emitted_warnings and record[0] not in GdsImport._unused_records:
                warnings.warn("[GDSPY] Record type {0} is not supported by gds_import.".format(GdsImport._record_name[record[0]]), stacklevel=2)
                emitted_warnings.append(record[0])
            record = self._read_record(infile)
        if close:
            infile.close()

    def __str__(self):
        return "GdsImport (" + ", ".join([c for c in self.cell_dict]) + ")"

    def _read_record(self, stream):
        """
        Read a complete record from a GDSII stream file.

        Parameters
        ----------
        stream : file
            GDSII stream file to be imported.

        Returns
        -------
        out : 2-tuple
            Record type and data (as a numpy.array)
        """
        header = stream.read(4)
        if len(header) < 4:
            return None
        size, rec_type = struct.unpack('>HH', header)
        data_type = (rec_type & 0x00ff)
        rec_type = rec_type // 256
        data = None
        if size > 4:
            if data_type == 0x01:
                data = numpy.array(struct.unpack('>{0}H'.format((size - 4) // 2), stream.read(size - 4)), dtype='uint')
            elif data_type == 0x02:
                data = numpy.array(struct.unpack('>{0}h'.format((size - 4) // 2), stream.read(size - 4)), dtype=int)
            elif data_type == 0x03:
                data = numpy.array(struct.unpack('>{0}l'.format((size - 4) // 4), stream.read(size - 4)), dtype=int)
            elif data_type == 0x05:
                data = numpy.array([_eight_byte_real_to_float(stream.read(8)) for _ in range((size - 4) // 8)])
            else:
                data = stream.read(size - 4)
                if data[-1] == '\0':
                    data = data[:-1]
        return [rec_type, data]

    def _create_polygon(self, layer, datatype, xy, verbose):
        return Polygon(xy[:-2].reshape((xy.size // 2 - 1, 2)), layer, datatype, verbose)

    def _create_path(self, **kwargs):
        xy = kwargs.pop('xy')
        kwargs['points'] = xy.reshape((xy.size // 2, 2))
        return PolyPath(**kwargs)

    def _create_label(self, xy, width=None, ends=None, **kwargs):
        kwargs['position'] = xy
        return Label(**kwargs)

    def _create_reference(self, **kwargs):
        kwargs['origin'] = kwargs.pop('xy')
        ref = CellReference(**kwargs)
        if not isinstance(ref.ref_cell, Cell):
            self._incomplete.append(ref)
        return ref

    def _create_array(self, **kwargs):
        xy = kwargs.pop('xy')
        kwargs['origin'] = xy[0:2]
        if 'x_reflection' in kwargs:
            if 'rotation' in kwargs:
                sa = -numpy.sin(kwargs['rotation'] * numpy.pi / 180.0)
                ca = numpy.cos(kwargs['rotation'] * numpy.pi / 180.0)
                x2 = (xy[2] - xy[0]) * ca - (xy[3] - xy[1]) * sa + xy[0]
                y3 = (xy[4] - xy[0]) * sa + (xy[5] - xy[1]) * ca + xy[1]
            else:
                x2 = xy[2]
                y3 = xy[5]
            if kwargs['x_reflection']:
                y3 = 2 * xy[1] - y3
            kwargs['spacing'] = ((x2 - xy[0]) / kwargs['columns'], (y3 - xy[1]) / kwargs['rows'])
        else:
            kwargs['spacing'] = ((xy[2] - xy[0]) / kwargs['columns'], (xy[5] - xy[1]) / kwargs['rows'])
        ref = CellArray(**kwargs)
        if not isinstance(ref.ref_cell, Cell):
            self._incomplete.append(ref)
        return ref

    def extract(self, cell):
        """
        Extract a cell from the imported GDSII file and include it in the
        current scope, including referenced dependencies.

        Parameters
        ----------
        cell : ``Cell`` or string
            Cell or name of the cell to be extracted from the imported file.
            Referenced cells will be automatically extracted as well.

        Returns
        -------
        out : ``Cell``
            The extracted cell.
        """
        cell = self.cell_dict.get(cell, cell)
        for c in cell.get_dependencies():
            if c not in Cell.cell_dict.values():
                self.extract(c)
        if cell not in Cell.cell_dict.values():
            Cell.cell_dict[cell.name] = cell
        return cell

    def top_level(self):
        """
        Output the top level cells from the GDSII data.  Top level cells are
        those that are not referenced by any other cells.

        Outputs
        ----------
        out: List
            List of top level cells.
        """
        top = list(self.cell_dict.itervalues())
        for cell in self.cell_dict.itervalues():
            for dependency in cell.get_dependencies():
                if dependency in top:
                    top.remove(dependency)
        return top


def _chop(polygon, position, axis):
    """
    Slice polygon at a given position along a given axis.

    Parameters
    ----------
    polygon : array-like[N][2]
        Coordinates of the vertices of the polygon.
    position : number
        Position to perform the slicing operation along the specified axis.
    axis : 0 or 1
        Axis along which the polygon will be sliced.

    Returns
    -------
    out : tuple[2]
        Each element is a list of polygons (array-like[N][2]). The first
        list contains the polygons left before the slicing position, and the
        second, the polygons left after that position.
    """
    out_polygons = ([], [])
    polygon = list(polygon)
    while polygon[-1][axis] == position:
        polygon = [polygon[-1]] + polygon[:-1]
    cross = list(numpy.sign(numpy.array(polygon)[:, axis] - position))
    bnd = ([], [])
    i = 0
    while i < len(cross):
        if cross[i - 1] * cross[i] < 0:
            if axis == 0:
                polygon.insert(i, [position, polygon[i - 1][1] + (position - polygon[i - 1][0]) * float(polygon[i][1] - polygon[i - 1][1]) / (polygon[i][0] - polygon[i - 1][0])])
            else:
                polygon.insert(i, [polygon[i - 1][0] + (position - polygon[i - 1][1]) * float(polygon[i][0] - polygon[i - 1][0]) / (polygon[i][1] - polygon[i - 1][1]), position])
            cross.insert(i, 0)
            bnd[1 * (cross[i + 1] > 0)].append(i)
            i += 2
        elif cross[i] == 0:
            j = i + 1
            while cross[j] == 0:
                j += 1
            if cross[i - 1] * cross[j] < 0:
                bnd[1 * (cross[j] > 0)].append(j - 1)
            i = j + 1
        else:
            i += 1
    if len(bnd[0]) == 0:
        out_polygons[1 * (numpy.sum(cross) > 0)].append(polygon)
        return out_polygons
    bnd = (numpy.array(bnd[0]), numpy.array(bnd[1]))
    bnd = (list(bnd[0][numpy.argsort(numpy.array(polygon)[bnd[0], 1 - axis])]),
           list(bnd[1][numpy.argsort(numpy.array(polygon)[bnd[1], 1 - axis])]))
    cross = numpy.ones(len(polygon), dtype=int)
    cross[bnd[0]] = -2
    cross[bnd[1]] = -1
    i = 0
    while i < len(polygon):
        if cross[i] > 0 and polygon[i][axis] != position:
            start = i
            side = 1 * (polygon[i][axis] > position)
            out_polygons[side].append([polygon[i]])
            cross[i] = 0
            nxt = i + 1
            if nxt == len(polygon):
                nxt = 0
            boundary = True
            while nxt != start:
                out_polygons[side][-1].append(polygon[nxt])
                if cross[nxt] > 0:
                    cross[nxt] = 0
                if cross[nxt] < 0 and boundary:
                    j = bnd[cross[nxt] + 2].index(nxt)
                    nxt = bnd[-cross[nxt] - 1][j]
                    boundary = False
                else:
                    nxt += 1
                    if nxt == len(polygon):
                        nxt = 0
                    boundary = True
        i += 1
    return out_polygons


def slice(objects, position, axis, layer=0, datatype=0):
    """
    Slice polygons and polygon sets at given positions along an axis.

    Parameters
    ----------
    objects : ``Polygon``, ``PolygonSet``, or list
        Operand of the slice operation.  If this is a list, each element
        must be a ``Polygon``, ``PolygonSet``, ``CellReference``,
        ``CellArray``, or an array-like[N][2] of vertices of a polygon.
    position : number or list of numbers
        Positions to perform the slicing operation along the specified axis.
    axis : 0 or 1
        Axis along which the polygon will be sliced.
    layer : integer, list
        The GDSII layer numbers for the elements between each division. If
        the number of layers in the list is less than the number of divided
        regions, the list is repeated.
    datatype : integer
        The GDSII datatype for the resulting element (between 0 and 255).

    Returns
    -------
    out : list[N] of PolygonSet
        Result of the slicing operation, with N = len(positions) + 1. Each
        PolygonSet comprises all polygons between 2 adjacent slicing
        positions, in crescent order.

    Examples
    --------
    >>> ring = gdspy.Round((0, 0), 10, inner_radius = 5)
    >>> result = gdspy.slice(ring, [-7, 7], 0)
    >>> cell.add(result[1])
    """
    if (layer.__class__ != [].__class__):
        layer = [layer]
    if (objects.__class__ != [].__class__):
        objects = [objects]
    if (position.__class__ != [].__class__):
        position = [position]
    position.sort()
    result = [[] for i in range(len(position) + 1)]
    polygons = []
    for obj in objects:
        if isinstance(obj, Polygon):
            polygons.append(obj.points)
        elif isinstance(obj, PolygonSet):
            polygons += obj.polygons
        elif isinstance(obj, CellReference) or isinstance(obj, CellArray):
            polygons += obj.get_polygons()
        else:
            polygons.append(obj)
    for i, p in enumerate(position):
        nxt_polygons = []
        for pol in polygons:
            (pol1, pol2) = _chop(pol, p, axis)
            result[i] += pol1
            nxt_polygons += pol2
        polygons = nxt_polygons
    result[-1] = polygons
    for i in range(len(result)):
        result[i] = PolygonSet(result[i], layer[i % len(layer)], datatype)
    return result


def offset(polygons, distance, join='miter', tolerance=2, precision=0.001, max_points=199, layer=0, datatype=0):
    """
    Shrink or expand a polygon or polygon set.

    Parameters
    ----------
    polygons : polygon or array-like
        Polygons to be offset. Must be a ``Polygon``, ``PolygonSet``,
        ``CellReference``, ``CellArray``, or an array. The array may
        contain any of the previous objects or an array-like[N][2] of
        vertices of a polygon.
    distance : number
        Offset distance. Positive to expand, negative to shrink.
    join : {'miter', 'bevel', 'round'}
        Type of join used to create the offset polygon.
    tolerance : integer or float
        For miter joints, this number must be at least 2 and it represents
        the maximun distance in multiples of offset betwen new vertices and
        their original position before beveling to avoid spikes at acute
        joints. For round joints, it indicates the curvature resolution in
        number of points per full circle.
    precision : float
        Desired precision for rounding vertice coordinates.
    max_points : integer
        If greater than 4, fracture the resulting polygons to ensure they
        have at most ``max_points`` vertices. This is not a tessellating
        function, so this number should be as high as possible. For example,
        it should be set to 199 for polygons being drawn in GDSII files.
    layer : integer
        The GDSII layer number for the resulting element.
    datatype : integer
        The GDSII datatype for the resulting element (between 0 and 255).

    Returns
    -------
    out : ``PolygonSet`` or ``None``
        Return the offset shape as a set of polygons.
    """
    poly = []
    if isinstance(polygons, Polygon):
        poly.append(polygons.points)
    elif isinstance(polygons, PolygonSet):
        poly += polygons.polygons
    elif isinstance(polygons, CellReference) or isinstance(polygons, CellArray):
        poly += polygons.get_polygons()
    else:
        for obj in polygons:
            if isinstance(obj, Polygon):
                poly.append(obj.points)
            elif isinstance(obj, PolygonSet):
                poly += obj.polygons
            elif isinstance(obj, CellReference) or isinstance(obj, CellArray):
                poly += obj.get_polygons()
            else:
                poly.append(obj)
    result = clipper.offset(poly, distance, join, tolerance, 1/precision)
    return None if result is None else PolygonSet(result, layer, datatype, False).fracture(max_points)



def boolean(polygons, operation, max_points=199, layer=0, datatype=0, eps=1e-13):
    """
    This function is deprecated in favor of 'fast_boolean'.

    Execute any generalized boolean operation on polygons and polygon sets.

    Parameters
    ----------
    polygons : array-like
        Operands of the boolean operation. Each element of this array must
        be a ``Polygon``, ``PolygonSet``, ``CellReference``, ``CellArray``,
        or an array-like[N][2] of vertices of a polygon.
    operation : function
        Function that accepts as input ``len(polygons)`` integers. Each
        integer represents the incidence of the corresponding ``polygon``.
        The function must return a bool or integer (interpreted as bool).
    max_points : integer
        If greater than 4, fracture the resulting polygons to ensure they
        have at most ``max_points`` vertices. This is not a tessellating
        function, so this number should be as high as possible. For example,
        it should be set to 199 for polygons being drawn in GDSII files.
    layer : integer
        The GDSII layer number for the resulting element.
    datatype : integer
        The GDSII datatype for the resulting element (between 0 and 255).
    eps : positive number
        Small number to be used as tolerance in intersection and overlap
        calculations.

    Returns
    -------
    out : PolygonSet or ``None``
        Result of the boolean operation.

    Notes
    -----
    Since ``operation`` receives a list of integers as input, it can be
    somewhat more general than boolean operations only. See the examples
    below.

    Because of roundoff errors there are a few cases when this function can
    cause segmentation faults. If that happens, increasing the value of
    ``eps`` might help.

    Examples
    --------
    >>> circle = gdspy.Round((0, 0), 10)
    >>> triangle = gdspy.Round((0, 0), 12, number_of_points=3)
    >>> bad_poly = gdspy.L1Path((0, 0), '+y', 2,
            [6, 4, 4, 8, 4, 5, 10], [-1, -1, -1, 1, 1, 1])
    >>> union = gdspy.boolean([circle, triangle],
            lambda cir, tri: cir or tri)
    >>> intersection = gdspy.boolean([circle, triangle],
            lambda cir, tri: cir and tri)
    >>> subtraction = gdspy.boolean([circle, triangle],
            lambda cir, tri: cir and not tri)
    >>> multi_xor = gdspy.boolean([badPath], lambda p: p % 2)
    """
    warnings.warn("[GDSPY] Function 'boolean' is deprecated and will be removed in a future version. Please use 'fast_boolean' instead.", stacklevel=2)
    poly = []
    indices = [0]
    special_function = False
    for obj in polygons:
        if isinstance(obj, Polygon):
            poly.append(obj.points)
            indices.append(indices[-1] + 1)
        elif isinstance(obj, PolygonSet):
            special_function = True
            poly += obj.polygons
            indices.append(indices[-1] + len(obj.polygons))
        elif isinstance(obj, CellReference) or isinstance(obj, CellArray):
            special_function = True
            a = obj.get_polygons()
            poly += a
            indices.append(indices[-1] + len(a))
        else:
            poly.append(obj)
            indices.append(indices[-1] + 1)
    if special_function:
        result = boolext.clip(poly, lambda *p: operation(*[sum(p[indices[ia]:indices[ia + 1]]) for ia in range(len(indices) - 1)]), eps)
    else:
        result = boolext.clip(poly, operation, eps)
    return None if result is None else PolygonSet(result, layer, datatype, False).fracture(max_points)


def fast_boolean(operandA, operandB, operation, precision=0.001, max_points=199, layer=0, datatype=0):
    """
    Execute any boolean operation between 2 polygons or polygon sets.

    Parameters
    ----------
    operandA : polygon or array-like
        First operand. Must be a ``Polygon``, ``PolygonSet``,
        ``CellReference``, ``CellArray``, or an array. The array may
        contain any of the previous objects or an array-like[N][2] of
        vertices of a polygon.
    operandB : polygon, array-like or ``None``
        Second operand. Must be ``None``, a ``Polygon``, ``PolygonSet``,
        ``CellReference``, ``CellArray``, or an array. The array may
        contain any of the previous objects or an array-like[N][2] of
        vertices of a polygon.
    operation : {'or', 'and', 'xor', 'not'}
        Boolean operation to be executed. The 'not' operation returns
        the difference ``operandA - operandB``.
    precision : float
        Desired precision for rounding vertice coordinates.
    max_points : integer
        If greater than 4, fracture the resulting polygons to ensure they
        have at most ``max_points`` vertices. This is not a tessellating
        function, so this number should be as high as possible. For example,
        it should be set to 199 for polygons being drawn in GDSII files.
    layer : integer
        The GDSII layer number for the resulting element.
    datatype : integer
        The GDSII datatype for the resulting element (between 0 and 255).

    Returns
    -------
    out : PolygonSet or ``None``
        Result of the boolean operation.
    """
    polyA = []
    polyB = []
    for poly,obj in zip((polyA, polyB),(operandA, operandB)):
        if isinstance(obj, Polygon):
            poly.append(obj.points)
        elif isinstance(obj, PolygonSet):
            poly += obj.polygons
        elif isinstance(obj, CellReference) or isinstance(obj, CellArray):
            poly += obj.get_polygons()
        elif obj is not  None:
            for inobj in obj:
                if isinstance(inobj, Polygon):
                    poly.append(inobj.points)
                elif isinstance(inobj, PolygonSet):
                    poly += inobj.polygons
                elif isinstance(inobj, CellReference) or isinstance(inobj, CellArray):
                    poly += inobj.get_polygons()
                else:
                    poly.append(inobj)
    if len(polyB) == 0:
        polyB.append(polyA.pop())
    result = clipper.clip(polyA, polyB, operation, 1/precision)
    return None if result is None else PolygonSet(result, layer, datatype, False).fracture(max_points)


def gds_print(outfile, cells=None, name='library', unit=1.0e-6, precision=1.0e-9):
    """
    Output a list of cells as a GDSII stream library.

    The dimensions actually written on the GDSII file will be the dimensions
    of the objects created times the ratio ``unit/precision``. For example,
    if a circle with radius 1.5 is created and we set ``unit=1.0e-6`` (1 um)
    and ``precision=1.0e-9`` (1 nm), the radius of the circle will be 1.5 um
    and the GDSII file will contain the dimension 1500 nm.

    Parameters
    ----------
    outfile : file or string
        The file (or path) where the GDSII stream will be written. It must
        be opened for writing operations in binary format.
    cells : array-like
        The list of cells or cell names to be included in the library. If
        ``None``, all cells listed in ``Cell.cell_dict`` are used.
    name : string
        Name of the GDSII library (file).
    unit : number
        Unit size for the objects in the library (in *meters*).
    precision : number
        Precision for the dimensions of the objects in the library (in
        *meters*).

    Examples
    --------
    >>> gdspy.gds_print('out-file.gds', unit=1.0e-6, precision=1.0e-9)
    """
    if outfile.__class__ == ''.__class__:
        outfile = open(outfile, 'wb')
        close = True
    else:
        close = False
    if cells == None:
        cells = Cell.cell_dict.itervalues()
    else:
        cells = [Cell.cell_dict.get(c, c) for c in cells]
        i = 0
        while i < len(cells):
            for cell in cells[i].get_dependencies():
                if cell not in cells:
                    cells.append(cell)
            i += 1
    now = datetime.datetime.today()
    if len(name)%2 != 0:
        name = name + '\0'
    outfile.write(struct.pack('>19h', 6, 0x0002, 0x0258, 28, 0x0102, now.year, now.month, now.day, now.hour, now.minute, now.second, now.year, now.month, now.day, now.hour, now.minute, now.second, 4+len(name), 0x0206) + name.encode('ascii') + struct.pack('>2h', 20, 0x0305) + _eight_byte_real(precision / unit) + _eight_byte_real(precision))
    for cell in cells:
        outfile.write(cell.to_gds(unit / precision))
    outfile.write(struct.pack('>2h', 4, 0x0400))
    if close:
        outfile.close()


class GdsPrint:
    """
    GDSII strem library printer.

    The dimensions actually written on the GDSII file will be the dimensions
    of the objects created times the ratio ``unit/precision``. For example,
    if a circle with radius 1.5 is created and we set ``unit=1.0e-6`` (1 um)
    and ``precision=1.0e-9`` (1 nm), the radius of the circle will be 1.5 um
    and the GDSII file will contain the dimension 1500 nm.

    Parameters
    ----------
    outfile : file or string
        The file (or path) where the GDSII stream will be written. It must
        be opened for writing operations in binary format.
    cells : array-like
        The list of cells or cell names to be included in the library. If
        ``None``, all cells listed in ``Cell.cell_dict`` are used.
    name : string
        Name of the GDSII library (file).
    unit : number
        Unit size for the objects in the library (in *meters*).
    precision : number
        Precision for the dimensions of the objects in the library (in
        *meters*).

    Notes
    -----

    This class can be used for incremental output of the geometry in case
    the complete layout is too large to be kept in memory all at once.

    Examples
    --------
    >>> prt = gdspy.GdsPrint('out-file.gds', unit=1.0e-6, precision=1.0e-9)
    >>> for i in range(10):
    ...     cell = gdspy.Cell('C{}'.format(i))
    ...     # Add the contents of this cell...
    ...     prt.write_cell(cell)
    ...     # Clear the memory: erase Cell objects and any other objects
    ...     # that won't be needed.
    ...     gdspy.Cell.cell_dict.clear()
    >>> prt.close()
    """

    def __init__(self, outfile, name='library', unit=1.0e-6, precision=1.0e-9):
        if outfile.__class__ == ''.__class__:
            self._outfile = open(outfile, 'wb')
            self._close = True
        else:
            self._outfile = outfile
            self._close = False
        self._res = unit / precision
        now = datetime.datetime.today()
        if len(name)%2 != 0:
            name = name + '\0'
        self._outfile.write(struct.pack('>19h', 6, 0x0002, 0x0258, 28, 0x0102, now.year, now.month, now.day, now.hour, now.minute, now.second, now.year, now.month, now.day, now.hour, now.minute, now.second, 4+len(name), 0x0206) + name.encode('ascii') + struct.pack('>2h', 20, 0x0305) + _eight_byte_real(precision / unit) + _eight_byte_real(precision))

    def write_cell(self, cell):
        """
        Write the specified cell to the file.

        Parameters
        ----------
        cell : ``Cell``
            Cell to be written.

        Notes
        -----
        Only the specified cell is written. Unlike in ``gds_print``, cell
        dependencies are not automatically included.

        Returns
        -------
        out : ``GdsPrint``
            This object.
        """
        self._outfile.write(cell.to_gds(self._res))
        return self

    def close(self):
        """
        Finalize the GDSII stream library.
        """
        self._outfile.write(struct.pack('>2h', 4, 0x0400))
        if self._close:
            self._outfile.close()
