######################################################################
#                                                                    #
#  Copyright 2009-2018 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################
"""
gdspy is a Python module that allows the creation of GDSII stream files.

Many features of the GDSII format are implemented, such as cell
references and arrays, but the support for fonts is quite limited.

GDSII format references:
- http://boolean.klaasholwerda.nl/interface/bnf/gdsformat.html
- http://www.artwork.com/gdsii/gdsii/
- http://www.buchanan1.net/stream_description.html
"""

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
import datetime
import warnings
import itertools
import numpy
import copy as libCopy
import hashlib

from gdspy import clipper
try:
    from gdspy.viewer import LayoutViewer
except ImportError as e:
    warnings.warn("[GDSPY] LayoutViewer not available: " + str(e), category=ImportWarning, stacklevel=2)

__version__ = '1.3.1'

_halfpi = 0.5 * numpy.pi
_directions_dict = {'+x': 0, '+y': 0.5, '-x': 1, '-y': -0.5}
_directions_list = ['+x', '+y', '-x', '-y']

_bounding_boxes = {}


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
    if value == 0:
        return b'\x00\x00\x00\x00\x00\x00\x00\x00'
    if value < 0:
        byte1 = 0x80
        value = -value
    else:
        byte1 = 0x00
    fexp = numpy.log2(value) / 4
    exponent = int(numpy.ceil(fexp))
    if fexp == exponent:
        exponent += 1
    mantissa = int(value * 16.**(14 - exponent))
    byte1 += exponent + 64
    byte2 = (mantissa // 281474976710656)
    short3 = (mantissa % 281474976710656) // 4294967296
    long4 = mantissa % 4294967296
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
    exponent = (short1 & 0x7f00) // 256 - 64
    mantissa = (((short1 & 0x00ff) * 65536 + short2) * 4294967296 + long3) / 72057594037927936.0
    if short1 & 0x8000:
        return -mantissa * 16.**exponent
    return mantissa * 16.**exponent


def _hobby(points, angles=None, curl_start=1, curl_end=1, t_in=1, t_out=1, cycle=False):
    """
    Calculate control points for a cubic Bezier curve interpolating the
    given points.

    Reference: Hobby, J.D.  *Discrete Comput. Geom.* (1986) 1: 123.
    https://doi.org/10.1007/BF02187690

    Parameters
    ----------
    points : Numpy array[N, 2]
        Vertices in the interpolating curve.
    angles : array-like[N] or ``None``
        Tangent angles at each point (in *radians*).  Any angles defined
        as ``None`` are automatically calculated.
    curl_start : number
        Ratio between the mock curvatures at the first point and at its
        neighbor.  A value of 1 renders the first segment a good
        approximation for a circular arc.  A value of 0 will better
        approximate a straight segment.  It has no effect for closed
        curves or when an angle is defined for the first point.
    curl_end : number
        Ratio between the mock curvatures at the last point and at its
        neighbor.  It has no effect for closed curves or when an angle
        is defined for the last point.
    t_in : number or array-like[N]
        Tension parameter when arriving at each point.  One value per
        point or a single value used for all points.
    t_out : number or array-like[N]
        Tension parameter when leaving each point.  One value per point
        or a single value used for all points.
    cycle : bool
        If ``True``, calculates control points for a closed curve, with
        an additional segment connecting the first and last points.

    Returns
    -------
    out : 2-tuple of Numpy array[M, 2]
        Pair of control points for each segment in the interpolating
        curve.  For a closed curve (``cycle = True``), M = N.  For an
        open curve (``cycle = False``), M = N - 1.
    """
    z = points[:, 0] + 1j * points[:, 1]
    if not hasattr(t_in, '__iter__'):
        t_in = t_in * numpy.ones(len(z))
    if not hasattr(t_out, '__iter__'):
        t_out = t_out * numpy.ones(len(z))
    if angles is None:
        angles = [None for _ in range(len(z))]
    v = numpy.roll(z, -1) - z
    d = numpy.abs(v)
    delta = numpy.angle(v)
    psi = (delta - numpy.roll(delta, 1) + numpy.pi) % (2 * numpy.pi) - numpy.pi
    theta = numpy.zeros(len(z))
    phi = numpy.zeros(len(z))
    for i in range(len(angles)):
        if angles[i] is not None:
            theta[i] = angles[i] - delta[i]
            phi[i] = delta[i - 1] - angles[i]
    i = 0
    while i < len(angles):
        while i < len(angles) and angles[i] is not None:
            i += 1
        j = i
        while j < len(angles) and angles[j] is None:
            j += 1
        if j > i:
            n = 2 * (j - i)
            m = numpy.zeros((n, n))
            coeff = numpy.zeros(n)
            coeff[::2] = -psi[i:j]
            if i == 0:
                if cycle:
                    m[1, n - 2] = d[0] * t_out[0]**2 * t_out[1]
                    m[1, 0] = -(1 - 3 * t_out[1]) * d[-1] * t_in[0]**2 * t_in[-1]
                    m[1, 1] = (1 - 3 * t_in[-1]) * d[0] * t_out[0]**2 * t_out[1]
                    m[1, 3] = -d[-1] * t_in[0]**2 * t_in[-1]
                else:
                    m[1, 0] = t_in[0]**3 * (1 - 3 * t_out[1]) - curl_start * t_out[1]**3
                    m[1, 3] = t_in[0]**3 - curl_start * t_out[1]**3 * (1 - 3 * t_in[0])
            else:
                m[1, 0] = -(1 - 3 * t_out[i + 1]) * d[i - 1] * t_in[i]**2 * t_in[i - 1]
                m[1, 1] = (1 - 3 * t_in[i - 1]) * d[i] * t_out[i]**2 * t_out[i + 1]
                m[1, 3] = -d[i - 1] * t_in[i]**2 * t_in[i - 1]
                coeff[1] = -d[i] * t_out[i]**2 * t_out[i + 1] * theta[i - 1]
            if j == len(angles):
                if cycle:
                    m[n - 1, n - 4] = d[j - 1] * t_out[j - 1]**2 * t_out[0]
                    m[n - 1, n - 2] = -(1 - 3 * t_out[0]) * d[j - 2] * t_in[j - 1]**2 * t_in[j - 2]
                    m[n - 1, n - 1] = (1 - 3 * t_in[j - 2]) * d[j - 1] * t_out[j - 1]**2 * t_out[0]
                    m[n - 1, 1] = -d[j - 2] * t_in[j - 1]**2 * t_in[j - 2]
                else:
                    m[n - 1, n - 4] = t_out[j - 1]**3 - curl_end * t_in[j - 2]**3 * (1 - 3 * t_out[j - 1])
                    m[n - 1, n - 1] = t_out[j - 1]**3 * (1 - 3 * t_in[j - 2]) - curl_end * t_in[j - 2]**3
            else:
                m[n - 1, n - 4] = d[j - 1] * t_out[j - 1]**2 * t_out[j]
                m[n - 1, n - 2] = -(1 - 3 * t_out[j]) * d[j - 2] * t_in[j - 1]**2 * t_in[j - 2]
                m[n - 1, n - 1] = (1 - 3 * t_in[j - 2]) * d[j - 1] * t_out[j - 1]**2 * t_out[j]
                coeff[n - 1] = d[j - 2] * t_in[j - 1]**2 * t_in[j - 2] * phi[j]
            m[0, 0] = 1
            m[0, 1] = 1
            m[n - 2, n - 2] = 1
            m[n - 2, n - 1] = 1
            for k in range(i + 1, j - 1):
                row = 2 * (k - i)
                m[row, row] = 1
                m[row, row + 1] = 1
                m[row + 1, row - 2] = d[k] * t_out[k]**2 * t_out[k + 1]
                m[row + 1, row] = -(1 - 3 * t_out[k + 1]) * d[k - 1] * t_in[k]**2 * t_in[k - 1]
                m[row + 1, row + 1] = (1 - 3 * t_in[k - 1]) * d[k] * t_out[k]**2 * t_out[k + 1]
                m[row + 1, row + 3] = -d[k - 1] * t_in[k]**2 * t_in[k - 1]
            sol = numpy.linalg.solve(m, coeff)
            theta[i:j] = sol[::2]
            phi[i:j] = sol[1::2]
            i = j
    w = numpy.exp(1j * (theta + delta))
    a = numpy.sqrt(2)
    b = 1.0 / 16
    c = (3 - numpy.sqrt(5)) / 2
    sintheta = numpy.sin(theta)
    costheta = numpy.cos(theta)
    sinphi = numpy.sin(numpy.roll(phi, -1))
    cosphi = numpy.cos(numpy.roll(phi, -1))
    alpha = a * (sintheta - b * sinphi) * (sinphi - b * sintheta) * (costheta - cosphi)
    cta = z + w * d * ((2 + alpha) / (1 + (1 - c) * costheta + c * cosphi)) / (3 * t_in)
    ctb = numpy.roll(z, -1) - numpy.roll(w, -1) * d * (
        (2 - alpha) / (1 + (1 - c) * cosphi + c * costheta)) / (3 * numpy.roll(t_out, -1))
    if not cycle:
        cta = cta[:-1]
        ctb = ctb[:-1]
    return (numpy.vstack((cta.real, cta.imag)).transpose(), numpy.vstack((ctb.real, ctb.imag)).transpose())


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

    Notes
    -----
    The last point should not be equal to the first (polygons are
    automatically closed).

    The GDSII specification supports only a maximum of 199 vertices per
    polygon.
    """

    __slots__ = 'layers', 'datatypes', 'polygons'

    def __init__(self, polygons, layer=0, datatype=0):
        self.polygons = [numpy.array(p) for p in polygons]
        self.layers = [layer] * len(self.polygons)
        self.datatypes = [datatype] * len(self.polygons)
        if poly_warnings and any(p.shape[0] > 199 for p in self.polygons):
            warnings.warn("[GDSPY] A polygon with more than 199 points was created (not officially supported by the GDSII format).", stacklevel=2)

    def __str__(self):
        return ("PolygonSet ({} polygons, {} vertices, layers {}, datatypes {})").format(len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))

    def get_bounding_box(self):
        """
        Returns the bounding box of the polygons.

        Returns
        -------
        out : Numpy array[2, 2] or ``None``
            Bounding box of this polygon in the form [[x_min, y_min],
            [x_max, y_max]], or ``None`` if the polygon is empty.
        """
        if len(self.polygons) == 0:
            return None
        return numpy.array(((min(pts[:, 0].min() for pts in self.polygons),
                             min(pts[:, 1].min() for pts in self.polygons)),
                            (max(pts[:, 0].max() for pts in self.polygons),
                             max(pts[:, 1].max() for pts in self.polygons))))

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
        self.polygons = [(points - c0) * ca + (points - c0)[:, ::-1] * sa + c0 for points in self.polygons]
        return self

    def scale(self, scalex, scaley=None, center=(0, 0)):
        """
        Scale this object.

        Parameters
        ----------
        scalex : number
            Scaling factor along the first axis.
        scaley : number or ``None``
            Scaling factor along the second axis.  If ``None``, same as
            ``scalex``.
        center : array-like[2]
            Center point for the scaling operation.

        Returns
        -------
        out : ``PolygonSet``
            This object.
        """
        c0 = numpy.array(center)
        s = scalex if scaley is None else numpy.array((scalex, scaley))
        self.polygons = [(points - c0) * s + c0 for points in self.polygons]
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
            if poly_warnings and len(self.polygons[ii]) > 8190:
                warnings.warn("[GDSPY] Polygons with more than 8190 are not supported by the official GDSII specification.  This extension might not be compatible with all GDSII readers.", stacklevel=3)
                data.append(struct.pack('>4Hh2Hh', 4, 0x0800, 6, 0x0D02, self.layers[ii], 6, 0x0E02,
                                        self.datatypes[ii]))
                xy = numpy.empty((self.polygons[ii].shape[0] + 1, 2), dtype='>i4')
                xy[:-1, :] = numpy.round(self.polygons[ii] * multiplier)
                xy[-1, :] = xy[0, :]
                i0 = 0
                while i0 < xy.shape[0]:
                    i1 = min(i0 + 8191, xy.shape[0])
                    data.append(struct.pack('>2H', 4 + 8 * (i1 - i0), 0x1003))
                    data.append(xy[i0:i1].tostring())
                    i0 = i1
                data.append(struct.pack('>2H', 4, 0x1100))
            else:
                data.append(struct.pack('>4Hh2Hh2H', 4, 0x0800, 6, 0x0D02, self.layers[ii], 6, 0x0E02,
                                        self.datatypes[ii], 12 + 8 * len(self.polygons[ii]), 0x1003))
                xy = numpy.round(self.polygons[ii] * multiplier).astype('>i4')
                data.append(xy.tostring())
                data.append(xy[0].tostring())
                data.append(struct.pack('>2H', 4, 0x1100))
        return b''.join(data)

    def area(self, by_spec=False):
        """
        Calculate the total area of the path(s).

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with
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
                    poly_area += (poly[0][0] - poly[ii + 1][0]) * (poly[ii][1] - poly[0][1]) - (poly[0][1] - poly[ii + 1][1]) * (poly[ii][0] - poly[0][0])
                if key in path_area:
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

    def fracture(self, max_points=199, precision=1e-3):
        """
        Slice these polygons in the horizontal and vertical directions
        so that each resulting piece has at most ``max_points``.  This
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
        out : ``PolygonSet``
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
                        cuts = [pts0[int(i * len(pts0) / (ncuts + 1.0) + 0.5)] for i in range(1, ncuts + 1)]
                        chopped = clipper._chop(self.polygons[ii], cuts, 0, 1 / precision)
                    else:
                        # Horizontal cuts
                        cuts = [pts1[int(i * len(pts1) / (ncuts + 1.0) + 0.5)] for i in range(1, ncuts + 1)]
                        chopped = clipper._chop(self.polygons[ii], cuts, 1, 1 / precision)
                    self.polygons.pop(ii)
                    layer = self.layers.pop(ii)
                    datatype = self.datatypes.pop(ii)
                    self.polygons.extend(numpy.array(x) for x in itertools.chain.from_iterable(chopped))
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
        radius : number, list
            Radius of the corners.  If number: all corners filleted by
            that amount.  If list: specify fillet radii on a per-polygon
            basis (list length must be equal to the number of polygons
            in this ``PolygonSet``).  Each element in the list can be a
            number (all corners filleted by the same amount) or a list
            of numbers, one per polygon vertex.  Alternatively, the list
            can be flattened to have one radius per ``PolygonSet``
            vertex, instead of being a list of lists.
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
        out : ``PolygonSet``
            This object.
        """
        two_pi = 2 * numpy.pi
        fracture = False
        if hasattr(radius, '__iter__'):
            if len(radius) == len(self.polygons):
                radii = []
                for r, p in zip(radius, self.polygons):
                    if hasattr(r, '__iter__'):
                        if len(r) != p.shape[0]:
                            raise ValueError("[GDSPY] Wrong length in fillet radius list.  Expected lengths are {} or {}; got {}.".format(len(self.polygons), total, len(radius)))
                        radii.append(r)
                    else:
                        radii.append([r] * p.shape[0])
            else:
                total = sum(p.shape[0] for p in self.polygons)
                if len(radius) != total:
                    raise ValueError("[GDSPY] Wrong length in fillet radius list.  Expected lengths are {} or {}; got {}.".format(len(self.polygons), total, len(radius)))
                radii = []
                n = 0
                for p in self.polygons:
                    radii.append(radius[n:n + p.shape[0]])
                    n += p.shape[0]
        else:
            radii = [[radius] * p.shape[0] for p in self.polygons]

        for jj in range(len(self.polygons)):
            vec = self.polygons[jj].astype(float) - numpy.roll(self.polygons[jj], 1, 0)
            length = numpy.sqrt(numpy.sum(vec**2, 1))
            ii = numpy.flatnonzero(length)
            if len(ii) < len(length):
                self.polygons[jj] = self.polygons[jj][ii]
                vec = self.polygons[jj] - numpy.roll(self.polygons[jj], 1, 0)
                length = numpy.sqrt(numpy.sum(vec**2, 1))
            vec[:, 0] = vec[:, 0] / length
            vec[:, 1] = vec[:, 1] / length
            dvec = numpy.roll(vec, -1, 0) - vec
            norm = numpy.sqrt(numpy.sum(dvec**2, 1))
            ii = numpy.flatnonzero(norm)
            dvec[ii, 0] = dvec[ii, 0] / norm[ii]
            dvec[ii, 1] = dvec[ii, 1] / norm[ii]
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
                    ll = radii[jj][ii] * tt[ii]
                    if ll > 0.49 * length[ii]:
                        r = 0.49 * length[ii] / tt[ii]
                        ll = 0.49 * length[ii]
                    else:
                        r = radii[jj][ii]
                    if ll > 0.49 * length[ii + 1]:
                        r = 0.49 * length[ii + 1] / tt[ii]
                    new_points.extend(r * dvec[ii] / ct[ii] + self.polygons[jj][ii] +
                                      numpy.vstack((r * numpy.cos(a), r * numpy.sin(a))).transpose())
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
        Move the polygons from one place to another

        Parameters
        ----------
        dx : number
            distance to move in the x-direction
        dy : number
            distance to move in the y-direction

        Returns
        -------
        out : ``PolygonSet``
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
        out : ``PolygonSet``
            This object.
        """
        origin = numpy.array(p1)
        vec = numpy.array(p2) - origin
        vec_r = vec * (2 / numpy.inner(vec, vec))
        self.polygons = [numpy.outer(numpy.inner(points - origin, vec_r), vec) - points + 2 * origin
                         for points in self.polygons]
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

    The GDSII specification supports only a maximum of 199 vertices per
    polygon.

    Examples
    --------
    >>> triangle_pts = [(0, 40), (15, 40), (10, 50)]
    >>> triangle = gdspy.Polygon(triangle_pts)
    >>> myCell.add(triangle)
    """

    __slots__ = 'layers', 'datatypes', 'polygons'

    def __init__(self, points, layer=0, datatype=0):
        if poly_warnings and len(points) > 199:
            warnings.warn("[GDSPY] A polygon with more than 199 points was created (not officially supported by the GDSII format).", stacklevel=2)
        self.layers = [layer]
        self.datatypes = [datatype]
        self.polygons = [numpy.array(points)]

    def __str__(self):
        return "Polygon ({} vertices, layer {}, datatype {})".format(len(self.polygons[0]), self.layers[0], self.datatypes[0])


class Rectangle(PolygonSet):
    """
    Rectangular geometric object.

    Parameters
    ----------
    point1 : array-like[2]
        Coordinates of a corner of the rectangle.
    point2 : array-like[2]
        Coordinates of the corner of the rectangle opposite to
        ``point1``.
    layer : integer
        The GDSII layer number for this element.
    datatype : integer
        The GDSII datatype for this element (between 0 and 255).

    Examples
    --------
    >>> rectangle = gdspy.Rectangle((0, 0), (10, 20))
    >>> myCell.add(rectangle)
    """

    __slots__ = 'layers', 'datatypes', 'polygons'

    def __init__(self, point1, point2, layer=0, datatype=0):
        self.layers = [layer]
        self.datatypes = [datatype]
        self.polygons = [numpy.array([[point1[0], point1[1]], [point1[0], point2[1]],
                                      [point2[0], point2[1]], [point2[0], point1[1]]])]

    def __str__(self):
        return ("Rectangle (({0[0]}, {0[1]}) to ({1[0]}, {1[1]}), layer {2}, datatype {3})").format(self.polygons[0][0], self.polygons[0][2], self.layers[0], self.datatypes[0])

    def __repr__(self):
        return "Rectangle(({0[0]}, {0[1]}), ({1[0]}, {1[1]}), {2}, {3})".format(self.polygons[0][0], self.polygons[0][2], self.layers[0], self.datatypes[0])


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
    number_of_points : integer or ``None``
        Manually define the number of vertices that form the object
        (polygonal approximation).  Overrides ``tolerance``.
    max_points : integer
        If the number of points in the element is greater than
        ``max_points``, it will be fractured in smaller polygons with
        at most ``max_points`` each.  If ``max_points = 0`` no fracture
        will occur.
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
    >>> ell_ring = gdspy.Round((50, 5), (8, 7), inner_radius=(5, 4))
    >>> pie_slice = gdspy.Round((30, 25), 8, initial_angle=0,
    ...                             final_angle=-5.0*numpy.pi/6.0)
    >>> arc = gdspy.Round((50, 25), 8, inner_radius=5,
    ...                       initial_angle=-5.0*numpy.pi/6.0,
    ...                       final_angle=0)
    """

    __slots__ = 'layers', 'datatypes', 'polygons'

    def __init__(self, center, radius, inner_radius=0, initial_angle=0, final_angle=0,
                 tolerance=0.01, number_of_points=None, max_points=199, layer=0, datatype=0):
        if hasattr(radius, '__iter__'):
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

        if hasattr(inner_radius, '__iter__'):
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
            warnings.warn("[GDSPY] Use of a floating number as number_of_points is deprecated in favor of tolerance.", category=DeprecationWarning, stacklevel=2)
            tolerance = number_of_points
            number_of_points = None

        if number_of_points is None:
            if inner_radius <= 0:
                if final_angle == initial_angle:
                    number_of_points = int(2 * radius * numpy.pi / tolerance + 0.5)
                else:
                    number_of_points = int(abs(final_angle - initial_angle) * radius / tolerance + 0.5) + 2
            else:
                if final_angle == initial_angle:
                    number_of_points = 2 * int(2 * radius * numpy.pi / tolerance + 0.5) + 2
                else:
                    number_of_points = 2 * int(abs(final_angle - initial_angle) * radius / tolerance + 0.5) + 2

        number_of_points = max(number_of_points, 3)
        pieces = 1 if max_points == 0 else int(numpy.ceil(number_of_points / float(max_points)))
        number_of_points = number_of_points // pieces
        self.layers = [layer] * pieces
        self.datatypes = [datatype] * pieces
        self.polygons = [numpy.zeros((number_of_points, 2)) for _ in range(pieces)]
        if final_angle == initial_angle and pieces > 1:
            final_angle += 2 * numpy.pi
        angles = numpy.linspace(initial_angle, final_angle, pieces + 1)
        oang = outer_transform(angles)
        iang = inner_transform(angles)
        for ii in range(pieces):
            if oang[ii + 1] == oang[ii]:
                if inner_radius <= 0:
                    t = numpy.arange(number_of_points) * 2.0 * numpy.pi / number_of_points
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
                    self.polygons[ii][0, :] += center
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
        return ("Round ({} polygons, {} vertices, layers {}, datatypes {})").format(len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))


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
        If ``True``, the text is written from left to right; if
        ``False``, from top to bottom.
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

    __slots__ = 'layers', 'datatypes', 'polygons'

    def __init__(self, text, size, position=(0, 0), horizontal=True, angle=0, layer=0, datatype=0):
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
                if text[jj] in Text._font:
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
        return ("Text ({} polygons, {} vertices, layers {}, datatypes {})").format(len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))


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
    """
    __slots__ = ('layers', 'datatypes', 'polygons', 'x', 'y', 'w', 'n', 'direction', 'distance', 'length')

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
        if isinstance(self.direction, basestring):
            self.direction = _directions_dict[self.direction] * numpy.pi
        self.direction += angle
        cur = numpy.array((self.x, self.y)) - c0
        self.x, self.y = cur * ca + cur[::-1] * sa + c0
        self.polygons = [(points - c0) * ca + (points - c0)[:, ::-1] * sa + c0 for points in self.polygons]
        return self

    def segment(self, length, direction=None, final_width=None, final_distance=None, axis_offset=0,
                layer=0, datatype=0):
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
        out : ``Path``
            This object.
        """
        if direction is None:
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
                    numpy.array([(old_x + (old_d0 - old_w) * sa, old_y - (old_d0 - old_w) * ca),
                                 (old_x + (old_d0 + old_w) * sa, old_y - (old_d0 + old_w) * ca),
                                 (self.x + (d0 + self.w) * sa, self.y - (d0 + self.w) * ca),
                                 (self.x + (d0 - self.w) * sa, self.y - (d0 - self.w) * ca)]))
                if self.w == 0:
                    self.polygons[-1] = self.polygons[-1][:-1, :]
                if old_w == 0:
                    self.polygons[-1] = self.polygons[-1][1:, :]
            self.length += numpy.sqrt(length**2 + axis_offset**2)
            if isinstance(layer, list):
                self.layers.extend((layer * (self.n // len(layer) + 1))[:self.n])
            else:
                self.layers.extend(layer for _ in range(self.n))
            if isinstance(datatype, list):
                self.datatypes.extend((datatype * (self.n // len(datatype) + 1))[:self.n])
            else:
                self.datatypes.extend(datatype for _ in range(self.n))
        return self

    def arc(self, radius, initial_angle, final_angle, tolerance=0.01, number_of_points=None,
            max_points=199, final_width=None, final_distance=None, layer=0, datatype=0):
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
        number_of_points : integer or ``None``
            Manually define the number of vertices that form the object
            (polygonal approximation).  Overrides ``tolerance``.
        max_points : integer
            If the number of points in the element is greater than
            ``max_points``, it will be fractured in smaller polygons
            with at most ``max_points`` each.  If ``max_points = 0`` no
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
        out : ``Path``
            This object.

        Notes
        -----
        The GDSII specification supports only a maximum of 199 vertices
        per polygon.
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
            warnings.warn("[GDSPY] Use of a floating number as number_of_points is deprecated in favor of tolerance.", category=DeprecationWarning, stacklevel=2)
            tolerance = number_of_points
            number_of_points = None
        if number_of_points is None:
            number_of_points = 2 * int(abs((final_angle - initial_angle) * (radius + max(old_distance, self.distance)
                                           * (self.n - 1) * 0.5 + max(old_w, self.w)) / tolerance) + 0.5) + 2
        number_of_points = max(number_of_points, 3)
        pieces = 1 if max_points == 0 else int(numpy.ceil(number_of_points / float(max_points)))
        number_of_points = number_of_points // pieces
        widths = numpy.linspace(old_w, self.w, pieces + 1)
        distances = numpy.linspace(old_distance, self.distance, pieces + 1)
        angles = numpy.linspace(initial_angle, final_angle, pieces + 1)
        if (self.w != 0) or (old_w != 0):
            for jj in range(pieces):
                for ii in range(self.n):
                    self.polygons.append(numpy.zeros((number_of_points, 2)))
                    r0 = radius + ii * distances[jj + 1] - (self.n - 1) * distances[jj + 1] * 0.5
                    old_r0 = radius + ii * distances[jj] - (self.n - 1) * distances[jj] * 0.5
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
                        self.polygons[-1][:pts1 - 1, :] = numpy.array(self.polygons[-1][1:pts1, :])
                        pts1 -= 1
                        pts2 += 1
                    ang = numpy.linspace(angles[jj + 1], angles[jj], pts2)
                    rad = numpy.linspace(r0 - widths[jj + 1], old_r0 - widths[jj], pts2)
                    if poly_warnings and (rad[0] <= 0 or rad[-1] <= 0):
                        warnings.warn("[GDSPY] Path arc with width larger than radius created: possible self-intersecting polygon.", stacklevel=2)
                    self.polygons[-1][pts1:, 0] = numpy.cos(ang) * rad + cx
                    self.polygons[-1][pts1:, 1] = numpy.sin(ang) * rad + cy
                self.length += abs((angles[jj + 1] - angles[jj]) * radius)
                if isinstance(layer, list):
                    self.layers.extend((layer * (self.n // len(layer) + 1))[:self.n])
                else:
                    self.layers.extend(layer for _ in range(self.n))
                if isinstance(datatype, list):
                    self.datatypes.extend((datatype * (self.n // len(datatype) + 1))[:self.n])
                else:
                    self.datatypes.extend(datatype for _ in range(self.n))
        return self

    def turn(self, radius, angle, tolerance=0.01, number_of_points=None, max_points=199,
             final_width=None, final_distance=None, layer=0, datatype=0):
        """
        Add a curved section to the path.

        Parameters
        ----------
        radius : number
            Central radius of the section.
        angle : {'r', 'l', 'rr', 'll'} or number
            Angle (in *radians*) of rotation of the path.  The values
            'r' and 'l' represent 90-degree turns cw and ccw,
            respectively; the values 'rr' and 'll' represent analogous
            180-degree turns.
        tolerance : float
            Approximate curvature resolution.  The number of points is
            automatically calculated.
        number_of_points : integer or ``None``
            Manually define the number of vertices that form the object
            (polygonal approximation).  Overrides ``tolerance``.
        max_points : integer
            If the number of points in the element is greater than
            ``max_points``, it will be fractured in smaller polygons
            with at most ``max_points`` each.  If ``max_points = 0`` no
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
        out : ``Path``
            This object.

        Notes
        -----
        The GDSII specification supports only a maximum of 199 vertices
        per polygon.
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
        self.arc(radius, self.direction + delta_i, self.direction + delta_f, tolerance,
                 number_of_points, max_points, final_width, final_distance, layer, datatype)
        if exact:
            self.direction = _directions_list[int(round(self.direction / _halfpi)) % 4]
        return self

    def parametric(self, curve_function, curve_derivative=None, tolerance=0.01, number_of_evaluations=5,
                   max_points=199, final_width=None, final_distance=None, layer=0, datatype=0):
        """
        Add a parametric curve to the path.

        The ``curve_function`` will be evaluated homogeneously in the
        interval [0, 1] at least ``number_of_points`` times.  More
        points will be added to the curve at the midpoint between
        evaluations if that points presents error larger than
        ``tolerance``.

        Parameters
        ----------
        curve_function : function
            Function that defines the curve.  Must be a function of one
            argument (that varies from 0 to 1) that returns a 2-element
            list, tuple or array (x, y).
        curve_derivative : function
            If set, it should be the derivative of the curve function.
            Must be a function of one argument (that varies from 0 to 1)
            that returns a 2-element list, tuple or array (x,y).  If
            ``None``, the derivative will be calculated numerically.
        tolerance : number
            Acceptable tolerance for the approximation of the curve
            function by a finite number of evaluations.
        number_of_evaluations : integer
            Initial number of points where the curve function will be
            evaluated.  According to ``tolerance``, more evaluations
            will be performed.
        max_points : integer
            Elements will be fractured until each polygon has at most
            ``max_points``.  If ``max_points = 0`` no fracture will
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
        out : ``Path``
            This object.

        Notes
        -----
        The norm of the vector returned by ``curve_derivative`` is not
        important.  Only the direction is used.

        The GDSII specification supports only a maximum of 199 vertices
        per polygon.

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
        values = list(numpy.array(curve_function(u)) for u in points)
        delta = points[1]
        i = 1
        while i < len(points):
            midpoint = 0.5 * (points[i] + points[i - 1])
            midvalue = numpy.array(curve_function(midpoint))
            if (((values[i] + values[i - 1]) / 2 - midvalue)**2).sum() > err:
                delta = min(delta, points[i] - midpoint)
                points.insert(i, midpoint)
                values.insert(i, midvalue)
            else:
                i += 1
        points = numpy.array(points)
        values = numpy.array(values)
        self.length += numpy.sqrt(((values[1:, :] - values[:-1, :])**2).sum(1)).sum()

        delta *= 0.5
        if curve_derivative is None:
            derivs = numpy.vstack((
                numpy.array(curve_function(delta)) - values[0],
                [numpy.array(curve_function(u + delta)) - numpy.array(curve_function(u - delta))
                    for u in points[1:-1]],
                values[-1] - numpy.array(curve_function(1 - delta))
            ))
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
        x0 = values + numpy.array((self.x, self.y))
        dx = derivs[:, ::-1] * numpy.array((-1, 1)) / numpy.sqrt((derivs**2).sum(1)).reshape(sh)
        width = width.reshape(sh)
        dist = dist.reshape(sh)

        self.x = x0[-1, 0]
        self.y = x0[-1, 1]
        self.direction = numpy.arctan2(-dx[-1, 0], dx[-1, 1])

        max_points = max(np, max_points // 2)
        i0 = 0
        while i0 < np - 1:
            i1 = min(i0 + max_points, np)
            for ii in range(self.n):
                p1 = x0[i0:i1] + dx[i0:i1] * (dist[i0:i1] * (ii - (self.n - 1) * 0.5) + width[i0:i1])
                p2 = (x0[i0:i1] + dx[i0:i1] * (dist[i0:i1] * (ii - (self.n - 1) * 0.5) - width[i0:i1]))[::-1, :]
                if width[i1 - 1, 0] == 0:
                    p2 = p2[1:]
                if width[i0, 0] == 0:
                    p1 = p1[1:]
                self.polygons.append(numpy.concatenate((p1, p2)))
            if isinstance(layer, list):
                self.layers.extend((layer * (self.n // len(layer) + 1))[:self.n])
            else:
                self.layers.extend(layer for _ in range(self.n))
            if isinstance(datatype, list):
                self.datatypes.extend((datatype * (self.n // len(datatype) + 1))[:self.n])
            else:
                self.datatypes.extend(datatype for _ in range(self.n))
            i0 = i1 - 1
        return self

    def bezier(self, points, tolerance=0.01, number_of_evaluations=5, max_points=199,
               final_width=None, final_distance=None, layer=0, datatype=0):
        """
        Add a Bezier curve to the path.

        Parameters
        ----------
        points : array-like[N][2]
            Control points defining the Bezier curve.
        tolerance : number
            Acceptable tolerance for the approximation of the curve
            function by a finite number of evaluations.
        number_of_evaluations : integer
            Initial number of points where the curve function will be
            evaluated.  According to ``tolerance``, more evaluations
            will be performed.
        max_points : integer
            Elements will be fractured until each polygon has at most
            ``max_points``.  If ``max_points = 0`` no fracture will
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
        out : ``Path``
            This object.

        Notes
        -----
        The GDSII specification supports only a maximum of 199 vertices
        per polygon.
        """

        def bez(ctrl):
            def _bez(u):
                p = ctrl
                for _ in range(ctrl.shape[0] - 1):
                    p = p[1:] * u + p[:-1] * (1 - u)
                return p[0]

            return _bez

        pts = numpy.array(points)
        dpts = (pts.shape[0] - 1) * (pts[1:] - pts[:-1])
        self.parametric(bez(pts), bez(dpts), tolerance, number_of_evaluations, max_points,
                        final_width, final_distance, layer, datatype)
        return self

    def smooth(self, points, angles=None, curl_start=1, curl_end=1, t_in=1, t_out=1, cycle=False,
               tolerance=0.01, number_of_evaluations=5, max_points=199, final_widths=None,
               final_distances=None, layer=0, datatype=0):
        """
        Add a smooth interpolating curve that passes through the given
        points.

        Uses the Hobby algorithm to calculate a smooth interpolating
        curve made of cubic Bezier segments between each pair of points.

        Reference: Hobby, J.D.  *Discrete Comput. Geom.* (1986) 1: 123.
        https://doi.org/10.1007/BF02187690

        Parameters
        ----------
        points : array-like[N][2]
            Vertices in the interpolating curve.
        angles : array-like[N] or ``None``
            Tangent angles at each point (in *radians*).  Any angles
            defined as ``None`` are automatically calculated.
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
        t_in : number or array-like[N]
            Tension parameter when arriving at each point.  One value
            per point or a single value used for all points.
        t_out : number or array-like[N]
            Tension parameter when leaving each point.  One value per
            point or a single value used for all points.
        cycle : bool
            If ``True``, calculates control points for a closed curve,
            with an additional segment connecting the first and last
            points.
        tolerance : number
            Acceptable tolerance for the approximation of the curve
            function by a finite number of evaluations.
        number_of_evaluations : integer
            Initial number of points where the curve function will be
            evaluated.  According to ``tolerance``, more evaluations
            will be performed.
        max_points : integer
            Elements will be fractured until each polygon has at most
            ``max_points``.  If ``max_points = 0`` no fracture will
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
        out : ``Path``
            This object.

        Notes
        -----
        The GDSII specification supports only a maximum of 199 vertices
        per polygon.
        """
        if not isinstance(points, numpy.ndarray):
            points = numpy.array(points)
        cta, ctb = _hobby(points, angles, curl_start, curl_end, t_in, t_out, cycle)
        if final_widths is None:
            final_widths = [None for _ in range(cta.shape[0])]
        if final_distances is None:
            final_distances = [None for _ in range(cta.shape[0])]
        x = self.x
        y = self.y
        for i in range(points.shape[0] - 1):
            self.x = x
            self.y = y
            self.bezier([points[i], cta[i], ctb[i], points[i + 1]], tolerance, number_of_evaluations,
                        max_points, final_widths[i], final_distances[i], layer, datatype)
        if cycle:
            self.x = x
            self.y = y
            self.bezier([points[-1], cta[-1], ctb[-1], points[0]], tolerance, number_of_evaluations,
                        max_points, final_widths[-1], final_distances[-1], layer, datatype)
        return self


class L1Path(PolygonSet):
    """
    Series of geometric objects that form a path or a collection of
    parallel paths with Manhattan geometry.

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
        have 1 element less then ``length``.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number
        Distance between the centers of adjacent paths.
    max_points : integer
        The paths will be fractured in polygons with at most
        ``max_points`` (must be at least 6).  If ``max_points = 0`` no
        fracture will occur.
    layer : integer, list
        The GDSII layer numbers for the elements of each path.  If the
        number of layers in the list is less than the number of paths,
        the list is repeated.
    datatype : integer, list
        The GDSII datatype for the elements of each path (between 0 and
        255).  If the number of datatypes in the list is less than the
        number of paths, the list is repeated.

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
    direction : '+x', '-x', '+y', '-y' or number
        Direction or angle (in *radians*) the path points to.  The
        numerical angle is returned only after a rotation of the object.

    Examples
    --------
    >>> length = [10, 30, 15, 15, 15, 15, 10]
    >>> turn = [1, -1, -1, 3, -1, 1]
    >>> l1path = gdspy.L1Path((0, 0), '+x', 2, length, turn)
    >>> myCell.add(l1path)
    """
    __slots__ = 'layers', 'datatypes', 'polygons', 'direction', 'x', 'y'

    def __init__(self, initial_point, direction, width, length, turn, number_of_paths=1, distance=0,
                 max_points=199, layer=0, datatype=0):
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
        self.direction = ['+x', '+y', '-x', '-y'][direction]
        self.polygons.extend(numpy.array(p[0] + p[1][::-1]) for p in paths)
        self.layers.extend(layer)
        self.datatypes.extend(datatype)

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
        if isinstance(self.direction, basestring):
            self.direction = _directions_dict[self.direction] * numpy.pi
        self.direction += angle
        cur = numpy.array((self.x, self.y)) - c0
        self.x, self.y = cur * ca + cur[::-1] * sa + c0
        self.polygons = [(points - c0) * ca + (points - c0)[:, ::-1] * sa + c0 for points in self.polygons]
        return self


class PolyPath(PolygonSet):
    """
    Series of geometric objects that form a polygonal path or a
    collection of parallel polygonal paths.

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
        ``max_points`` (must be at least 4).  If ``max_points = 0`` no
        fracture will occur.
    layer : integer, list
        The GDSII layer numbers for the elements of each path.  If the
        number of layers in the list is less than the number of paths,
        the list is repeated.
    datatype : integer, list
        The GDSII datatype for the elements of each path (between 0 and
        255).  If the number of datatypes in the list is less than the
        number of paths, the list is repeated.

    Returns
    -------
    out : ``PolyPath``
        This object.

    Notes
    -----
    The bevel join will give strange results if the number of paths is
    greater than 1.
    """
    __slots__ = 'layers', 'datatypes', 'polygons'

    def __init__(self, points, width, number_of_paths=1, distance=0, corners='miter', ends='flush',
                 max_points=199, layer=0, datatype=0):
        if not isinstance(layer, list):
            layer = [layer]
        if not isinstance(datatype, list):
            datatype = [datatype]
        if hasattr(width, '__iter__'):
            width = numpy.array(width) * 0.5
        else:
            width = numpy.array([width * 0.5])
        len_w = len(width)
        if not hasattr(distance, '__iter__'):
            distance = (distance, )
        len_d = len(distance)
        self.polygons = []
        self.layers = []
        self.datatypes = []
        points = numpy.array(points, dtype=float)
        if corners not in ['miter', 'bevel']:
            if corners in [0, 1]:
                corners = ['miter', 'bevel'][corners]
                warnings.warn("[GDSPY] Argument corners must be one of 'miter' or 'bevel'.", category=DeprecationWarning, stacklevel=2)
            else:
                raise ValueError("[GDSPY] Argument corners must be one of 'miter' or 'bevel'.")
        bevel = corners == 'bevel'
        if ends not in ['flush', 'round', 'extended']:
            if ends in [0, 1, 2]:
                ends = ['flush', 'round', 'extended'][ends]
                warnings.warn("[GDSPY] Argument ends must be one of 'flush', 'round', or 'extended'.", category=DeprecationWarning, stacklevel=2)
            else:
                raise ValueError("[GDSPY] Argument ends must be one of 'flush', 'round', or 'extended'.")
        if ends == 'extended':
            v = points[0, :] - points[1, :]
            v = v / numpy.sqrt(numpy.sum(v * v))
            points[0, :] = points[0, :] + v * width[0]
            v = points[-1, :] - points[-2, :]
            v = v / numpy.sqrt(numpy.sum(v * v))
            points[-1, :] = points[-1, :] + v * width[(points.shape[0] - 1) % len_w]
        elif ends == 'round':
            v0 = points[1, :] - points[0, :]
            angle0 = numpy.arctan2(v0[1], v0[0]) + _halfpi
            v0 = numpy.array((-v0[1], v0[0])) / numpy.sqrt(numpy.sum(v0 * v0))
            d0 = 0.5 * (number_of_paths - 1) * distance[0]
            v1 = points[-1, :] - points[-2, :]
            angle1 = numpy.arctan2(v1[1], v1[0]) - _halfpi
            v1 = numpy.array((-v1[1], v1[0])) / numpy.sqrt(numpy.sum(v1 * v1))
            j1w = (points.shape[0] - 1) % len_w
            j1d = (points.shape[0] - 1) % len_d
            d1 = 0.5 * (number_of_paths - 1) * distance[j1d]
            self.polygons.extend((Round(points[0, :] + v0 * (ii * distance[0] - d0), width[0],
                                        initial_angle=angle0, final_angle=angle0 + numpy.pi,
                                        number_of_points=33).polygons[0] for ii in range(number_of_paths)))
            self.polygons.extend((Round(points[-1, :] + v1 * (ii * distance[j1d] - d1), width[j1w],
                                        initial_angle=angle1, final_angle=angle1 + numpy.pi,
                                        number_of_points=33).polygons[0]) for ii in range(number_of_paths))
            self.layers.extend(((layer * (number_of_paths // len(layer) + 1))[:number_of_paths]) * 2)
            self.datatypes.extend(((datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths]) * 2)
        v = points[1, :] - points[0, :]
        v = numpy.array((-v[1], v[0])) / numpy.sqrt(numpy.sum(v * v))
        d0 = 0.5 * (number_of_paths - 1) * distance[0]
        d1 = 0.5 * (number_of_paths - 1) * distance[1 % len_d]
        paths = [[[points[0, :] + (ii * distance[0] - d0 - width[0]) * v],
                  [points[0, :] + (ii * distance[0] - d0 + width[0]) * v]] for ii in range(number_of_paths)]
        p1 = [(points[1, :] + (ii * distance[1 % len_d] - d1 - width[1 % len_w]) * v,
               points[1, :] + (ii * distance[1 % len_d] - d1 + width[1 % len_w]) * v) for ii in range(number_of_paths)]
        for jj in range(1, points.shape[0] - 1):
            j0d = jj % len_d
            j0w = jj % len_w
            j1d = (jj + 1) % len_d
            j1w = (jj + 1) % len_w
            v = points[jj + 1, :] - points[jj, :]
            v = numpy.array((-v[1], v[0])) / numpy.sqrt(numpy.sum(v * v))
            d0 = d1
            d1 = 0.5 * (number_of_paths - 1) * distance[j1d]
            p0 = p1
            p1 = []
            pp = []
            for ii in range(number_of_paths):
                pp.append((points[jj, :] + (ii * distance[j0d] - d0 - width[j0w]) * v,
                           points[jj, :] + (ii * distance[j0d] - d0 + width[j0w]) * v))
                p1.append((points[jj + 1, :] + (ii * distance[j1d] - d1 - width[j1w]) * v,
                           points[jj + 1, :] + (ii * distance[j1d] - d1 + width[j1w]) * v))
                for kk in (0, 1):
                    p0m = paths[ii][kk][-1] - p0[ii][kk]
                    p1p = pp[ii][kk] - p1[ii][kk]
                    vec = p0m[0] * p1p[1] - p1p[0] * p0m[1]
                    if abs(vec) > 1e-30:
                        p = numpy.array((1, -1)) * (p0m * p1p[::-1] * p1[ii][kk] - p1p * p0m[::-1] * p0[ii][kk] + p0m * p1p * (p0[ii][kk][::-1] - p1[ii][kk][::-1])) / vec
                        if bevel and numpy.sum((p - pp[ii][kk]) * p1p) > 0 and numpy.sum((p - p0[ii][kk]) * p0m) < 0:
                            paths[ii][kk].append(p0[ii][kk])
                            paths[ii][kk].append(pp[ii][kk])
                        else:
                            paths[ii][kk].append(p)
                if max_points > 0 and len(paths[ii][0]) + len(paths[ii][1]) + 3 > max_points:
                    if numpy.sum((paths[ii][0][0] - paths[ii][1][0])**2) == 0:
                        paths[ii][1] = paths[ii][1][1:]
                    if numpy.sum((paths[ii][0][-1] - paths[ii][1][-1])**2) == 0:
                        self.polygons.append(numpy.array(paths[ii][0] + paths[ii][1][-2::-1]))
                    else:
                        self.polygons.append(numpy.array(paths[ii][0] + paths[ii][1][::-1]))
                    paths[ii][0] = paths[ii][0][-1:]
                    paths[ii][1] = paths[ii][1][-1:]
                    self.layers.append(layer[ii % len(layer)])
                    self.datatypes.append(datatype[ii % len(datatype)])
        for ii in range(number_of_paths):
            if numpy.sum((paths[ii][0][0] - paths[ii][1][0])**2) == 0:
                paths[ii][1] = paths[ii][1][1:]
            if numpy.sum((p1[ii][0] - p1[ii][1])**2) != 0:
                paths[ii][0].append(p1[ii][0])
            paths[ii][1].append(p1[ii][1])
        self.polygons.extend(numpy.array(pol[0] + pol[1][::-1]) for pol in paths)
        self.layers.extend((layer * (number_of_paths // len(layer) + 1))[:number_of_paths])
        self.datatypes.extend((datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths])

    def __str__(self):
        return "PolyPath ({} polygons, {} vertices, layers {}, datatypes {})".format(len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))


class UPath:
    """
    Path object according to GDSII specification.

    This class is stored internaly as a sequence of points, instead of
    polygonal boundaries that compose a path.  It can be used when the
    width of the path is constant.


    offsets instead of distances, list of (number|callable|array), 1 per path
        number = constant
        array = 1 per vertex (extarpolate constant)
        callable = func(u = 0 thru #subpaths)

    path(u) = position
    grad(u) = velocity

    all methods from Path: subpaths
    PolyPath/L1Path = N subpaths

    get_polygons(corners, ends, tolerance, max_points)

    to_gds(as_path=True/False)

    from_svg()


    Parameters
    ----------
    points : array-like[N][2]
        Points along the center of the path.
    width : number, list
        Width of each parallel path being created.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number or array-like[N]
        Distance between the centers of adjacent paths.  If an array is
        given, distance at each endpoint.
    ends : 'flush', 'round', 'extended', or array-like[2]
        Type of end caps for the paths.  An array represents the start
        and end extensions to the paths.
    layer : integer, list
        The GDSII layer numbers for the elements of each path.  If the
        number of layers in the list is less than the number of paths,
        the list is repeated.
    datatype : integer, list
        The GDSII datatype for the elements of each path (between 0 and
        255).  If the number of datatypes in the list is less than the
        number of paths, the list is repeated.

    Returns
    -------
    out : ``SimplePath``
        This object.

    Notes
    -----
    Acute joins are automatically bevelled when drawing the path.  The
    amount of beveling is implementation-dependent (the GDSII file does
    not store this information).
    """
    __slots__ = 'layers', 'datatypes', 'paths', 'widths', 'ends'

    def __init__(self, points, width, number_of_paths=1, distance=0, ends='flush', layer=0, datatype=0):
        self.ends = ends

        if isinstance(width, list):
            self.widths = (width * (number_of_paths // len(width) + 1))[:number_of_paths]
        else:
            self.widths = [width for _ in range(number_of_paths)]

        if isinstance(layer, list):
            self.layers = (layer * (number_of_paths // len(layer) + 1))[:number_of_paths]
        else:
            self.layers = [layer for _ in range(number_of_paths)]

        if isinstance(datatype, list):
            self.datatypes = (datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths]
        else:
            self.datatypes = [datatype for _ in range(number_of_paths)]

        if not hasattr(distance, '__iter__'):
            distance = (distance, )
        len_d = len(distance)

        self.paths = []

        points = numpy.array(points, dtype=float)

    def __str__(self):
        pass
    def __repr__(self):
        pass
    def polygons():
        pass
    def to_gds():
        pass


class Label(object):
    """
    Text that can be used to label parts of the geometry or display
    messages.  The text does not create additional geometry, it's meant
    for display and labeling purposes only.

    Parameters
    ----------
    text : string
        The text of this label.
    position : array-like[2]
        Text anchor position.
    anchor : 'n', 's', 'e', 'w', 'o', 'ne', 'nw'...
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

    Attributes
    ----------
    text : string
        The text of this label.
    position : array-like[2]
        Text anchor position.
    anchor : int
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

    _anchor = {
        'nw':            0,
        'top left':      0,
        'upper left':    0,
        'n':             1,
        'top center':    1,
        'upper center':  1,
        'ne':            2,
        'top right':     2,
        'upper right':   2,
        'w':             4,
        'middle left':   4,
        'o':             5,
        'middle center': 5,
        'e':             6,
        'middle right':  6,
        'sw':            8,
        'bottom left':   8,
        'lower left':    8,
        's':             9,
        'bottom center': 9,
        'lower center':  9,
        'se':           10,
        'bottom right': 10,
        'lower right':  10,
    }

    __slots__ = ('layer', 'texttype', 'text', 'position', 'anchor', 'rotation', 'magnification', 'x_reflection')

    def __init__(self, text, position, anchor='o', rotation=None, magnification=None,
                 x_reflection=False, layer=0, texttype=0):
        self.layer = layer
        self.text = text
        self.position = numpy.array(position)
        self.anchor = Label._anchor.get(anchor.lower(), None)
        if self.anchor is None:
            raise ValueError("[GDSPY] Label anchors must be one of: '" + "', '".join(Label._anchor.keys()) + "'.")
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection
        self.texttype = texttype

    def __repr__(self):
        return "Label(\"{0}\", ({1[0]}, {1[1]}), {2}, {3}, {4}, {5}, {6})".format(self.text, self.position, self.rotation, self.magnification, self.x_reflection, self.layer, self.texttype)

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
        if len(text) % 2 != 0:
            text = text + '\0'
        data = struct.pack('>4Hh2Hh2Hh', 4, 0x0C00, 6, 0x0D02, self.layer, 6, 0x1602, self.texttype, 6, 0x1701,
                           self.anchor)
        if (self.rotation is not None) or (self.magnification is not None) or self.x_reflection:
            word = 0
            values = b''
            if self.x_reflection:
                word += 0x8000
            if not (self.magnification is None):
                # This flag indicates that the magnification is absolute, not
                # relative (not supported).
                # word += 0x0004
                values += struct.pack('>2H', 12, 0x1B05) + _eight_byte_real(self.magnification)
            if not (self.rotation is None):
                # This flag indicates that the rotation is absolute, not
                # relative (not supported).
                # word += 0x0002
                values += struct.pack('>2H', 12, 0x1C05) + _eight_byte_real(self.rotation)
            data += struct.pack('>3H', 6, 0x1A01, word) + values
        return data + struct.pack('>2H2l2H', 12, 0x1003, int(round(self.position[0] * multiplier)),
                                  int(round(self.position[1] * multiplier)), 4 + len(text),
                                  0x1906) + text.encode('ascii') + struct.pack('>2H', 4, 0x1100)

    def translate(self, dx, dy):
        """
        Move the text from one place to another

        Parameters
        ----------
        dx : float
            distance to move in the x-direction
        dy : float
            distance to move in the y-direction

        Returns
        -------
        out : ``Label``
            This object.

        Examples
        --------
        >>> text = gdspy.Label((0, 0), (10, 20))
        >>> text = text.translate(2, 0)
        >>> myCell.add(text)
        """
        self.position = numpy.array((dx + self.position[0], dy + self.position[1]))

        return self


class Cell(object):
    """
    Collection of elements, both geometric objects and references to
    other cells.

    Parameters
    ----------
    name : string
        The name of the cell.
    exclude_from_current : bool
        If ``True``, the cell will not be automatically included in the
        current library.

    Attributes
    ----------
    name : string
        The name of this cell.
    elements : list
        List of cell elements (``PolygonSet``, ``CellReference``,
        ``CellArray``).
    labels : list
        List of ``Label``.
    """
    __slots__ = 'name', 'elements', 'labels', '_bb_valid'

    def __init__(self, name, exclude_from_current=False):
        self.name = name
        self.elements = []
        self.labels = []
        self._bb_valid = False
        if not exclude_from_current:
            current_library.add(self)

    def __str__(self):
        return "Cell (\"{}\", {} elements, {} labels)".format(self.name, len(self.elements), len(self.labels))

    def to_gds(self, multiplier, timestamp=None):
        """
        Convert this cell to a GDSII structure.

        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            structure.
        timestamp : datetime object
            Sets the GDSII timestamp.  If ``None``, the current time is
            used.

        Returns
        -------
        out : string
            The GDSII binary string that represents this cell.
        """
        now = datetime.datetime.today() if timestamp is None else timestamp
        name = self.name
        if len(name) % 2 != 0:
            name = name + '\0'
        return (struct.pack('>2H12h2H', 28, 0x0502, now.year, now.month, now.day, now.hour, now.minute,
                            now.second, now.year, now.month, now.day, now.hour, now.minute, now.second,
                            4 + len(name), 0x0606)
                + name.encode('ascii')
                + b''.join(element.to_gds(multiplier) for element in self.elements)
                + b''.join(label.to_gds(multiplier) for label in self.labels)
                + struct.pack('>2H', 4, 0x0700))

    def copy(self, name, exclude_from_current=False, deep_copy=False):
        """
        Creates a copy of this cell.

        Parameters
        ----------
        name : string
            The name of the cell.
        exclude_from_current : bool
            If ``True``, the cell will not be included in the global
            list of cells maintained by ``gdspy``.
        deep_copy : bool
            If ``False``, the new cell will contain only references to
            the existing elements.  If ``True``, copies of all elements
            are also created.

        Returns
        -------
        out : ``Cell``
            The new copy of this cell.
        """
        new_cell = Cell(name, exclude_from_current)
        if deep_copy:
            new_cell.elements = libCopy.deepcopy(self.elements)
            new_cell.labels = libCopy.deepcopy(self.labels)
            for ref in new_cell.get_dependencies(True):
                if ref._bb_valid:
                    ref._bb_valid = False
        else:
            new_cell.elements = list(self.elements)
            new_cell.labels = list(self.labels)
        return new_cell

    def add(self, element):
        """
        Add a new element or list of elements to this cell.

        Parameters
        ----------
        element : ``PolygonSet``, ``CellReference``, ``CellArray`` or iterable
            The element or iterable of elements to be inserted in this
            cell.

        Returns
        -------
        out : ``Cell``
            This cell.
        """
        if isinstance(element, PolygonSet) or isinstance(element, CellReference) or isinstance(element, CellArray):
            self.elements.append(element)
        elif isinstance(element, Label):
            self.labels.append(element)
        else:
            for e in element:
                if isinstance(e, Label):
                    self.labels.append(e)
                else:
                    self.elements.append(e)
        self._bb_valid = False
        return self

    def remove_polygons(self, test):
        """
        Remove polygons from this cell.

        The function or callable ``test`` is called for each polygon in
        the cell.  If its return value evaluates to ``True``, the
        corresponding polygon is removed from the cell.

        Parameters
        ----------
        test : callable
            Test function to query whether a polygon should be removed.
            The function is called with arguments:
            ``(points, layer, datatype)``

        Returns
        -------
        out : ``Cell``
            This cell.

        Examples
        --------
        Remove polygons in layer 1:

        >>> cell.remove_polygons(lambda pts, layer, datatype:
        ...                      layer == 1)

        Remove polygons with negative x coordinates:

        >>> cell.remove_polygons(lambda pts, layer, datatype:
        ...                      any(pts[:, 0] < 0))
        """
        empty = []
        for element in self.elements:
            if isinstance(element, PolygonSet):
                ii = 0
                while ii < len(element.polygons):
                    if test(element.polygons[ii], element.layers[ii], element.datatypes[ii]):
                        element.polygons.pop(ii)
                        element.layers.pop(ii)
                        element.datatypes.pop(ii)
                    else:
                        ii += 1
                if len(element.polygons) == 0:
                    empty.append(element)
        for element in empty:
            self.elements.remove(element)
        return self

    def remove_labels(self, test):
        """
        Remove labels from this cell.

        The function or callable ``test`` is called for each label in
        the cell.  If its return value evaluates to ``True``, the
        corresponding label is removed from the cell.

        Parameters
        ----------
        test : callable
            Test function to query whether a label should be removed.
            The function is called with the label as the only argument.

        Returns
        -------
        out : ``Cell``
            This cell.

        Examples
        --------
        Remove labels in layer 1:

        >>> cell.remove_labels(lambda lbl: lbl.layer == 1)
        """
        ii = 0
        while ii < len(self.labels):
            if test(self.labels[ii]):
                self.labels.pop(ii)
            else:
                ii += 1
        return self

    def area(self, by_spec=False):
        """
        Calculate the total area of the elements on this cell, including
        cell references and arrays.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the areas
            of each individual pair (layer, datatype).

        Returns
        -------
        out : number, dictionary
            Area of this cell.
        """
        if by_spec:
            cell_area = {}
            for element in self.elements:
                element_area = element.area(True)
                for ll in element_area.keys():
                    if ll in cell_area:
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
        Returns a set of layers in this cell.

        Returns
        -------
        out : set
            Set of the layers used in this cell.
        """
        layers = set()
        for element in self.elements:
            if isinstance(element, PolygonSet):
                layers.update(element.layers)
            elif isinstance(element, CellReference) or isinstance(element, CellArray):
                layers.update(element.ref_cell.get_layers())
        for label in self.labels:
            layers.add(label.layer)
        return layers

    def get_datatypes(self):
        """
        Returns a set of datatypes in this cell.

        Returns
        -------
        out : set
            Set of the datatypes used in this cell.
        """
        datatypes = set()
        for element in self.elements:
            if isinstance(element, PolygonSet):
                datatypes.update(element.datatypes)
            elif isinstance(element, CellReference) or isinstance(element, CellArray):
                datatypes.update(element.ref_cell.get_datatypes())
        return datatypes

    def get_bounding_box(self):
        """
        Returns the bounding box for this cell.

        Returns
        -------
        out : Numpy array[2, 2] or ``None``
            Bounding box of this cell [[x_min, y_min], [x_max, y_max]],
            or ``None`` if the cell is empty.
        """
        if len(self.elements) == 0:
            return None
        if not (self._bb_valid and all(ref._bb_valid for ref in self.get_dependencies(True))):
            bb = numpy.array(((1e300, 1e300), (-1e300, -1e300)))
            all_polygons = []
            for element in self.elements:
                if isinstance(element, PolygonSet):
                    all_polygons.extend(element.polygons)
                elif isinstance(element, CellReference) or isinstance(element, CellArray):
                    element_bb = element.get_bounding_box()
                    if element_bb is not None:
                        bb[0, 0] = min(bb[0, 0], element_bb[0, 0])
                        bb[0, 1] = min(bb[0, 1], element_bb[0, 1])
                        bb[1, 0] = max(bb[1, 0], element_bb[1, 0])
                        bb[1, 1] = max(bb[1, 1], element_bb[1, 1])
            if len(all_polygons) > 0:
                all_points = numpy.concatenate(all_polygons).transpose()
                bb[0, 0] = min(bb[0, 0], all_points[0].min())
                bb[0, 1] = min(bb[0, 1], all_points[1].min())
                bb[1, 0] = max(bb[1, 0], all_points[0].max())
                bb[1, 1] = max(bb[1, 1], all_points[1].max())
            self._bb_valid = True
            _bounding_boxes[self] = bb
        return _bounding_boxes[self]

    def get_polygons(self, by_spec=False, depth=None):
        """
        Returns a list of polygons in this cell.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype).
        depth : integer or ``None``
            If not ``None``, defines from how many reference levels to
            retrieve polygons.  References below this level will result
            in a bounding box.  If ``by_spec`` is ``True`` the key will
            be the name of this cell.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons (if
            ``by_spec`` is ``True``).
        """
        if depth is not None and depth < 0:
            bb = self.get_bounding_box()
            if bb is None:
                return {} if by_spec else []
            pts = [numpy.array([(bb[0, 0], bb[0, 1]), (bb[0, 0], bb[1, 1]),
                                (bb[1, 0], bb[1, 1]), (bb[1, 0], bb[0, 1])])]
            polygons = {self.name: pts} if by_spec else pts
        else:
            if by_spec:
                polygons = {}
                for element in self.elements:
                    if isinstance(element, PolygonSet):
                        for ii in range(len(element.polygons)):
                            key = (element.layers[ii], element.datatypes[ii])
                            if key in polygons:
                                polygons[key].append(numpy.array(element.polygons[ii]))
                            else:
                                polygons[key] = [numpy.array(element.polygons[ii])]
                    else:
                        if depth is None:
                            next_depth = None
                        else:
                            next_depth = depth - 1
                        cell_polygons = element.get_polygons(True, next_depth)
                        for kk in cell_polygons.keys():
                            if kk in polygons:
                                polygons[kk].extend(cell_polygons[kk])
                            else:
                                polygons[kk] = cell_polygons[kk]
            else:
                polygons = []
                for element in self.elements:
                    if isinstance(element, PolygonSet):
                        for points in element.polygons:
                            polygons.append(numpy.array(points))
                    else:
                        if depth is None:
                            next_depth = None
                        else:
                            next_depth = depth - 1
                        polygons.extend(element.get_polygons(depth=next_depth))
        return polygons

    def get_labels(self, depth=None):
        """
        Returns a list with a copy of the labels in this cell.

        Parameters
        ----------
        depth : integer or ``None``
            If not ``None``, defines from how many reference levels to
            retrieve labels from.

        Returns
        -------
        out : list of ``Label``
            List containing the labels in this cell and its references.
        """
        labels = libCopy.deepcopy(self.labels)
        if depth is None or depth > 0:
            for element in self.elements:
                if depth is None:
                    next_depth = None
                else:
                    next_depth = depth - 1
                if isinstance(element, CellReference):
                    labels.extend(element.get_labels(next_depth))
                elif isinstance(element, CellArray):
                    labels.extend(element.get_labels(next_depth))
        return labels

    def get_dependencies(self, recursive=False):
        """
        Returns a list of the cells included in this cell as references.

        Parameters
        ----------
        recursive : bool
            If True returns cascading dependencies.

        Returns
        -------
        out : set of ``Cell``
            List of the cells referenced by this cell.
        """
        dependencies = set()
        for element in self.elements:
            if isinstance(element, CellReference) or isinstance(element, CellArray):
                if recursive:
                    dependencies.update(element.ref_cell.get_dependencies(True))
                dependencies.add(element.ref_cell)
        return dependencies

    def flatten(self, single_layer=None, single_datatype=None, single_texttype=None):
        """
        Flatten all ``CellReference`` and ``CellArray`` elements in this
        cell into real polygons and labels, instead of references.

        Parameters
        ----------
        single_layer : integer or None
            If not ``None``, all polygons will be transfered to the
            layer indicated by this number.
        single_datatype : integer or None
            If not ``None``, all polygons will be transfered to the
            datatype indicated by this number.
        single_datatype : integer or None
            If not ``None``, all labels will be transfered to the
            texttype indicated by this number.

        Returns
        -------
        out : ``Cell``
            This cell.
        """
        self.labels = self.get_labels()
        if single_layer is not None:
            for lbl in self.labels:
                lbl.layer = single_layer
        if single_texttype is not None:
            for lbl in self.labels:
                lbl.texttype = single_texttype
        if single_layer is None or single_datatype is None:
            poly_dic = self.get_polygons(True)
            self.elements = []
            if single_layer is None and single_datatype is None:
                for ld in poly_dic.keys():
                    self.add(PolygonSet(poly_dic[ld], *ld))
            elif single_layer is None:
                for ld in poly_dic.keys():
                    self.add(PolygonSet(poly_dic[ld], ld[0], single_datatype))
            else:
                for ld in poly_dic.keys():
                    self.add(PolygonSet(poly_dic[ld], single_layer, ld[1]))
        else:
            polygons = self.get_polygons()
            self.elements = []
            self.add(PolygonSet(polygons, single_layer, single_datatype))
        return self


class CellReference(object):
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
        If ``True`` the reference is reflected parallel to the x
        direction before being rotated.
    ignore_missing : bool
        If ``False`` a warning is issued when the referenced cell is not
        found.
    """
    __slots__ = ('ref_cell', 'origin', 'rotation', 'magnification', 'x_reflection')

    def __init__(self, ref_cell, origin=(0, 0), rotation=None, magnification=None,
                 x_reflection=False, ignore_missing=False):
        self.origin = origin
        self.ref_cell = current_library.cell_dict.get(ref_cell, ref_cell)
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection
        if not isinstance(self.ref_cell, Cell) and not ignore_missing:
            warnings.warn("[GDSPY] Cell {0} not found; operations on this CellReference may not work.".format(self.ref_cell), stacklevel=2)

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
        if len(name) % 2 != 0:
            name = name + '\0'
        data = struct.pack('>4H', 4, 0x0A00, 4 + len(name), 0x1206) + name.encode('ascii')
        if (self.rotation is not None) or (self.magnification is not None) or self.x_reflection:
            word = 0
            values = b''
            if self.x_reflection:
                word += 0x8000
            if not (self.magnification is None):
                # This flag indicates that the magnification is absolute, not
                # relative (not supported).
                # word += 0x0004
                values += struct.pack('>2H', 12, 0x1B05) + _eight_byte_real(self.magnification)
            if not (self.rotation is None):
                # This flag indicates that the rotation is absolute, not
                # relative (not supported).
                # word += 0x0002
                values += struct.pack('>2H', 12, 0x1C05) + _eight_byte_real(self.rotation)
            data += struct.pack('>3H', 6, 0x1A01, word) + values
        return data + struct.pack('>2H2l2H', 12, 0x1003, int(round(self.origin[0] * multiplier)),
                                  int(round(self.origin[1] * multiplier)), 4, 0x1100)

    def area(self, by_spec=False):
        """
        Calculate the total area of the referenced cell with the
        magnification factor included.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the areas
            of each individual pair (layer, datatype).

        Returns
        -------
        out : number, dictionary
            Area of this cell.
        """
        if not isinstance(self.ref_cell, Cell):
            return dict() if by_spec else 0
        if self.magnification is None:
            return self.ref_cell.area(by_spec)
        else:
            if by_spec:
                factor = self.magnification**2
                cell_area = self.ref_cell.area(True)
                for kk in cell_area.keys():
                    cell_area[kk] *= factor
                return cell_area
            else:
                return self.ref_cell.area() * self.magnification**2

    def get_polygons(self, by_spec=False, depth=None):
        """
        Returns a list of polygons created by this reference.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype).
        depth : integer or ``None``
            If not ``None``, defines from how many reference levels to
            retrieve polygons.  References below this level will result
            in a bounding box.  If ``by_spec`` is ``True`` the key will
            be the name of the referenced cell.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons (if
            ``by_spec`` is ``True``).
        """
        if not isinstance(self.ref_cell, Cell):
            return dict() if by_spec else []
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0)
            st = numpy.array([-st, st])
        if self.x_reflection:
            xrefl = numpy.array([1, -1], dtype='int')
        if self.magnification is not None:
            mag = numpy.array([self.magnification, self.magnification])
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        if by_spec:
            polygons = self.ref_cell.get_polygons(True, depth)
            for kk in polygons.keys():
                for ii in range(len(polygons[kk])):
                    if self.x_reflection:
                        polygons[kk][ii] = polygons[kk][ii] * xrefl
                    if self.magnification is not None:
                        polygons[kk][ii] = polygons[kk][ii] * mag
                    if self.rotation is not None:
                        polygons[kk][ii] = (polygons[kk][ii] * ct + polygons[kk][ii][:, ::-1] * st)
                    if self.origin is not None:
                        polygons[kk][ii] = polygons[kk][ii] + orgn
        else:
            polygons = self.ref_cell.get_polygons(depth=depth)
            for ii in range(len(polygons)):
                if self.x_reflection:
                    polygons[ii] = polygons[ii] * xrefl
                if self.magnification is not None:
                    polygons[ii] = polygons[ii] * mag
                if self.rotation is not None:
                    polygons[ii] = (polygons[ii] * ct + polygons[ii][:, ::-1] * st)
                if self.origin is not None:
                    polygons[ii] = polygons[ii] + orgn
        return polygons

    def get_labels(self, depth=None):
        """
        Returns a list of labels created by this reference.

        Parameters
        ----------
        depth : integer or ``None``
            If not ``None``, defines from how many reference levels to
            retrieve labels from.

        Returns
        -------
        out : list of ``Label``
            List containing the labels in this cell and its references.
        """
        if not isinstance(self.ref_cell, Cell):
            return []
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0)
            st = numpy.array([-st, st])
        if self.x_reflection:
            xrefl = numpy.array([1, -1], dtype='int')
        if self.magnification is not None:
            mag = numpy.array([self.magnification, self.magnification])
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        labels = self.ref_cell.get_labels(depth=depth)
        for lbl in labels:
            if self.x_reflection:
                lbl.position = lbl.position * xrefl
            if self.magnification is not None:
                lbl.position = lbl.position * mag
            if self.rotation is not None:
                lbl.position = lbl.position * ct + lbl.position[::-1] * st
            if self.origin is not None:
                lbl.position = lbl.position + orgn
        return labels

    def get_bounding_box(self):
        """
        Returns the bounding box for this reference.

        Returns
        -------
        out : Numpy array[2, 2] or ``None``
            Bounding box of this cell [[x_min, y_min], [x_max, y_max]],
            or ``None`` if the cell is empty.
        """
        if not isinstance(self.ref_cell, Cell):
            return None
        if self.rotation is None and self.magnification is None and self.x_reflection is None:
            key = self
        else:
            key = (self.ref_cell, self.rotation, self.magnification, self.x_reflection)
        deps = self.ref_cell.get_dependencies(True)
        if not (self.ref_cell._bb_valid and all(ref._bb_valid for ref in deps) and key in _bounding_boxes):
            for ref in deps:
                ref.get_bounding_box()
            self.ref_cell.get_bounding_box()
            tmp = self.origin
            self.origin = None
            polygons = self.get_polygons()
            self.origin = tmp
            if len(polygons) == 0:
                bb = None
            else:
                all_points = numpy.concatenate(polygons).transpose()
                bb = numpy.array(((all_points[0].min(), all_points[1].min()),
                                  (all_points[0].max(), all_points[1].max())))
            _bounding_boxes[key] = bb
        else:
            bb = _bounding_boxes[key]
        if self.origin is None or bb is None:
            return bb
        else:
            return bb + numpy.array(((self.origin[0], self.origin[1]), (self.origin[0], self.origin[1])))

    def translate(self, dx, dy):
        """
        Move the reference from one place to another

        Parameters
        ----------
        dx : float
            distance to move in the x-direction
        dy : float
            distance to move in the y-direction

        Returns
        -------
        out : ``CellReference``
            This object.
        """
        self.origin = (self.origin[0] + dx, self.origin[1] + dy)
        return self


class CellArray(object):
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
        If ``True``, the reference is reflected parallel to the x
        direction before being rotated.
    ignore_missing : bool
        If ``False`` a warning is issued when the referenced cell is not
        found.
    """
    __slots__ = ('ref_cell', 'origin', 'rotation', 'magnification', 'x_reflection', 'columns', 'rows', 'spacing')

    def __init__(self, ref_cell, columns, rows, spacing, origin=(0, 0), rotation=None,
                 magnification=None, x_reflection=False, ignore_missing=False):
        self.columns = columns
        self.rows = rows
        self.spacing = spacing
        self.origin = origin
        self.ref_cell = current_library.cell_dict.get(ref_cell, ref_cell)
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection
        if not isinstance(self.ref_cell, Cell) and not ignore_missing:
            warnings.warn("[GDSPY] Cell {0} not found; operations on this CellArray may not work.".format(self.ref_cell), stacklevel=2)

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
        if len(name) % 2 != 0:
            name = name + '\0'
        data = struct.pack('>4H', 4, 0x0B00, 4 + len(name), 0x1206) + name.encode('ascii')
        x2 = self.origin[0] + self.columns * self.spacing[0]
        y2 = self.origin[1]
        x3 = self.origin[0]
        y3 = self.origin[1] + self.rows * self.spacing[1]
        if (self.rotation is not None) or (self.magnification is not None) or self.x_reflection:
            word = 0
            values = b''
            if self.x_reflection:
                word += 0x8000
                y3 = 2 * self.origin[1] - y3
            if not (self.magnification is None):
                # This flag indicates that the magnification is absolute, not
                # relative (not supported).
                # word += 0x0004
                values += struct.pack('>2H', 12, 0x1B05) + _eight_byte_real(self.magnification)
            if not (self.rotation is None):
                # This flag indicates that the rotation is absolute, not
                # relative (not supported).
                # word += 0x0002
                sa = numpy.sin(self.rotation * numpy.pi / 180.0)
                ca = numpy.cos(self.rotation * numpy.pi / 180.0)
                tmp = (x2 - self.origin[0]) * ca - (y2 - self.origin[1]) * sa + self.origin[0]
                y2 = (x2 - self.origin[0]) * sa + (y2 - self.origin[1]) * ca + self.origin[1]
                x2 = tmp
                tmp = (x3 - self.origin[0]) * ca - (y3 - self.origin[1]) * sa + self.origin[0]
                y3 = (x3 - self.origin[0]) * sa + (y3 - self.origin[1]) * ca + self.origin[1]
                x3 = tmp
                values += struct.pack('>2H', 12, 0x1C05) + _eight_byte_real(self.rotation)
            data += struct.pack('>3H', 6, 0x1A01, word) + values
        return data + struct.pack('>2H2h2H6l2H', 8, 0x1302, self.columns, self.rows, 28, 0x1003,
                                  int(round(self.origin[0] * multiplier)), int(round(self.origin[1] * multiplier)),
                                  int(round(x2 * multiplier)), int(round(y2 * multiplier)), int(round(x3 * multiplier)),
                                  int(round(y3 * multiplier)), 4, 0x1100)

    def area(self, by_spec=False):
        """
        Calculate the total area of the cell array with the
        magnification factor included.

        Parameters
        ----------
        by_spec : bool
            If ``True``, the return value is a dictionary with the areas
            of each individual pair (layer, datatype).

        Returns
        -------
        out : number, dictionary
            Area of this cell.
        """
        if not isinstance(self.ref_cell, Cell):
            return dict() if by_spec else 0
        if self.magnification is None:
            factor = self.columns * self.rows
        else:
            factor = self.columns * self.rows * self.magnification**2
        if by_spec:
            cell_area = self.ref_cell.area(True)
            for kk in cell_area.keys():
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
            If ``True``, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype).
        depth : integer or ``None``
            If not ``None``, defines from how many reference levels to
            retrieve polygons.  References below this level will result
            in a bounding box.  If ``by_spec`` is ``True`` the key will
            be name of the referenced cell.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons (if
            ``by_spec`` is ``True``).
        """
        if not isinstance(self.ref_cell, Cell):
            return dict() if by_spec else []
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0)
            st = numpy.array([-st, st])
        if self.magnification is not None:
            mag = numpy.array([self.magnification, self.magnification])
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        if self.x_reflection:
            xrefl = numpy.array([1, -1], dtype='int')
        if by_spec:
            cell_polygons = self.ref_cell.get_polygons(True, depth)
            polygons = {}
            for kk in cell_polygons.keys():
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
                                polygons[kk][-1] = polygons[kk][-1] * xrefl
                            if self.rotation is not None:
                                polygons[kk][-1] = (polygons[kk][-1] * ct + polygons[kk][-1][:, ::-1] * st)
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
                            polygons[-1] = polygons[-1] * xrefl
                        if self.rotation is not None:
                            polygons[-1] = (polygons[-1] * ct + polygons[-1][:, ::-1] * st)
                        if self.origin is not None:
                            polygons[-1] = polygons[-1] + orgn
        return polygons

    def get_labels(self, depth=None):
        """
        Returns a list of labels created by this reference.

        Parameters
        ----------
        depth : integer or ``None``
            If not ``None``, defines from how many reference levels to
            retrieve labels from.

        Returns
        -------
        out : list of ``Label``
            List containing the labels in this cell and its references.
        """
        if not isinstance(self.ref_cell, Cell):
            return []
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0)
            st = numpy.array([-st, st])
        if self.magnification is not None:
            mag = numpy.array([self.magnification, self.magnification])
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        if self.x_reflection:
            xrefl = numpy.array([1, -1], dtype='int')
        cell_labels = self.ref_cell.get_labels(depth=depth)
        labels = []
        for ii in range(self.columns):
            for jj in range(self.rows):
                spc = numpy.array([self.spacing[0] * ii, self.spacing[1] * jj])
                for clbl in cell_labels:
                    lbl = libCopy.deepcopy(clbl)
                    if self.magnification:
                        lbl.position = lbl.position * mag + spc
                    else:
                        lbl.position = lbl.position + spc
                    if self.x_reflection:
                        lbl.position = lbl.position * xrefl
                    if self.rotation is not None:
                        lbl.position = (lbl.position * ct + lbl.position[::-1] * st)
                    if self.origin is not None:
                        lbl.position = lbl.position + orgn
                    labels.append(lbl)
        return labels

    def get_bounding_box(self):
        """
        Returns the bounding box for this reference.

        Returns
        -------
        out : Numpy array[2, 2] or ``None``
            Bounding box of this cell [[x_min, y_min], [x_max, y_max]],
            or ``None`` if the cell is empty.
        """
        if not isinstance(self.ref_cell, Cell):
            return None
        key = (self.ref_cell, self.rotation, self.magnification, self.x_reflection, self.columns, self.rows,
               self.spacing[0], self.spacing[1])
        deps = self.ref_cell.get_dependencies(True)
        if not (self.ref_cell._bb_valid and all(ref._bb_valid for ref in deps) and key in _bounding_boxes):
            for ref in deps:
                ref.get_bounding_box()
            self.ref_cell.get_bounding_box()
            tmp = self.origin
            self.origin = None
            polygons = self.get_polygons()
            self.origin = tmp
            if len(polygons) == 0:
                bb = None
            else:
                all_points = numpy.concatenate(polygons).transpose()
                bb = numpy.array(((all_points[0].min(), all_points[1].min()),
                                  (all_points[0].max(), all_points[1].max())))
            _bounding_boxes[key] = bb
        else:
            bb = _bounding_boxes[key]
        if self.origin is None or bb is None:
            return bb
        else:
            return bb + numpy.array(((self.origin[0], self.origin[1]), (self.origin[0], self.origin[1])))

    def translate(self, dx, dy):
        """
        Move the reference from one place to another

        Parameters
        ----------
        dx : float
            distance to move in the x-direction
        dy : float
            distance to move in the y-direction

        Returns
        -------
        out : ``CellArray``
            This object.
        """
        self.origin = (self.origin[0] + dx, self.origin[1] + dy)
        return self


def _record_reader(stream):
    """
    Iterator over complete records from a GDSII stream file.

    Parameters
    ----------
    stream : file
        GDSII stream file to be read.

    Returns
    -------
    out : list
        Record type and data (as a numpy.array, string or None)
    """
    while True:
        header = stream.read(4)
        if len(header) < 4:
            return
        size, rec_type = struct.unpack('>HH', header)
        data_type = (rec_type & 0x00ff)
        rec_type = rec_type // 256
        data = None
        if size > 4:
            if data_type == 0x01:
                data = numpy.array(struct.unpack('>{0}H'.format((size - 4) // 2), stream.read(size - 4)), dtype='uint')
            elif data_type == 0x02:
                data = numpy.array(struct.unpack('>{0}h'.format((size - 4) // 2), stream.read(size - 4)), dtype='int')
            elif data_type == 0x03:
                data = numpy.array(struct.unpack('>{0}l'.format((size - 4) // 4), stream.read(size - 4)), dtype='int')
            elif data_type == 0x05:
                data = numpy.array([_eight_byte_real_to_float(stream.read(8)) for _ in range((size - 4) // 8)])
            else:
                data = stream.read(size - 4)
                if str is not bytes:
                    if data[-1] == 0:
                        data = data[:-1].decode('ascii')
                    else:
                        data = data.decode('ascii')
                elif data[-1] == '\0':
                    data = data[:-1]
        yield [rec_type, data]


class GdsLibrary(object):
    """
    GDSII library (file).

    Represent a GDSII library containing a dictionary of cells.

    Parameters
    ----------
    name : string
        Name of the GDSII library.  Ignored if a name is defined in
        `infile`.
    infile : file or string
        GDSII stream file (or path) to be imported.  It must be opened
        for reading in binary format.
    kwargs : keyword arguments
        Arguments passed to `read_gds`.

    Attributes
    ----------
    name : string
        Name of the GDSII library.
    cell_dict : dictionary
        Dictionary of cells in this library, indexed by name.
    unit : number
        Unit size for the objects in the library (in *meters*).
    precision : number
        Precision for the dimensions of the objects in the library (in
        *meters*).
    """

    _record_name = ('HEADER', 'BGNLIB', 'LIBNAME', 'UNITS', 'ENDLIB', 'BGNSTR', 'STRNAME', 'ENDSTR',
                    'BOUNDARY', 'PATH', 'SREF', 'AREF', 'TEXT', 'LAYER', 'DATATYPE', 'WIDTH', 'XY',
                    'ENDEL', 'SNAME', 'COLROW', 'TEXTNODE', 'NODE', 'TEXTTYPE', 'PRESENTATION',
                    'SPACING', 'STRING', 'STRANS', 'MAG', 'ANGLE', 'UINTEGER', 'USTRING', 'REFLIBS',
                    'FONTS', 'PATHTYPE', 'GENERATIONS', 'ATTRTABLE', 'STYPTABLE', 'STRTYPE', 'ELFLAGS',
                    'ELKEY', 'LINKTYPE', 'LINKKEYS', 'NODETYPE', 'PROPATTR', 'PROPVALUE', 'BOX',
                    'BOXTYPE', 'PLEX', 'BGNEXTN', 'ENDTEXTN', 'TAPENUM', 'TAPECODE', 'STRCLASS',
                    'RESERVED', 'FORMAT', 'MASK', 'ENDMASKS', 'LIBDIRSIZE', 'SRFNAME', 'LIBSECUR')
    _unused_records = (0x05, 0x00, 0x01, 0x02, 0x034, 0x38)
    _import_anchors = ['nw', 'n', 'ne', None, 'w', 'o', 'e', None, 'sw', 's', 'se']

    __slots__ = 'name', 'cell_dict', 'unit', 'precision', '_references'

    def __init__(self, name='library', infile=None, unit=1e-6, precision=1e-9, **kwargs):
        self.name = name
        self.cell_dict = {}
        self.unit = unit
        self.precision = precision
        if infile is not None:
            self.read_gds(infile, **kwargs)

    def __str__(self):
        return "GdsLibrary (" + ", ".join([c for c in self.cell_dict]) + ")"

    def add(self, cell, overwrite_duplicate=False):
        """
        Add one or more cells to the library.

        Parameters
        ----------
        cell : ``Cell`` or iterable
            Cells to be included in the library.
        overwrite_duplicate : bool
            If True an existing cell with the same name in the library
            will be overwritten.

        Returns
        -------
        out : ``GdsLibrary``
            This object.

        Notes
        -----
        ``CellReference`` or ``CellArray`` instances that referred to an
        overwritten cell are not automatically updated.
        """
        if isinstance(cell, Cell):
            if not overwrite_duplicate and cell.name in self.cell_dict and self.cell_dict[cell.name] is not cell:
                raise ValueError("[GDSPY] cell named {0} already present in library.".format(cell.name))
            self.cell_dict[cell.name] = cell
        else:
            for c in cell:
                if not overwrite_duplicate and c.name in self.cell_dict and self.cell_dict[c.name] is not c:
                    raise ValueError("[GDSPY] cell named {0} already present in library.".format(c.name))
                self.cell_dict[c.name] = c
        return self

    def write_gds(self, outfile, cells=None, timestamp=None, binary_cells=None):
        """
        Write the GDSII library to a file.

        The dimensions actually written on the GDSII file will be the
        dimensions of the objects created times the ratio
        ``unit/precision``.  For example, if a circle with radius 1.5 is
        created and we set ``unit=1.0e-6`` (1 um) and
        ``precision=1.0e-9`` (1 nm), the radius of the circle will be
        1.5 um and the GDSII file will contain the dimension 1500 nm.

        Parameters
        ----------
        outfile : file or string
            The file (or path) where the GDSII stream will be written.
            It must be opened for writing operations in binary format.
        cells : iterable
            The cells or cell names to be included in the library.  If
            ``None``, all cells are used.
        timestamp : datetime object
            Sets the GDSII timestamp.  If ``None``, the current time is
            used.
        binary_cells : iterable of bytes
            Iterable with binary data for GDSII cells (from
            ``get_binary_cells``, for example).

        Notes
        -----
        Only the specified cells are written.  The user is responsible
        for ensuring all cell dependencies are satisfied.
        """
        if isinstance(outfile, basestring):
            outfile = open(outfile, 'wb')
            close = True
        else:
            close = False
        now = datetime.datetime.today() if timestamp is None else timestamp
        name = self.name if len(self.name) % 2 == 0 else (self.name + '\0')
        outfile.write(struct.pack('>5H12h2H', 6, 0x0002, 0x0258, 28, 0x0102, now.year, now.month,
                                  now.day, now.hour, now.minute, now.second, now.year, now.month,
                                  now.day, now.hour, now.minute, now.second, 4 + len(name), 0x0206)
                      + name.encode('ascii') + struct.pack('>2H', 20, 0x0305)
                      + _eight_byte_real(self.precision / self.unit) + _eight_byte_real(self.precision))
        if cells is None:
            cells = self.cell_dict.values()
        else:
            cells = [self.cell_dict.get(c, c) for c in cells]
        for cell in cells:
            outfile.write(cell.to_gds(self.unit / self.precision))
        if binary_cells is not None:
            for bc in binary_cells:
                outfile.write(bc)
        outfile.write(struct.pack('>2H', 4, 0x0400))
        if close:
            outfile.close()

    def read_gds(self, infile, units='skip', rename={}, rename_template='{name}', layers={},
                 datatypes={}, texttypes={}):
        """
        Read a GDSII file into this library.

        Parameters
        ----------
        infile : file or string
            GDSII stream file (or path) to be imported.  It must be
            opened for reading in binary format.
        units : {'convert', 'import', 'skip'}
            Controls how to scale and use the units in the imported
            file.  ``'convert'``: the imported geometry is scaled to
            this library units. ``'import'``: the unit and precision in
            this library are replaced by those from the imported file.
            ``'skip'``: the imported geometry is not scaled and units
            are not replaced; the geometry is imported in the *user
            units* of the file.
        rename : dictionary
            Dictionary used to rename the imported cells.  Keys and
            values must be strings.
        rename_template : string
            Template string used to rename the imported cells. Appiled
            only if the cell name is not in the ``rename`` dictionary.
            Examples: ``'prefix-{name}'``, ``'{name}-suffix'``
        layers : dictionary
            Dictionary used to convert the layers in the imported cells.
            Keys and values must be integers.
        datatypes : dictionary
            Dictionary used to convert the datatypes in the imported
            cells.  Keys and values must be integers.
        texttypes : dictionary
            Dictionary used to convert the text types in the imported
            cells.  Keys and values must be integers.

        Returns
        -------
        out : ``GdsLibrary``
            This object.

        Notes
        -----
        Not all features from the GDSII specification are currently
        supported.  A warning will be produced if any unsuported
        features are found in the imported file.
        """
        self._references = []
        if isinstance(infile, basestring):
            infile = open(infile, 'rb')
            close = True
        else:
            close = False
        emitted_warnings = []
        kwargs = {}
        create_element = None
        factor = 1
        cell = None
        for record in _record_reader(infile):
            # LAYER
            if record[0] == 0x0d:
                kwargs['layer'] = layers.get(record[1][0], record[1][0])
            # DATATYPE
            elif record[0] == 0x0e:
                kwargs['datatype'] = datatypes.get(record[1][0], record[1][0])
            # TEXTTYPE
            elif record[0] == 0x16:
                kwargs['texttype'] = texttypes.get(record[1][0], record[1][0])
            # XY
            elif record[0] == 0x10:
                if 'xy' in kwargs:
                    kwargs['xy'] = numpy.concatenate((kwargs['xy'], factor * record[1]))
                else:
                    kwargs['xy'] = factor * record[1]
            # WIDTH
            elif record[0] == 0x0f:
                kwargs['width'] = factor * abs(record[1][0])
                if record[1][0] < 0 and record[0] not in emitted_warnings:
                    warnings.warn("[GDSPY] Paths with absolute width value are not supported.  Scaling these paths will also scale their width.", stacklevel=2)
                    emitted_warnings.append(record[0])
            # ENDEL
            elif record[0] == 0x11:
                if create_element is not None:
                    cell.add(create_element(**kwargs))
                    create_element = None
                kwargs = {}
            # BOUNDARY
            elif record[0] == 0x08:
                create_element = self._create_polygon
            # PATH
            elif record[0] == 0x09:
                create_element = self._create_path
            # TEXT
            elif record[0] == 0x0c:
                create_element = self._create_label
            # SNAME
            elif record[0] == 0x12:
                if record[1] in rename:
                    name = rename[record[1]]
                else:
                    name = rename_template.format(name=record[1])
                kwargs['ref_cell'] = name
            # COLROW
            elif record[0] == 0x13:
                kwargs['columns'] = record[1][0]
                kwargs['rows'] = record[1][1]
            # STRANS
            elif record[0] == 0x1a:
                kwargs['x_reflection'] = ((int(record[1][0]) & 0x8000) > 0)
                if (int(record[1][0]) & 0x0006) and record[0] not in emitted_warnings:
                    warnings.warn("[GDSPY] Absolute magnification or rotation of references is not supported.  Transformations will be interpreted as relative.", stacklevel=2)
                    emitted_warnings.append(record[0])
            # MAG
            elif record[0] == 0x1b:
                kwargs['magnification'] = record[1][0]
            # ANGLE
            elif record[0] == 0x1c:
                kwargs['rotation'] = record[1][0]
            # SREF
            elif record[0] == 0x0a:
                create_element = self._create_reference
            # AREF
            elif record[0] == 0x0b:
                create_element = self._create_array
            # STRNAME
            elif record[0] == 0x06:
                if record[1] in rename:
                    name = rename[record[1]]
                else:
                    name = rename_template.format(name=record[1])
                cell = Cell(name, exclude_from_current=True)
                self.cell_dict[name] = cell
            # STRING
            elif record[0] == 0x19:
                kwargs['text'] = record[1]
            # ENDSTR
            elif record[0] == 0x07:
                cell = None
            # UNITS
            elif record[0] == 0x03:
                if units == 'skip':
                    factor = record[1][0]
                elif units == 'import':
                    self.unit = record[1][1] / record[1][0]
                    self.precision = record[1][1]
                    factor = record[1][0]
                elif units == 'convert':
                    factor = record[1][1] / self.unit
                else:
                    raise ValueError("[GDSPY] units must be one of 'convert', 'import' or 'skip'.")
            # LIBNAME
            elif record[0] == 0x02:
                self.name = record[1]
            # PRESENTATION
            elif record[0] == 0x17:
                kwargs['anchor'] = GdsLibrary._import_anchors[int(record[1][0]) & 0x000f]
            # PATHTYPE
            elif record[0] == 0x21:
                if record[1][0] > 2:
                    # Custom path extension is manually added to the
                    # first and last path segments.
                    kwargs['ends'] = 0
                else:
                    kwargs['ends'] = record[1][0]
            # BGNEXTN
            elif record[0] == 0x30:
                kwargs['bgnextn'] = factor * record[1][0]
            # ENDEXTN
            elif record[0] == 0x31:
                kwargs['endextn'] = factor * record[1][0]
            # ENDLIB
            elif record[0] == 0x04:
                for ref in self._references:
                    if ref.ref_cell in self.cell_dict:
                        ref.ref_cell = self.cell_dict[ref.ref_cell]
                    elif ref.ref_cell in current_library.cell_dict:
                        ref.ref_cell = current_library.cell_dict[ref.ref_cell]
            # Not supported
            elif record[0] not in emitted_warnings and record[0] not in GdsLibrary._unused_records:
                warnings.warn("[GDSPY] Record type {0} ({1:02X}) is not supported.".format(GdsLibrary._record_name[record[0]], record[0]), stacklevel=2)
                emitted_warnings.append(record[0])
        if close:
            infile.close()
        return self

    def _create_polygon(self, layer, datatype, xy):
        return Polygon(xy[:-2].reshape((xy.size // 2 - 1, 2)), layer, datatype)

    def _create_path(self, **kwargs):
        xy = kwargs.pop('xy')
        if 'bgnextn' in kwargs or 'endextn' in kwargs:
            d = kwargs.pop('bgnextn', 0)
            v = xy[0:2] - xy[2:4]
            xy[0:2] = xy[0:2] + d * v / numpy.sqrt(v[0]**2 + v[1]**2)
            d = kwargs.pop('endextn', 0)
            v = xy[-2:] - xy[-4:-2]
            xy[-2:] = xy[-2:] + d * v / numpy.sqrt(v[0]**2 + v[1]**2)
        kwargs['points'] = xy.reshape((xy.size // 2, 2))
        return PolyPath(**kwargs)

    def _create_label(self, xy, width=None, ends=None, **kwargs):
        kwargs['position'] = xy
        return Label(**kwargs)

    def _create_reference(self, **kwargs):
        kwargs['origin'] = kwargs.pop('xy')
        kwargs['ignore_missing'] = True
        ref = CellReference(**kwargs)
        ref.ref_cell = kwargs['ref_cell']
        self._references.append(ref)
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
        kwargs['ignore_missing'] = True
        ref = CellArray(**kwargs)
        ref.ref_cell = kwargs['ref_cell']
        self._references.append(ref)
        return ref

    def extract(self, cell, overwrite_duplicate=False):
        """
        Extract a cell from the this GDSII file and include it in the
        current global library, including referenced dependencies.

        Parameters
        ----------
        cell : ``Cell`` or string
            Cell or name of the cell to be extracted from the imported
            file.  Referenced cells will be automatically extracted as
            well.
        overwrite_duplicate : bool
            If True an existing cell with the same name in the current
            global library will be overwritten.

        Returns
        -------
        out : ``Cell``
            The extracted cell.

        Notes
        -----
        ``CellReference`` or ``CellArray`` instances that referred to an
        overwritten cell are not automatically updated.
        """
        cell = self.cell_dict.get(cell, cell)
        current_library.add(cell, overwrite_duplicate=overwrite_duplicate)
        current_library.add(cell.get_dependencies(True), overwrite_duplicate=overwrite_duplicate)
        return cell

    def top_level(self):
        """
        Output the top level cells from the GDSII data.

        Top level cells are those that are not referenced by any other
        cells.

        Returns
        -------
        out : list
            List of top level cells.
        """
        top = list(self.cell_dict.values())
        for cell in self.cell_dict.values():
            for dependency in cell.get_dependencies():
                if dependency in top:
                    top.remove(dependency)
        return top


class GdsWriter(object):
    """
    GDSII strem library writer.

    The dimensions actually written on the GDSII file will be the
    dimensions of the objects created times the ratio
    ``unit/precision``. For example, if a circle with radius 1.5 is
    created and we set ``unit=1.0e-6`` (1 um) and ``precision=1.0e-9``
    (1 nm), the radius of the circle will be 1.5 um and the GDSII file
    will contain the dimension 1500 nm.

    Parameters
    ----------
    outfile : file or string
        The file (or path) where the GDSII stream will be written.  It
        must be opened for writing operations in binary format.
    name : string
        Name of the GDSII library (file).
    unit : number
        Unit size for the objects in the library (in *meters*).
    precision : number
        Precision for the dimensions of the objects in the library (in
        *meters*).
    timestamp : datetime object
        Sets the GDSII timestamp.  If ``None``, the current time is
        used.

    Notes
    -----

    This class can be used for incremental output of the geometry in
    case the complete layout is too large to be kept in memory all at
    once.

    Examples
    --------
    >>> writer = gdspy.GdsWriter('out-file.gds', unit=1.0e-6,
    ...                          precision=1.0e-9)
    >>> for i in range(10):
    ...     cell = gdspy.Cell('C{}'.format(i), True)
    ...     # Add the contents of this cell...
    ...     writer.write_cell(cell)
    ...     # Clear the memory: erase Cell objects and any other objects
    ...     # that won't be needed.
    ...     del cell
    >>> writer.close()
    """
    __slots__ = '_outfile', '_close', '_res'

    def __init__(self, outfile, name='library', unit=1.0e-6, precision=1.0e-9, timestamp=None):
        if isinstance(outfile, basestring):
            self._outfile = open(outfile, 'wb')
            self._close = True
        else:
            self._outfile = outfile
            self._close = False
        self._res = unit / precision
        now = datetime.datetime.today() if timestamp is None else timestamp
        if len(name) % 2 != 0:
            name = name + '\0'
        self._outfile.write(struct.pack('>5H12h2H', 6, 0x0002, 0x0258, 28, 0x0102, now.year, now.month,
                                        now.day, now.hour, now.minute, now.second, now.year, now.month,
                                        now.day, now.hour, now.minute, now.second, 4 + len(name), 0x0206)
                            + name.encode('ascii') + struct.pack('>2H', 20, 0x0305)
                            + _eight_byte_real(precision / unit) + _eight_byte_real(precision))

    def write_cell(self, cell, timestamp=None):
        """
        Write the specified cell to the file.

        Parameters
        ----------
        cell : ``Cell``
            Cell to be written.
        timestamp : datetime object
            Sets the GDSII timestamp.  If ``None``, the current time is
            used.

        Notes
        -----
        Only the specified cell is written.  Dependencies must be
        manually included.

        Returns
        -------
        out : ``GdsWriter``
            This object.
        """
        self._outfile.write(cell.to_gds(self._res, timestamp))
        return self

    def write_binary_cells(self, binary_cells):
        """
        Write the specified binary cells to the file.

        Parameters
        ----------
        binary_cells : iterable of bytes
            Iterable with binary data for GDSII cells (from
            ``get_binary_cells``, for example).

        Returns
        -------
        out : ``GdsWriter``
            This object.
        """
        for bc in binary_cells:
            self._outfile.write(bc)
        return self

    def close(self):
        """
        Finalize the GDSII stream library.
        """
        self._outfile.write(struct.pack('>2H', 4, 0x0400))
        if self._close:
            self._outfile.close()


def _raw_record_reader(stream):
    """
    Iterator over complete records from a GDSII stream file.

    Parameters
    ----------
    stream : file
        GDSII stream file to be read.

    Returns
    -------
    out : 2-tuple
        Record type and binary data (including header)
    """
    while True:
        header = stream.read(4)
        if len(header) < 4:
            return
        size, rec_type = struct.unpack('>HH', header)
        rec_type = rec_type // 256
        yield (rec_type, header + stream.read(size - 4))


def get_gds_units(infile):
    """
    Return the unit and precision used in the GDS stream file.

    Parameters
    ----------
    infile : file or string
        GDSII stream file to be queried.

    Returns
    -------
    out : 2-tuple
        Return ``(unit, precision)`` from the file.
    """
    if isinstance(infile, basestring):
        infile = open(infile, 'rb')
        close = True
    else:
        close = False
    unit = precision = None
    for rec_type, data in _raw_record_reader(infile):
        # UNITS
        if rec_type == 0x03:
            db_user = _eight_byte_real_to_float(data[4:12])
            db_meters = _eight_byte_real_to_float(data[12:])
            unit = db_meters / db_user
            precision = db_meters
            break
    if close:
        infile.close()
    return (unit, precision)


def get_binary_cells(infile):
    """
    Load all cells from a GDSII stream file in binary format.

    Parameters
    ----------
    infile : file or string
        GDSII stream file (or path) to be loaded.  It must be opened for
        reading in binary format.

    Returns
    -------
    out : dictionary
        Dictionary of binary cell representations indexed by name.

    Notes
    -----
    The returned cells inherit the units of the loaded file.  If they
    are used in a new library, the new library must use compatible
    units.
    """
    if isinstance(infile, basestring):
        infile = open(infile, 'rb')
        close = True
    else:
        close = False
    cells = {}
    name = None
    cell_data = None
    for rec_type, data in _raw_record_reader(infile):
        # BGNSTR
        if rec_type == 0x05:
            cell_data = [data]
        # STRNAME
        elif rec_type == 0x06:
            cell_data.append(data)
            if str is not bytes:
                if data[-1] == 0:
                    name = data[4:-1].decode('ascii')
                else:
                    name = data[4:].decode('ascii')
            else:
                if data[-1] == '\0':
                    name = data[4:-1]
                else:
                    name = data[4:]
        # ENDSTR
        elif rec_type == 0x07:
            cell_data.append(data)
            cells[name] = b''.join(cell_data)
            cell_data = None
        elif cell_data is not None:
            cell_data.append(data)
    if close:
        infile.close()
    return cells


def _gather_polys(args):
    """
    Gather polygons from different argument types into a list.

    Parameters
    ----------
    args : ``None``, ``PolygonSet``, ``CellReference``, ``CellArray`` or iterable
        Polygon types.  If this is an iterable, each element must be a
        ``PolygonSet``, ``CellReference``, ``CellArray``, or an
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
    if isinstance(args, CellReference) or isinstance(args, CellArray):
        return args.get_polygons()
    polys = []
    for p in args:
        if isinstance(p, PolygonSet):
            polys.extend(p.polygons)
        elif isinstance(p, CellReference) or isinstance(p, CellArray):
            polys.extend(p.get_polygons())
        else:
            polys.append(p)
    return polys


def slice(polygons, position, axis, precision=1e-3, layer=0, datatype=0):
    """
    Slice polygons and polygon sets at given positions along an axis.

    Parameters
    ----------
    polygons : ``PolygonSet``, ``CellReference``, ``CellArray`` or iterable
        Operand of the slice operation.  If this is an iterable, each
        element must be a ``PolygonSet``, ``CellReference``,
        ``CellArray``, or an array-like[N][2] of vertices of a polygon.
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
    out : list[N] of ``PolygonSet`` or ``None``
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
            result[i] = PolygonSet(result[i], layer[i % len(layer)], datatype[i % len(datatype)])
    return result


def offset(polygons, distance, join='miter', tolerance=2, precision=0.001, join_first=False,
           max_points=199, layer=0, datatype=0):
    """
    Shrink or expand a polygon or polygon set.

    Parameters
    ----------
    polygons : ``PolygonSet``, ``CellReference``, ``CellArray`` or iterable
        Polygons to be offset.  If this is an iterable, each element
        must be a ``PolygonSet``, ``CellReference``, ``CellArray``, or
        an array-like[N][2] of vertices of a polygon.
    distance : number
        Offset distance.  Positive to expand, negative to shrink.
    join : 'miter', 'bevel', 'round'
        Type of join used to create the offset polygon.
    tolerance : number
        For miter joints, this number must be at least 2 and it
        represents the maximun distance in multiples of offset betwen
        new vertices and their original position before beveling to
        avoid spikes at acute joints.  For round joints, it indicates
        the curvature resolution in number of points per full circle.
    precision : float
        Desired precision for rounding vertice coordinates.
    join_first : bool
        Join all paths before offseting to avoid unecessary joins in
        adjacent polygon sides.
    max_points : integer
        If greater than 4, fracture the resulting polygons to ensure
        they have at most ``max_points`` vertices.  This is not a
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
    out : ``PolygonSet`` or ``None``
        Return the offset shape as a set of polygons.
    """
    result = clipper.offset(_gather_polys(polygons), distance, join, tolerance, 1 / precision,
                            1 if join_first else 0)
    if len(result) == 0:
        return None
    return PolygonSet(result, layer, datatype).fracture(max_points, precision)


def boolean(operand1, operand2, operation, precision=0.001, max_points=199, layer=0, datatype=0):
    """
    Execute any boolean operation between 2 polygons or polygon sets.

    Parameters
    ----------
    operand1 : ``PolygonSet``, ``CellReference``, ``CellArray`` or iterable
        First operand.  If this is an iterable, each element must be a
        ``PolygonSet``, ``CellReference``, ``CellArray``, or an
        array-like[N][2] of vertices of a polygon.
    operand2 : ``None``, ``PolygonSet``, ``CellReference``, ``CellArray`` or iterable
        Second operand.  If this is an iterable, each element must be a
        ``PolygonSet``, ``CellReference``, ``CellArray``, or an
        array-like[N][2] of vertices of a polygon.
    operation : {'or', 'and', 'xor', 'not'}
        Boolean operation to be executed.  The 'not' operation returns
        the difference ``operand1 - operand2``.
    precision : float
        Desired precision for rounding vertice coordinates.
    max_points : integer
        If greater than 4, fracture the resulting polygons to ensure
        they have at most ``max_points`` vertices.  This is not a
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
    out : PolygonSet or ``None``
        Result of the boolean operation.
    """
    poly1 = _gather_polys(operand1)
    poly2 = _gather_polys(operand2)
    if len(poly2) == 0:
        poly2.append(poly1.pop())
    result = clipper.clip(poly1, poly2, operation, 1 / precision)
    if len(result) == 0:
        return None
    return PolygonSet(result, layer, datatype).fracture(max_points, precision)


fast_boolean = boolean


def inside(points, polygons, short_circuit='any', precision=0.001):
    """
    Test whether each of the points is within the given set of polygons.

    Parameters
    ----------
    points : array-like[N][2] or sequence of array-like[N][2]
        Coordinates of the points to be tested or groups of points to be
        tested together.
    polygons : ``PolygonSet``, ``CellReference``, ``CellArray`` or iterable
        Polygons to be tested against.  If this is an iterable, each
        element must be a ``PolygonSet``, ``CellReference``,
        ``CellArray``, or an array-like[N][2] of vertices of a polygon.
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
    if hasattr(points[0][0], '__iter__'):
        pts = points
        sc = 1 if short_circuit == 'any' else -1
    else:
        pts = (points, )
        sc = 0
    return clipper.inside(pts, polys, sc, 1 / precision)


def copy(obj, dx, dy):
    """
    Creates a copy of ``obj`` and translates the new object to a new
    location.

    Parameters
    ----------
    obj : ``obj``
        any translatable geometery object.
    dx  : float
        distance to move in the x-direction
    dy  : float
        distance to move in the y-direction


    Returns
    -------
    out : ``obj``
        Translated copy of original ``obj``

    Examples
    --------
    >>> rectangle = gdspy.Rectangle((0, 0), (10, 20))
    >>> rectangle2 = gdspy.copy(rectangle, 2,0)
    >>> myCell.add(rectangle)
    >>> myCell.add(rectangle2)
    """

    newObj = libCopy.deepcopy(obj)
    newObj.translate(dx, dy)
    return newObj


def write_gds(outfile, cells=None, name='library', unit=1.0e-6, precision=1.0e-9, timestamp=None,
              binary_cells=None):
    """
    Write the current GDSII library to a file.

    The dimensions actually written on the GDSII file will be the
    dimensions of the objects created times the ratio
    ``unit/precision``.  For example, if a circle with radius 1.5 is
    created and we set ``unit=1.0e-6`` (1 um) and ``precision=1.0e-9``
    (1 nm), the radius of the circle will be 1.5 um and the GDSII file
    will contain the dimension 1500 nm.

    Parameters
    ----------
    outfile : file or string
        The file (or path) where the GDSII stream will be written.  It
        must be opened for writing operations in binary format.
    cells : array-like
        The sequence of cells or cell names to be included in the
        library.  If ``None``, all cells are used.
    name : string
        Name of the GDSII library.
    unit : number
        Unit size for the objects in the library (in *meters*).
    precision : number
        Precision for the dimensions of the objects in the library (in
        *meters*).
    timestamp : datetime object
        Sets the GDSII timestamp.  If ``None``, the current time is
        used.
    binary_cells : iterable of bytes
        Iterable with binary data for GDSII cells (from
        ``get_binary_cells``, for example).
    """
    current_library.name = name
    current_library.unit = unit
    current_library.precision = precision
    current_library.write_gds(outfile, cells, timestamp, binary_cells)


def gdsii_hash(filename, engine=None):
    """
    Calculate the a hash value for a GDSII file.

    The hash is generated based only on the contents of the cells in the
    GDSII library, ignoring any timestamp records present in the file
    structure.

    Parameters
    ----------
    filename : string
        Full path to the GDSII file.
    engine : hashlib-like engine
        The engine that executes the hashing algorithm.  It must provide
        the methods ``update`` and ``hexdigest`` as defined in the
        hashlib module.  If ``None``, the dafault ``hashlib.sha1()`` is
        used.

    Returns
    -------
    out : string
        The hash correponding to the library contents in hex format.
    """
    with open(filename, 'rb') as fin:
        data = fin.read()
    contents = []
    start = pos = 0
    while pos < len(data):
        size, rec = struct.unpack('>HH', data[pos:pos + 4])
        if rec == 0x0502:
            start = pos + 28
        elif rec == 0x0700:
            contents.append(data[start:pos])
        pos += size
    h = hashlib.sha1() if engine is None else engine
    for x in sorted(contents):
        h.update(x)
    return h.hexdigest()


current_library = GdsLibrary()
"""
Current ``GdsLibrary`` instance for automatic creation of GDSII files.

This variable can be freely overwritten by the user with a new instance
of ``GdsLibrary``.
"""

poly_warnings = True
"""
Flag controlling the emission of warnings relative to polygon creation.
"""
