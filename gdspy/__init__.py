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
import copy as libcopy
import hashlib

from gdspy import clipper
try:
    from gdspy.viewer import LayoutViewer
except ImportError as e:
    warnings.warn(
        "[GDSPY] LayoutViewer not available: " + str(e),
        category=ImportWarning,
        stacklevel=2)

__version__ = '1.4'

_halfpi = 0.5 * numpy.pi
_zero = numpy.array((0.0, 0.0))
_one = numpy.array((1.0, 1.0))
_mpone = numpy.array((-1.0, 1.0))
_pmone = numpy.array((1.0, -1.0))
_pmone_int = numpy.array((1, -1))

_directions_dict = {'+x': 0, '+y': 0.5, '-x': 1, '-y': -0.5}
_directions_list = ['+x', '+y', '-x', '-y']
_angle_dic = {'l': _halfpi, 'r': -_halfpi, 'll': numpy.pi, 'rr': -numpy.pi}

_bounding_boxes = {}


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
        The GDSII binary string that represents `value`.
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
        The number represented by `value`.
    """
    short1, short2, long3 = struct.unpack('>HHL', value)
    exponent = (short1 & 0x7f00) // 256 - 64
    mantissa = (((short1 & 0x00ff) * 65536 + short2) * 4294967296 + long3) / 72057594037927936.0
    if short1 & 0x8000:
        return -mantissa * 16.**exponent
    return mantissa * 16.**exponent


def _hobby(points, angles=None, curl_start=1, curl_end=1, t_in=1, t_out=1, cycle=False):
    """
    Calculate control points for a smooth interpolating curve.

    Uses the Hobby algorithm [1]_ to calculate a smooth interpolating
    curve made of cubic Bezier segments between each pair of points.

    Parameters
    ----------
    points : Numpy array[N, 2]
        Vertices in the interpolating curve.
    angles : array-like[N] or None
        Tangent angles at each point (in *radians*).  Any angles defined
        as None are automatically calculated.
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
        If True, calculates control points for a closed curve, with
        an additional segment connecting the first and last points.

    Returns
    -------
    out : 2-tuple of Numpy array[M, 2]
        Pair of control points for each segment in the interpolating
        curve.  For a closed curve (`cycle` True), M = N.  For an open
        curve (`cycle` False), M = N - 1.

    References
    ----------
    .. [1] Hobby, J.D.  *Discrete Comput. Geom.* (1986) 1: 123.
       `DOI: 10.1007/BF02187690 <https://doi.org/10.1007/BF02187690>`_
    """
    z = points[:, 0] + 1j * points[:, 1]
    n = z.size
    if numpy.isscalar(t_in):
        t_in = t_in * numpy.ones(n)
    else:
        t_in = numpy.array(t_in)
    if numpy.isscalar(t_out):
        t_out = t_out * numpy.ones(n)
    else:
        t_out = numpy.array(t_out)
    if angles is None:
        angles = [None] * n
    rotate = 0
    if cycle and any(a is not None for a in angles):
        while angles[rotate] is None:
            rotate += 1
        angles = [angles[(rotate + j) % n] for j in range(n + 1)]
        z = numpy.hstack((numpy.roll(z, -rotate), z[rotate:rotate + 1]))
        t_in = numpy.hstack((numpy.roll(t_in, -rotate), t_in[rotate:rotate + 1]))
        t_out = numpy.hstack((numpy.roll(t_out, -rotate), t_out[rotate:rotate + 1]))
        cycle = False
    if cycle:
        # Closed curve
        v = numpy.roll(z, -1) - z
        d = numpy.abs(v)
        delta = numpy.angle(v)
        psi = (delta - numpy.roll(delta, 1) + numpy.pi) % (2 * numpy.pi) - numpy.pi
        coef = numpy.zeros(2 * n)
        coef[:n] = -psi
        m = numpy.zeros((2 * n, 2 * n))
        i = numpy.arange(n)
        i1 = (i + 1) % n
        i2 = (i + 2) % n
        ni = n + i
        m[i, i] = 1
        m[i, n + (i - 1) % n] = 1
        # A_i
        m[ni, i] = d[i1] * t_in[i2] * t_in[i1]**2
        # B_{i+1}
        m[ni, i1] = -d[i] * t_out[i] * t_out[i1]**2 * (1 - 3 * t_in[i2])
        # C_{i+1}
        m[ni, ni] = d[i1] * t_in[i2] * t_in[i1]**2 * (1 - 3 * t_out[i])
        # D_{i+2}
        m[ni, n + i1] = -d[i] * t_out[i] * t_out[i1]**2
        sol = numpy.linalg.solve(m, coef)
        theta = sol[:n]
        phi = sol[n:]
        w = numpy.exp(1j * (theta + delta))
        a = 2**0.5
        b = 1.0 / 16
        c = (3 - 5**0.5) / 2
        sintheta = numpy.sin(theta)
        costheta = numpy.cos(theta)
        sinphi = numpy.sin(phi)
        cosphi = numpy.cos(phi)
        alpha = a * (sintheta - b * sinphi) * (sinphi - b * sintheta) * (costheta - cosphi)
        cta = z + w * d * ((2 + alpha) / (1 + (1 - c) * costheta + c * cosphi)) / (3 * t_out)
        ctb = numpy.roll(z, -1) - numpy.roll(w, -1) * d * (
            (2 - alpha) / (1 + (1 - c) * cosphi + c * costheta)) / (3 * numpy.roll(t_in, -1))
    else:
        # Open curve(s)
        n = z.size - 1
        v = z[1:] - z[:-1]
        d = numpy.abs(v)
        delta = numpy.angle(v)
        psi = (delta[1:] - delta[:-1] + numpy.pi) % (2 * numpy.pi) - numpy.pi
        theta = numpy.empty(n)
        phi = numpy.empty(n)
        i = 0
        if angles[0] is not None:
            theta[0] = angles[0] - delta[0]
        while i < n:
            j = i + 1
            while j < n + 1 and angles[j] is None:
                j += 1
            if j == n + 1:
                j -= 1
            else:
                phi[j - 1] = delta[j - 1] - angles[j]
                if j < n:
                    theta[j] = angles[j] - delta[j]
            # Solve open curve z_i, ..., z_j
            nn = j - i
            coef = numpy.zeros(2 * nn)
            coef[1:nn] = -psi[i:j - 1]
            m = numpy.zeros((2 * nn, 2 * nn))
            if nn > 1:
                ii = numpy.arange(nn - 1) #[0 .. nn-2]
                i0 = i + ii               #[i .. j-1]
                i1 = 1 + i0               #[i+1 .. j]
                i2 = 2 + i0               #[i+2 .. j+1]
                ni = nn + ii              #[nn .. 2*nn-2]
                ii1 = 1 + ii              #[1 .. nn-1]
                m[ii1, ii1] = 1
                m[ii1, ni] = 1
                # A_ii
                m[ni, ii] = d[i1] * t_in[i2] * t_in[i1]**2
                # B_{ii+1}
                m[ni, ii1] = -d[i0] * t_out[i0] * t_out[i1]**2 * (1 - 3 * t_in[i2])
                # C_{ii+1}
                m[ni, ni] = d[i1] * t_in[i2] * t_in[i1]**2 * (1 - 3 * t_out[i0])
                # D_{ii+2}
                m[ni, ni + 1] = -d[i0] * t_out[i0] * t_out[i1]**2
            if angles[i] is None:
                to3 = t_out[0]**3
                cti3 = curl_start * t_in[1]**3
                # B_0
                m[0, 0] = to3 * (1 - 3 * t_in[1]) - cti3
                # D_1
                m[0, nn] = to3 - cti3 * (1 - 3 * t_out[0])
            else:
                coef[0] = theta[i]
                m[0, 0] = 1
                m[0, nn] = 0
            if angles[j] is None:
                ti3 = t_in[n]**3
                cto3 = curl_end * t_out[n - 1]**3
                # A_{nn-1}
                m[2 * nn - 1, nn - 1] = ti3 - cto3 * (1 - 3 * t_in[n])
                # C_nn
                m[2 * nn - 1, 2 * nn - 1] = ti3 * (1 - 3 * t_out[n - 1]) - cto3
            else:
                coef[2 * nn - 1] = phi[j - 1]
                m[2 * nn - 1, nn - 1] = 0
                m[2 * nn - 1, 2 * nn - 1] = 1
            if nn > 1 or angles[i] is None or angles[j] is None:
                sol = numpy.linalg.solve(m, coef)
                theta[i:j] = sol[:nn]
                phi[i:j] = sol[nn:]
            i = j
        w = numpy.hstack((numpy.exp(1j * (delta + theta)),
                          numpy.exp(1j * (delta[-1:] - phi[-1:]))))
        a = 2**0.5
        b = 1.0 / 16
        c = (3 - 5**0.5) / 2
        sintheta = numpy.sin(theta)
        costheta = numpy.cos(theta)
        sinphi = numpy.sin(phi)
        cosphi = numpy.cos(phi)
        alpha = a * (sintheta - b * sinphi) * (sinphi - b * sintheta) * (costheta - cosphi)
        cta = z[:-1] + w[:-1] * d * (
            (2 + alpha) / (1 + (1 - c) * costheta + c * cosphi)) / (3 * t_out[:-1])
        ctb = z[1:] - w[1:] * d * (
            (2 - alpha) / (1 + (1 - c) * cosphi + c * costheta)) / (3 * t_in[1:])
        if rotate > 0:
            cta = numpy.roll(cta, rotate)
            ctb = numpy.roll(ctb, rotate)
    return (numpy.vstack((cta.real, cta.imag)).transpose(), numpy.vstack((ctb.real, ctb.imag)).transpose())


def _func_const(c, nargs=1):
    if nargs == 1:
        return lambda u: c
    elif nargs == 2:
        return lambda u, h: c
    return lambda *args: c


def _func_linear(c0, c1):
    return lambda u: c0 * (1 - u) + c1 * u


def _func_multadd(f, m=None, a=None, nargs=1):
    if nargs == 1:
        if m is None:
            return lambda u: f(u) + a
        if a is None:
            return lambda u: f(u) * m
        return lambda u: f(u) * m + a
    elif nargs == 2:
        if m is None:
            return lambda u, h: f(u, h) + a
        if a is None:
            return lambda u, h: f(u, h) * m
        return lambda u, h: f(u, h) * m + a
    if m is None:
        return lambda *args: f(*args) + a
    if a is None:
        return lambda *args: f(*args) * m
    return lambda *args: f(*args) * m + a


def _func_rotate(f, cos, sin, center=0, nargs=1):
    if nargs == 1:
        def _f(u):
            x = f(u) - center
            return x * cos + x[::-1] * sin + center
    elif nargs == 2:
        def _f(u, h):
            x = f(u, h) - center
            return x * cos + x[::-1] * sin + center
    return _f


def _func_trafo(f, translation, rotation, scale, x_reflection, array_trans, nargs=1):
    if translation is None:
        translation = _zero
    if array_trans is None:
        array_trans = _zero
    if rotation is None:
        cos = _one
        sin = _zero
    else:
        cos = numpy.cos(rotation) * _one
        sin = numpy.sin(rotation) * _mpone
    if scale is not None:
        cos = cos * scale
        sin = sin * scale
        array_trans = array_trans / scale
    if x_reflection:
        cos[1] = -cos[1]
        sin[0] = -sin[0]
    if nargs == 1:
        def _f(u):
            x = f(u) + array_trans
            return x * cos + x[::-1] * sin + translation
    elif nargs == 2:
        def _f(u, h):
            x = f(u, h) + array_trans
            return x * cos + x[::-1] * sin + translation
    return _f


def _func_bezier(ctrl, nargs=1):
    if nargs == 1:
        def _f(u):
            p = ctrl
            for _ in range(ctrl.shape[0] - 1):
                p = p[1:] * u + p[:-1] * (1 - u)
            return p[0]
    elif nargs == 2:
        def _f(u, h):
            p = ctrl
            for _ in range(ctrl.shape[0] - 1):
                p = p[1:] * u + p[:-1] * (1 - u)
            return p[0]
    return _f
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
    if isinstance(args, RobustPath) or isinstance(args, FlexPath) or isinstance(args, CellReference) or isinstance(args, CellArray):
        return args.get_polygons()
    polys = []
    for p in args:
        if isinstance(p, PolygonSet):
            polys.extend(p.polygons)
        elif isinstance(p, RobustPath) or isinstance(args, FlexPath) or isinstance(p, CellReference) or isinstance(p, CellArray):
            polys.extend(p.get_polygons())
        else:
            polys.append(p)
    return polys


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

    The original GDSII specification supports only a maximum of 199
    vertices per polygon.
    """

    __slots__ = 'layers', 'datatypes', 'polygons'

    def __init__(self, polygons, layer=0, datatype=0):
        self.polygons = [numpy.array(p) for p in polygons]
        self.layers = [layer] * len(self.polygons)
        self.datatypes = [datatype] * len(self.polygons)

    def __str__(self):
        return ("PolygonSet ({} polygons, {} vertices, layers {}, datatypes {})").format(len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))

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
            if len(self.polygons[ii]) > 8190:
                warnings.warn("[GDSPY] Polygons with more than 8190 are not supported by the official GDSII specification.  This GDSII file might not be compatible with all readers.", stacklevel=4)
                data.append(struct.pack('>4Hh2Hh', 4, 0x0800, 6, 0x0D02, self.layers[ii], 6, 0x0E02,
                                        self.datatypes[ii]))
                xy = numpy.empty((self.polygons[ii].shape[0] + 1, 2), dtype='>i4')
                xy[:-1] = numpy.round(self.polygons[ii] * multiplier)
                xy[-1] = xy[0]
                i0 = 0
                while i0 < xy.shape[0]:
                    i1 = min(i0 + 8190, xy.shape[0])
                    data.append(struct.pack('>2H', 4 + 8 * (i1 - i0), 0x1003))
                    data.append(xy[i0:i1].tostring())
                    i0 = i1
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
                            raise ValueError("[GDSPY] Wrong length in fillet radius list.  Expected lengths are {} or {}; got {}.".format(len(self.polygons), total, len(radius)))
                        radii.append(r)
            else:
                total = sum(p.shape[0] for p in self.polygons)
                if len(radius) != total:
                    raise ValueError("[GDSPY] Wrong length in fillet radius list.  Expected lengths are {} or {}; got {}.".format(len(self.polygons), total, len(radius)))
                radii = []
                n = 0
                for p in self.polygons:
                    radii.append(radius[n:n + p.shape[0]])
                    n += p.shape[0]

        for jj in range(len(self.polygons)):
            vec = self.polygons[jj].astype(float) - numpy.roll(self.polygons[jj], 1, 0)
            length = (vec[:, 0]**2 + vec[:, 1]**2)**0.5
            ii = numpy.flatnonzero(length)
            if len(ii) < len(length):
                self.polygons[jj] = numpy.array(self.polygons[jj][ii])
                radii[jj] = [radii[jj][i] for i in ii]
                vec = self.polygons[jj].astype(float) - numpy.roll(self.polygons[jj], 1, 0)
                length = (vec[:, 0]**2 + vec[:, 1]**2)**0.5
            vec[:, 0] = vec[:, 0] / length
            vec[:, 1] = vec[:, 1] / length
            dvec = numpy.roll(vec, -1, 0) - vec
            norm = (dvec[:, 0]**2 + dvec[:, 1]**2)**0.5
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

    The original GDSII specification supports only a maximum of 199
    vertices per polygon.

    Examples
    --------
    >>> triangle_pts = [(0, 40), (15, 40), (10, 50)]
    >>> triangle = gdspy.Polygon(triangle_pts)
    >>> myCell.add(triangle)
    """

    __slots__ = 'layers', 'datatypes', 'polygons'

    def __init__(self, points, layer=0, datatype=0):
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
            full_angle = 2 * numpy.pi if final_angle == initial_angle else abs(final_angle - initial_angle)
            number_of_points = max(3, 1 + int(0.5 * full_angle / numpy.arccos(1 - tolerance / radius) + 0.5))
            if inner_radius > 0:
                number_of_points *= 2

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
        self.polygons = [(points - c0) * ca + (points - c0)[:, ::-1] * sa + c0 for points in self.polygons]
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
        self.polygons = [numpy.outer(numpy.inner(points - origin, vec_r), vec) - points + 2 * origin
                         for points in self.polygons]
        dot = (self.x - origin[0]) * vec_r[0] + (self.y - origin[1]) * vec_r[1]
        self.x = dot * vec[0] - self.x + 2 * origin[0]
        self.y = dot * vec[1] - self.y + 2 * origin[1]
        if isinstance(self.direction, basestring):
            self.direction = _directions_dict[self.direction] * numpy.pi
        self.direction = 2 * numpy.arctan2(vec[1], vec[0]) - self.direction
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
        out : `Path`
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
                    self.polygons[-1] = self.polygons[-1][:-1]
                if old_w == 0:
                    self.polygons[-1] = self.polygons[-1][1:]
            self.length += (length**2 + axis_offset**2)**0.5
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
            warnings.warn("[GDSPY] Use of a floating number as number_of_points is deprecated in favor of tolerance.", category=DeprecationWarning, stacklevel=2)
            tolerance = number_of_points
            number_of_points = None
        if number_of_points is None:
            r = radius + max(old_distance, self.distance) * (self.n - 1) * 0.5 + max(old_w, self.w)
            number_of_points = max(6, 2 + 2 * int(0.5 * abs(final_angle - initial_angle) / numpy.arccos(1 - tolerance / r) + 0.5))
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
                        self.polygons[-1][:pts1 - 1] = numpy.array(self.polygons[-1][1:pts1])
                        pts1 -= 1
                        pts2 += 1
                    ang = numpy.linspace(angles[jj + 1], angles[jj], pts2)
                    rad = numpy.linspace(r0 - widths[jj + 1], old_r0 - widths[jj], pts2)
                    if rad[0] <= 0 or rad[-1] <= 0:
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

    def parametric(self, curve_function, curve_derivative=None, tolerance=0.01,
                   number_of_evaluations=5, max_points=199, final_width=None, final_distance=None,
                   relative=True, layer=0, datatype=0):
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
            if test_err[0]**2 + test_err[1]**2 > err:
                delta = min(delta, points[i] - midpoint)
                points.insert(i, midpoint)
                values.insert(i, midvalue)
            else:
                i += 1
        points = numpy.array(points)
        values = numpy.array(values)
        dvs = values[1:] - values[:-1]
        self.length += ((dvs[:, 0]**2 + dvs[:, 1]**2)**0.5).sum()

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
        if relative:
            x0 = values + numpy.array((self.x, self.y))
        else:
            x0 = values
        dx = derivs[:, ::-1] * _mpone / ((derivs[:, 0]**2 + derivs[:, 1]**2)**0.5).reshape(sh)
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
                p2 = (x0[i0:i1] + dx[i0:i1] * (dist[i0:i1] * (ii - (self.n - 1) * 0.5) - width[i0:i1]))[::-1]
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
               final_width=None, final_distance=None, relative=True, layer=0, datatype=0):
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
        self.parametric(_func_bezier(pts), _func_bezier(dpts), tolerance, number_of_evaluations,
                        max_points, final_width, final_distance, relative, layer, datatype)
        return self

    def smooth(self, points, angles=None, curl_start=1, curl_end=1, t_in=1, t_out=1, cycle=False,
               tolerance=0.01, number_of_evaluations=5, max_points=199, final_widths=None,
               final_distances=None, relative=True, layer=0, datatype=0):
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
            points = numpy.vstack(([(0.0, 0.0)], points)) + numpy.array((self.x, self.y))
        else:
            points = numpy.vstack(([(self.x, self.y)], points))
        cta, ctb = _hobby(points, angles, curl_start, curl_end, t_in, t_out, cycle)
        if final_widths is None:
            final_widths = [None] * cta.shape[0]
        if final_distances is None:
            final_distances = [None] * cta.shape[0]
        for i in range(points.shape[0] - 1):
            self.bezier([cta[i], ctb[i], points[i + 1]], tolerance, number_of_evaluations, max_points,
                        final_widths[i], final_distances[i], False, layer, datatype)
        if cycle:
            self.bezier([cta[-1], ctb[-1], points[0]], tolerance, number_of_evaluations, max_points,
                        final_widths[-1], final_distances[-1], False, layer, datatype)
        return self


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
        warnings.warn("[GDSPY] L1Path is deprecated favor of FlexPath and will be removed in a future version of Gdspy.", category=DeprecationWarning, stacklevel=2)
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
        self.polygons = [(points - c0) * ca + (points - c0)[:, ::-1] * sa + c0 for points in self.polygons]
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
    __slots__ = 'layers', 'datatypes', 'polygons'

    def __init__(self, points, width, number_of_paths=1, distance=0, corners='miter', ends='flush',
                 max_points=199, layer=0, datatype=0):
        warnings.warn("[GDSPY] PolyPath is deprecated favor of FlexPath and will be removed in a future version of Gdspy.", category=DeprecationWarning, stacklevel=2)
        if not isinstance(layer, list):
            layer = [layer]
        if not isinstance(datatype, list):
            datatype = [datatype]
        if hasattr(width, '__iter__'):
            width = numpy.array(width) * 0.5
        else:
            width = numpy.array([width * 0.5])
        len_w = len(width)
        if hasattr(distance, '__iter__'):
            distance = numpy.array(distance)
        else:
            distance = numpy.array([distance])
        len_d = len(distance)
        points = numpy.array(points, dtype=float)
        self.polygons = []
        self.layers = []
        self.datatypes = []
        if points.shape[0] == 2 and number_of_paths == 1:
            v = points[1] - points[0]
            v = v / (v[0]**2 + v[1]**2)**0.5
            w0 = width[0]
            w1 = width[1 % len_w]
            if ends == 'round':
                a = numpy.arctan2(v[1], v[0]) + _halfpi
                self.polygons.append(Round(points[0], w0, initial_angle=a,
                                           final_angle=a + numpy.pi, number_of_points=33).polygons[0])
                self.polygons.append(Round(points[1], w1, initial_angle=a - numpy.pi,
                                           final_angle=a, number_of_points=33).polygons[0])
                self.layers.extend(layer[:1] * 2)
                self.datatypes.extend(datatype[:1] * 2)
            if ends == 'extended':
                points[0] = points[0] - v * w0
                points[1] = points[1] + v * w1
            u = numpy.array((-v[1], v[0]))
            if w0 == 0:
                self.polygons.append(numpy.array((points[0], points[1] - u * w1, points[1] + u * w1)))
            elif w1 == 0:
                self.polygons.append(numpy.array((points[0] + u * w0, points[0] - u * w0, points[1])))
            else:
                self.polygons.append(numpy.array((points[0] + u * w0, points[0] - u * w0, points[1] - u * w1, points[1] + u * w1)))
            self.layers.append(layer[0])
            self.datatypes.append(datatype[0])
            return
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
            v = points[0] - points[1]
            v = v / (v[0]**2 + v[1]**2)**0.5
            points[0] = points[0] + v * width[0]
            v = points[-1] - points[-2]
            v = v / (v[0]**2 + v[1]**2)**0.5
            points[-1] = points[-1] + v * width[(points.shape[0] - 1) % len_w]
        elif ends == 'round':
            v0 = points[1] - points[0]
            angle0 = numpy.arctan2(v0[1], v0[0]) + _halfpi
            v0 = numpy.array((-v0[1], v0[0])) / (v0[0]**2 + v0[1]**2)**0.5
            d0 = 0.5 * (number_of_paths - 1) * distance[0]
            v1 = points[-1] - points[-2]
            angle1 = numpy.arctan2(v1[1], v1[0]) - _halfpi
            v1 = numpy.array((-v1[1], v1[0])) / (v1[0]**2 + v1[1]**2)**0.5
            j1w = (points.shape[0] - 1) % len_w
            j1d = (points.shape[0] - 1) % len_d
            d1 = 0.5 * (number_of_paths - 1) * distance[j1d]
            self.polygons.extend((Round(points[0] + v0 * (ii * distance[0] - d0), width[0],
                                        initial_angle=angle0, final_angle=angle0 + numpy.pi,
                                        number_of_points=33).polygons[0] for ii in range(number_of_paths)))
            self.polygons.extend((Round(points[-1] + v1 * (ii * distance[j1d] - d1), width[j1w],
                                        initial_angle=angle1, final_angle=angle1 + numpy.pi,
                                        number_of_points=33).polygons[0]) for ii in range(number_of_paths))
            self.layers.extend(((layer * (number_of_paths // len(layer) + 1))[:number_of_paths]) * 2)
            self.datatypes.extend(((datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths]) * 2)
        v = points[1] - points[0]
        v = numpy.array((-v[1], v[0])) / (v[0]**2 + v[1]**2)**0.5
        d0 = 0.5 * (number_of_paths - 1) * distance[0]
        d1 = 0.5 * (number_of_paths - 1) * distance[1 % len_d]
        paths = [[[points[0] + (ii * distance[0] - d0 - width[0]) * v],
                  [points[0] + (ii * distance[0] - d0 + width[0]) * v]] for ii in range(number_of_paths)]
        p1 = [(points[1] + (ii * distance[1 % len_d] - d1 - width[1 % len_w]) * v,
               points[1] + (ii * distance[1 % len_d] - d1 + width[1 % len_w]) * v) for ii in range(number_of_paths)]
        for jj in range(1, points.shape[0] - 1):
            j0d = jj % len_d
            j0w = jj % len_w
            j1d = (jj + 1) % len_d
            j1w = (jj + 1) % len_w
            v = points[jj + 1] - points[jj]
            v = numpy.array((-v[1], v[0])) / (v[0]**2 + v[1]**2)**0.5
            d0 = d1
            d1 = 0.5 * (number_of_paths - 1) * distance[j1d]
            p0 = p1
            p1 = []
            pp = []
            for ii in range(number_of_paths):
                pp.append((points[jj] + (ii * distance[j0d] - d0 - width[j0w]) * v,
                           points[jj] + (ii * distance[j0d] - d0 + width[j0w]) * v))
                p1.append((points[jj + 1] + (ii * distance[j1d] - d1 - width[j1w]) * v,
                           points[jj + 1] + (ii * distance[j1d] - d1 + width[j1w]) * v))
                for kk in (0, 1):
                    p0m = paths[ii][kk][-1] - p0[ii][kk]
                    p1p = pp[ii][kk] - p1[ii][kk]
                    vec = p0m[0] * p1p[1] - p1p[0] * p0m[1]
                    if abs(vec) > 1e-30:
                        p = _pmone * (p0m * p1p[::-1] * p1[ii][kk] - p1p * p0m[::-1] * p0[ii][kk] + p0m * p1p * (p0[ii][kk][::-1] - p1[ii][kk][::-1])) / vec
                        l0 = (p - pp[ii][kk]) * p1p
                        l1 = (p - p0[ii][kk]) * p0m
                        if bevel and l0[0] + l0[1] > 0 and l1[0] + l1[1] < 0:
                            paths[ii][kk].append(p0[ii][kk])
                            paths[ii][kk].append(pp[ii][kk])
                        else:
                            paths[ii][kk].append(p)
                if max_points > 0 and len(paths[ii][0]) + len(paths[ii][1]) + 3 > max_points:
                    diff = paths[ii][0][0] - paths[ii][1][0]
                    if diff[0]**2 + diff[1]**2 == 0:
                        paths[ii][1] = paths[ii][1][1:]
                    diff = paths[ii][0][-1] - paths[ii][1][-1]
                    if diff[0]**2 + diff[1]**2 == 0:
                        self.polygons.append(numpy.array(paths[ii][0] + paths[ii][1][-2::-1]))
                    else:
                        self.polygons.append(numpy.array(paths[ii][0] + paths[ii][1][::-1]))
                    paths[ii][0] = paths[ii][0][-1:]
                    paths[ii][1] = paths[ii][1][-1:]
                    self.layers.append(layer[ii % len(layer)])
                    self.datatypes.append(datatype[ii % len(datatype)])
        for ii in range(number_of_paths):
            diff = paths[ii][0][0] - paths[ii][1][0]
            if diff[0]**2 + diff[1]**2 == 0:
                paths[ii][1] = paths[ii][1][1:]
            diff = p1[ii][0] - p1[ii][1]
            if diff[0]**2 + diff[1]**2 != 0:
                paths[ii][0].append(p1[ii][0])
            paths[ii][1].append(p1[ii][1])
        self.polygons.extend(numpy.array(pol[0] + pol[1][::-1]) for pol in paths)
        self.layers.extend((layer * (number_of_paths // len(layer) + 1))[:number_of_paths])
        self.datatypes.extend((datatype * (number_of_paths // len(datatype) + 1))[:number_of_paths])

    def __str__(self):
        return "PolyPath ({} polygons, {} vertices, layers {}, datatypes {})".format(len(self.polygons), sum([len(p) for p in self.polygons]), list(set(self.layers)), list(set(self.datatypes)))


class _SubPath(object):
    """
    Single path component.
    """
    __slots__ = 'x', 'dx', 'off', 'wid', 'h', 'err', 'max_evals'

    def __init__(self, x, dx, off, wid, tolerance, max_evals):
        self.x = x
        self.dx = dx
        self.off = off
        self.wid = wid
        self.err = tolerance ** 2
        self.h = 0.5 / max_evals
        self.max_evals = max_evals

    def __str__(self):
        return 'SubPath ({} - {})'.format(self(0, 1e-6, 0), self(1, 1e-6, 0))

    def __call__(self, u, arm):
        v = self.dx(u, self.h)[::-1] * _pmone
        v /= (v[0]**2 + v[1]**2)**0.5
        x = self.x(u) + self.off(u) * v
        if arm == 0:
            return x
        u0 = max(0, u - self.h)
        u1 = min(1, u + self.h)
        w = (self(u1, 0) - self(u0, 0))[::-1] * _pmone
        w /= (w[0]**2 + w[1]**2)**0.5
        if arm < 0:
            return x - 0.5 * self.wid(u) * w
        return x + 0.5 * self.wid(u) * w

    def grad(self, u, arm):
        u0 = max(0, u - self.h)
        u1 = min(1, u + self.h)
        return (self(u1, arm) - self(u0, arm)) / (u1 - u0)

    def points(self, u0, u1, arm):
        u = [u0, u1]
        pts = [numpy.array(self(u[0], arm)), numpy.array(self(u[1], arm))]
        i = 1
        while i < len(pts) < self.max_evals:
            f = 0.2
            while f < 1:
                test_u = u[i - 1] * (1 - f) +  u[i] * f
                test_pt = numpy.array(self(test_u, arm))
                test_err = pts[i - 1] * (1 - f) +  pts[i] * f - test_pt
                if test_err[0]**2 + test_err[1]**2 > self.err:
                    u.insert(i, test_u)
                    pts.insert(i, test_pt)
                    f = 1
                    i -=1
                else:
                    f += 0.3
            i += 1
        return pts


class FlexPath(object):
    """
    Path object.

    This class keeps information about the constructive parameters of
    the path and calculates its boundaries only upon request.

    It can be stored as a proper path element in the GDSII format,
    unlike `Path`.  In this case, the width must be constant along the
    whole path.

    Parameters
    ----------
    points : array-like[N][2]
        Points along the center of the path.
    width : number, list
        Width of each parallel path being created.  The number of
        parallel paths being created is defined by the length of this
        list.
    offset : number, list
        Offsets of each parallel path from the center.  If `width` is
        not a list, the length of this list is used to determine the
        number of parallel paths being created.  Otherwise, offset must
        be a list with the same length as width, or a number, which is
        used as distance between adjacent paths.
    corners : 'natural', 'miter', 'bevel', 'round', 'smooth', 'circular bend', callable, list
        Type of joins.  A callable must receive 6 arguments (vertex and
        direction vector from both segments being joined, the center
        and width of the path) and return a list of vertices that make
        the join.  A list can be used to define the join for each
        parallel path.
    ends : 'flush', 'extended', 'round', 'smooth', 2-tuple, callable, list
        Type of end caps for the paths.  A 2-element tuple represents
        the start and end extensions to the paths.  A callable must
        receive 4 arguments (vertex and direction vectors from both
        sides of the path and return a list of vertices that make the
        end cap.  A list can be used to define the end type for each
        parallel path.
    bend_radius : number, list
        Bend radii for each path when `corners` is 'circular bend'.
        It has no effect for other corner types.
    tolerance : number
        Tolerance used to draw the paths and calculate joins.
    precision : number
        Precision for rounding the coordinates of vertices when
        fracturing the final polygonal boundary.
    max_points : integer
        If the number of points in the polygonal path boundary is
        greater than `max_points`, it will be fractured in smaller
        polygons with at most `max_points` each.  If `max_points` is
        zero no fracture will occur.
    gdsii_path : bool
        If True, treat this object as a GDSII path element.
        Otherwise, it will be converted into polygonal boundaries when
        required.
    width_transform : bool
        If `gdsii_path` is True, this flag indicates whether the width
        of the path should transform when scaling this object.  It has
        no effect when `gdsii_path` is False.
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
    The value of `tolerance` should not be smaller than `precision`,
    otherwise there would be wasted computational effort in calculating
    the paths.
    """
    __slots__ = ('n', 'ends', 'corners', 'points', 'offsets', 'widths', 'layers', 'datatypes',
                 'tolerance', 'precision', 'max_points', 'gdsii_path', 'width_transform',
                 'bend_radius', '_polygon_dict')

    _pathtype_dict = {'flush': 0, 'round': 1, 'extended': 2, 'smooth': 1}

    def __init__(self, points, width, offset=0, corners='natural', ends='flush', bend_radius=None,
                 tolerance=0.01, precision=1e-3, max_points=199, gdsii_path=False, width_transform=True,
                 layer=0, datatype=0):
        self._polygon_dict = None
        if isinstance(width, list):
            self.n = len(width)
            self.widths = width
            if isinstance(offset, list):
                self.offsets = offset
            else:
                self.offsets = [(i - 0.5 * (self.n - 1)) * offset for i in range(self.n)]
        else:
            if isinstance(offset, list):
                self.n = len(offset)
                self.offsets = offset
            else:
                self.n = 1
                self.offsets = [offset]
            self.widths = [width] * self.n
        self.widths = numpy.tile(self.widths, (len(points), 1))
        self.offsets = numpy.tile(self.offsets, (len(points), 1))
        self.points = numpy.array(points)
        if isinstance(ends, list):
            self.ends = [ends[i % len(ends)] for i in range(self.n)]
        else:
            self.ends = [ends for _ in range(self.n)]
        if isinstance(corners, list):
            self.corners = [corners[i % len(corners)] for i in range(self.n)]
        else:
            self.corners = [corners for _ in range(self.n)]
        if isinstance(bend_radius, list):
            self.bend_radius = [bend_radius[i % len(bend_radius)] for i in range(self.n)]
        else:
            self.bend_radius = [bend_radius for _ in range(self.n)]
        if isinstance(layer, list):
            self.layers = [layer[i % len(layer)] for i in range(self.n)]
        else:
            self.layers = [layer] * self.n
        if isinstance(datatype, list):
            self.datatypes = [datatype[i % len(datatype)] for i in range(self.n)]
        else:
            self.datatypes = [datatype] * self.n
        self.tolerance = tolerance
        self.precision = precision
        self.max_points = max_points
        self.gdsii_path = gdsii_path
        self.width_transform = width_transform
        if self.gdsii_path:
            if any(end == 'smooth' or callable(end) for end in self.ends):
                warnings.warn("[GDSPY] Smooth and custom end caps are not supported in `FlexPath` with `gdsii_path == True`.", stacklevel=3)
            if any(corner != 'natural' and corner != 'circular bend' for corner in self.corners):
                warnings.warn("[GDSPY] Corner specification not supported in `FlexPath` with `gdsii_path == True`.", stacklevel=3)

    def __str__(self):
        if self.n > 1:
            return "FlexPath (x{}, {} segments, layers {}, datatypes {})".format(self.n, self.points.shape[0], self.layers, self.datatypes)
        else:
            return "FlexPath ({} segments, layer {}, datatype {})".format(self.points.shape[0], self.layers[0], self.datatypes[0])

    def get_polygons(self, by_spec=False):
        """
        Calculate the polygonal boundaries described by this path.

        Parameters
        ----------
        by_spec : bool
            If True, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype).

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons (if
            `by_spec` is True).
        """
        if self._polygon_dict is None:
            self._polygon_dict = {}
            if self.points.shape[0] == 2:
                # Common case: single path with 2 points
                un = self.points[1] - self.points[0]
                un = un[::-1] * _mpone / (un[0]**2 + un[1]**2)**0.5
                for kk in range(self.n):
                    end = self.ends[kk]
                    pts = numpy.array((self.points[0] + un * self.offsets[0, kk],
                                       self.points[1] + un * self.offsets[1, kk]))
                    vn = pts[1] - pts[0]
                    vn = vn[::-1] * _mpone / (vn[0]**2 + vn[1]**2)**0.5
                    v = (vn * (0.5 * self.widths[0, kk]), vn * (0.5 * self.widths[1, kk]))
                    poly = numpy.array((pts[0] - v[0], pts[0] + v[0],
                                        pts[1] + v[1], pts[1] - v[1]))
                    if end != 'flush':
                        v0 = poly[3] - poly[0]
                        v1 = poly[2] - poly[1]
                        if callable(end):
                            cap0 = end(poly[0], -v0, poly[1], v1)
                            cap1 = end(poly[3], v0, poly[2], -v1)
                            poly = numpy.array(cap0[::-1] + cap1)
                        elif end == 'smooth':
                            angles = [numpy.arctan2(-v0[1], -v0[0]), numpy.arctan2(v1[1], v1[0])]
                            cta, ctb = _hobby(poly[:2], angles)
                            f = _func_bezier(numpy.array([poly[0], cta[0], ctb[0], poly[1]]))
                            tol = self.tolerance ** 2
                            uu = [0, 1]
                            fu = [f(0), f(1)]
                            iu = 1
                            while iu < len(fu):
                                test_u = 0.5 * (uu[iu - 1] +  uu[iu])
                                test_pt = f(test_u)
                                test_err = 0.5 * (fu[iu - 1] +  fu[iu]) - test_pt
                                if test_err[0]**2 + test_err[1]**2 > tol:
                                    uu.insert(iu, test_u)
                                    fu.insert(iu, test_pt)
                                else:
                                    iu += 1
                            cap = fu
                            cta, ctb = _hobby(poly[2:], [angles[1], angles[0]])
                            f = _func_bezier(numpy.array([poly[2], cta[0], ctb[0], poly[3]]))
                            tol = self.tolerance ** 2
                            uu = [0, 1]
                            fu = [f(0), f(1)]
                            iu = 1
                            while iu < len(fu):
                                test_u = 0.5 * (uu[iu - 1] +  uu[iu])
                                test_pt = f(test_u)
                                test_err = 0.5 * (fu[iu - 1] +  fu[iu]) - test_pt
                                if test_err[0]**2 + test_err[1]**2 > tol:
                                    uu.insert(iu, test_u)
                                    fu.insert(iu, test_pt)
                                else:
                                    iu += 1
                            cap.extend(fu)
                            poly = numpy.array(cap)
                        elif end == 'round':
                            v = pts[1] - pts[0]
                            r = 0.5 * self.widths[0, kk]
                            np = max(5, 1 + int(_halfpi / numpy.arccos(1 - self.tolerance / r) + 0.5))
                            ang = numpy.linspace(_halfpi, -_halfpi, np) + numpy.arctan2(-v[1], -v[0])
                            poly = pts[0] + r * numpy.vstack((numpy.cos(ang), numpy.sin(ang))).T
                            r = 0.5 * self.widths[1, kk]
                            np = max(5, 1 + int(_halfpi / numpy.arccos(1 - self.tolerance / r) + 0.5))
                            ang = numpy.linspace(_halfpi, -_halfpi, np) + numpy.arctan2(v[1], v[0])
                            poly = numpy.vstack((poly, pts[1] + r * numpy.vstack((numpy.cos(ang), numpy.sin(ang))).T))
                        else: # 'extended'/list
                            v = pts[1] - pts[0]
                            v /= (v[0]**2 + v[1]**2)**0.5
                            if end == 'extended':
                                v0 = 0.5 * self.widths[0, kk] * v
                                v1 = 0.5 * self.widths[1, kk] * v
                            else:
                                v0 = end[0] * v
                                v1 = end[1] * v
                            if self.widths[0, kk] == self.widths[1, kk]:
                                poly[0] -= v0
                                poly[1] -= v0
                                poly[2] += v1
                                poly[3] += v1
                            else:
                                poly = numpy.array((poly[0], poly[0] - v0, poly[1] - v0, poly[1],
                                                    poly[2], poly[2] + v1, poly[3] + v1, poly[3]))
                    polygons = [poly]
                    if self.max_points > 4 and poly.shape[0] > self.max_points:
                        ii = 0
                        while ii < len(polygons):
                            if len(polygons[ii]) > self.max_points:
                                pts0 = sorted(polygons[ii][:, 0])
                                pts1 = sorted(polygons[ii][:, 1])
                                ncuts = len(pts0) // self.max_points
                                if pts0[-1] - pts0[0] > pts1[-1] - pts1[0]:
                                    # Vertical cuts
                                    cuts = [pts0[int(i * len(pts0) / (ncuts + 1.0) + 0.5)] for i in range(1, ncuts + 1)]
                                    chopped = clipper._chop(polygons[ii], cuts, 0, 1 / self.precision)
                                else:
                                    # Horizontal cuts
                                    cuts = [pts1[int(i * len(pts1) / (ncuts + 1.0) + 0.5)] for i in range(1, ncuts + 1)]
                                    chopped = clipper._chop(polygons[ii], cuts, 1, 1 / self.precision)
                                polygons.pop(ii)
                                polygons.extend(numpy.array(x) for x in itertools.chain.from_iterable(chopped))
                            else:
                                ii += 1
                    key = (self.layers[kk], self.datatypes[kk])
                    if key in self._polygon_dict:
                        self._polygon_dict[key].extend(polygons)
                    else:
                        self._polygon_dict[key] = polygons
            else:
                # More than 1 path or more than 2 points
                un = self.points[1:] - self.points[:-1]
                un = un[:, ::-1] * _mpone / ((un[:, 0]**2 + un[:, 1]**2)**0.5).reshape((un.shape[0], 1))
                for kk in range(self.n):
                    corner = self.corners[kk]
                    end = self.ends[kk]
                    if any(self.offsets[:, kk] != 0):
                        pts = numpy.empty(self.points.shape)
                        sa = self.points[:-1] + un * self.offsets[:-1, kk:kk + 1]
                        sb = self.points[1:] + un * self.offsets[1:, kk:kk + 1]
                        vn = sb - sa
                        den = vn[1:, 0] * vn[:-1, 1] - vn[1:, 1] * vn[:-1, 0]
                        idx = numpy.nonzero(den**2 < 1e-12 * (vn[1:, 0]**2 + vn[1:, 1]**2) * (vn[:-1, 0]**2 + vn[:-1, 1]**2))[0]
                        if len(idx) > 0:
                            den[idx] = 1
                        ds = sb[:-1] - sa[1:]
                        u0 = (vn[1:, 1] * ds[:, 0] - vn[1:, 0] * ds[:, 1]) / den
                        u1 = (vn[:-1, 1] * ds[:, 0] - vn[:-1, 0] * ds[:, 1]) / den
                        if any(u0 < -1) or any(u1 > 1):
                            warnings.warn("[GDSPY] Possible inconsistency found in `FlexPath` due to sharp corner.")
                        pts[1:-1] = sb[:-1] + u0.reshape((u0.shape[0], 1)) * vn[:-1]
                        if len(idx) > 0:
                            pts[idx + 1] = 0.5 * (sa[idx + 1] + sb[idx])
                        pts[0] = sa[0]
                        pts[-1] = sb[-1]
                    else:
                        pts = self.points
                    vn = pts[1:] - pts[:-1]
                    vn = vn[:, ::-1] * _mpone / ((vn[:, 0]**2 + vn[:, 1]**2)**0.5).reshape((vn.shape[0], 1))
                    arms = [[], []]
                    caps = [[], []]
                    for ii in (0, 1):
                        sign = -1 if ii == 0 else 1
                        pa = pts[:-1] + vn * (sign * 0.5 * self.widths[:-1, kk:kk + 1])
                        pb = pts[1:] + vn * (sign * 0.5 * self.widths[1:, kk:kk + 1])
                        vec = pb - pa
                        caps[0].append(pa[0])
                        caps[1].append(pb[-1])
                        for jj in range(1, self.points.shape[0] - 1):
                            p0 = pb[jj - 1]
                            v0 = vec[jj - 1]
                            p1 = pa[jj]
                            v1 = vec[jj]
                            half_w = 0.5 * self.widths[jj, kk]
                            if corner == 'natural':
                                v0 = v0 * (half_w / (v0[0]**2 + v0[1]**2)**0.5)
                                v1 = v1 * (half_w / (v1[0]**2 + v1[1]**2)**0.5)
                                den = v1[1] * v0[0] - v1[0] * v0[1]
                                if den**2 < 1e-12 * half_w**4:
                                    u0 = u1 = 0
                                    p = 0.5 * (p0 + p1)
                                else:
                                    dx = p1[0] - p0[0]
                                    dy = p1[1] - p0[1]
                                    u0 = (v1[1] * dx - v1[0] * dy) / den
                                    u1 = (v0[1] * dx - v0[0] * dy) / den
                                    p = 0.5 * (p0 + v0 * u0 + p1 + v1 * u1)
                                if u0 < 0 and u1 > 0:
                                    arms[ii].append(p)
                                elif u0 <= 1 and u1 >= -1:
                                    arms[ii].append(0.5 * (p0 + min(1, u0) * v0 + p1 + max(-1, u1) * v1))
                                else:
                                    arms[ii].append(p0 + min(1, u0) * v0)
                                    arms[ii].append(p1 + max(-1, u1) * v1)
                            elif corner == 'circular bend':
                                v2 = p0 - pts[jj]
                                direction = v0[0] * v1[1] - v0[1] * v1[0]
                                if direction == 0:
                                    arms[ii].append(0.5 * (p0 + p1))
                                else:
                                    if direction > 0:
                                        a0 = numpy.arctan2(-v0[0], v0[1])
                                        a1 = numpy.arctan2(-v1[0], v1[1])
                                    else:
                                        a0 = numpy.arctan2(v0[0], -v0[1])
                                        a1 = numpy.arctan2(v1[0], -v1[1])
                                    if abs(a1 - a0) > numpy.pi:
                                        if a1 > a0:
                                            a0 += 2 * numpy.pi
                                        else:
                                            a1 += 2 * numpy.pi
                                    side = direction * (v0[0] * v2[1] - v0[1] * v2[0])
                                    if side > 0:
                                        r = self.bend_radius[kk] - half_w
                                    else:
                                        r = self.bend_radius[kk] + half_w
                                    da = 0.5 * abs(a1 - a0)
                                    d = self.bend_radius[kk] * numpy.tan(da) / (v0[0]**2 + v0[1]**2)**0.5
                                    np = max(2, 1 + int(da / numpy.arccos(1 - self.tolerance / r) + 0.5))
                                    angles = numpy.linspace(a0, a1, np)
                                    points = r * numpy.vstack((numpy.cos(angles), numpy.sin(angles))).T
                                    arms[ii].extend(points - points[0] + p0 - d * v0)
                            elif callable(corner):
                                arms[ii].extend(corner(p0, v0, p1, v1, pts[jj], self.widths[jj, kk]))
                            else:
                                den = v1[1] * v0[0] - v1[0] * v0[1]
                                lim = 1e-12 * (v0[0]**2 + v0[1]**2) * (v1[0]**2 + v1[1]**2)
                                if den**2 < lim:
                                    u0 = u1 = 0
                                    p = 0.5 * (p0 + p1)
                                else:
                                    dx = p1[0] - p0[0]
                                    dy = p1[1] - p0[1]
                                    u0 = (v1[1] * dx - v1[0] * dy) / den
                                    u1 = (v0[1] * dx - v0[0] * dy) / den
                                    p = 0.5 * (p0 + v0 * u0 + p1 + v1 * u1)
                                if corner == 'miter':
                                    arms[ii].append(p)
                                elif u0 <= 0 and u1 >= 0:
                                    arms[ii].append(p)
                                elif corner == 'bevel':
                                    arms[ii].append(p0)
                                    arms[ii].append(p1)
                                elif corner == 'round':
                                    if v0[1] * v1[0] - v0[0] * v1[1] < 0:
                                        a0 = numpy.arctan2(-v0[0], v0[1])
                                        a1 = numpy.arctan2(-v1[0], v1[1])
                                    else:
                                        a0 = numpy.arctan2(v0[0], -v0[1])
                                        a1 = numpy.arctan2(v1[0], -v1[1])
                                    if abs(a1 - a0) > numpy.pi:
                                        if a0 < a1:
                                            a0 += 2 * numpy.pi
                                        else:
                                            a1 += 2 * numpy.pi
                                    np = max(4, 1 + int(0.5 * abs(a1 - a0) / numpy.arccos(1 - self.tolerance / half_w) + 0.5))
                                    angles = numpy.linspace(a0, a1, np)
                                    arms[ii].extend(pts[jj] + half_w * numpy.vstack((numpy.cos(angles),
                                                                                numpy.sin(angles))).T)
                                elif corner == 'smooth':
                                    angles = [numpy.arctan2(v0[1], v0[0]), numpy.arctan2(v1[1], v1[0])]
                                    bezpts = numpy.vstack((p0, p1))
                                    cta, ctb = _hobby(bezpts, angles)
                                    f = _func_bezier(numpy.array([bezpts[0], cta[0], ctb[0], bezpts[1]]))
                                    tol = self.tolerance ** 2
                                    uu = [0, 1]
                                    fu = [f(0), f(1)]
                                    iu = 1
                                    while iu < len(fu):
                                        test_u = 0.5 * (uu[iu - 1] +  uu[iu])
                                        test_pt = f(test_u)
                                        test_err = 0.5 * (fu[iu - 1] +  fu[iu]) - test_pt
                                        if test_err[0]**2 + test_err[1]**2 > tol:
                                            uu.insert(iu, test_u)
                                            fu.insert(iu, test_pt)
                                        else:
                                            iu += 1
                                    arms[ii].extend(fu)
                    if end != 'flush':
                        for ii in (0, 1):
                            if callable(end):
                                vecs = [caps[ii][0] - arms[0][-ii], arms[1][-ii] - caps[ii][1]]
                                caps[ii] = end(caps[ii][0], vecs[0], caps[ii][1], vecs[1])
                            elif end == 'smooth':
                                points = numpy.array(caps[ii])
                                vecs = [caps[ii][0] - arms[0][-ii], arms[1][-ii] - caps[ii][1]]
                                angles = [numpy.arctan2(vecs[0][1], vecs[0][0]),
                                          numpy.arctan2(vecs[1][1], vecs[1][0])]
                                cta, ctb = _hobby(points, angles)
                                f = _func_bezier(numpy.array([points[0], cta[0], ctb[0], points[1]]))
                                tol = self.tolerance ** 2
                                uu = [0, 1]
                                fu = [f(0), f(1)]
                                iu = 1
                                while iu < len(fu):
                                    test_u = 0.5 * (uu[iu - 1] +  uu[iu])
                                    test_pt = f(test_u)
                                    test_err = 0.5 * (fu[iu - 1] +  fu[iu]) - test_pt
                                    if test_err[0]**2 + test_err[1]**2 > tol:
                                        uu.insert(iu, test_u)
                                        fu.insert(iu, test_pt)
                                    else:
                                        iu += 1
                                caps[ii] = fu
                            elif end == 'round':
                                v = pts[0] - pts[1] if ii == 0 else pts[-1] - pts[-2]
                                r = 0.5 * self.widths[-ii, kk]
                                np = max(5, 1 + int(_halfpi / numpy.arccos(1 - self.tolerance / r) + 0.5))
                                ang = (2 * ii - 1) * numpy.linspace(-_halfpi, _halfpi, np) + numpy.arctan2(v[1], v[0])
                                caps[ii] = list(pts[-ii] + r * numpy.vstack((numpy.cos(ang), numpy.sin(ang))).T)
                            else: # 'extended'/list
                                v = pts[0] - pts[1] if ii == 0 else pts[-1] - pts[-2]
                                v = v / (v[0]**2 + v[1]**2)**0.5
                                w = (2 * ii - 1) * v[::-1] * _pmone
                                r = 0.5 * self.widths[-ii, kk]
                                d = r if end == 'extended' else end[ii]
                                caps[ii] = [pts[-ii] + r * w, pts[-ii] + r * w + d * v,
                                            pts[-ii] - r * w + d * v, pts[-ii] - r * w]
                    poly = caps[0][::-1]
                    poly.extend(arms[0])
                    poly.extend(caps[1])
                    poly.extend(arms[1][::-1])
                    polygons = [numpy.array(poly)]
                    if self.max_points > 4 and polygons[0].shape[0] > self.max_points:
                        ii = 0
                        while ii < len(polygons):
                            if len(polygons[ii]) > self.max_points:
                                pts0 = sorted(polygons[ii][:, 0])
                                pts1 = sorted(polygons[ii][:, 1])
                                ncuts = len(pts0) // self.max_points
                                if pts0[-1] - pts0[0] > pts1[-1] - pts1[0]:
                                    # Vertical cuts
                                    cuts = [pts0[int(i * len(pts0) / (ncuts + 1.0) + 0.5)] for i in range(1, ncuts + 1)]
                                    chopped = clipper._chop(polygons[ii], cuts, 0, 1 / self.precision)
                                else:
                                    # Horizontal cuts
                                    cuts = [pts1[int(i * len(pts1) / (ncuts + 1.0) + 0.5)] for i in range(1, ncuts + 1)]
                                    chopped = clipper._chop(polygons[ii], cuts, 1, 1 / self.precision)
                                polygons.pop(ii)
                                polygons.extend(numpy.array(x) for x in itertools.chain.from_iterable(chopped))
                            else:
                                ii += 1
                    key = (self.layers[kk], self.datatypes[kk])
                    if key in self._polygon_dict:
                        self._polygon_dict[key].extend(polygons)
                    else:
                        self._polygon_dict[key] = polygons
        if by_spec:
            return libcopy.deepcopy(self._polygon_dict)
        else:
            return list(itertools.chain.from_iterable(self._polygon_dict.values()))

    def to_polygonset(self):
        """
        Create a `PolygonSet` representation of this object.

        The resulting object will be fractured according to the
        parameter `max_points` used when instantiating this object.

        Returns
        -------
        out : `PolygonSet` or None
            A `PolygonSet` that contains all boundaries for this path.
            If the path is empty, returns None.
        """
        if self.points.shape[0] < 2:
            return None
        polygons = self.get_polygons(True)
        pol = PolygonSet([])
        for k, v in polygons.items():
            pol.layers.extend([k[0]] * len(v))
            pol.datatypes.extend([k[1]] * len(v))
            pol.polygons.extend(v)
        return pol.fracture(self.max_points, self.precision)

    def to_gds(self, multiplier):
        """
        Convert this object to a series of GDSII elements.

        If `FlexPath.gdsii_path` is True, GDSII path elements are
        created instead of boundaries.  Such paths do not support
        variable widths, but their memeory footprint is smaller than
        full polygonal boundaries.

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
        if len(self.points) == 0:
            return b''
        if self.gdsii_path:
            sign = 1 if self.width_transform else -1
        else:
            return self.to_polygonset().to_gds(multiplier)
        data = []
        un = self.points[1:] - self.points[:-1]
        un = un[:, ::-1] * _mpone / ((un[:, 0]**2 + un[:, 1]**2)**0.5).reshape((un.shape[0], 1))
        for ii in range(self.n):
            pathtype = 0 if callable(self.ends[ii]) else FlexPath._pathtype_dict.get(self.ends[ii], 4)
            data.append(struct.pack('>4Hh2Hh2Hh2Hl', 4, 0x0900, 6, 0x0D02, self.layers[ii],
                                    6, 0x0E02, self.datatypes[ii], 6, 0x2102, pathtype,
                                    8, 0x0F03, sign * int(round(self.widths[0, ii] * multiplier))))
            if pathtype == 4:
                data.append(struct.pack('>2Hl2Hl', 8, 0x3003, int(round(self.ends[ii][0] * multiplier)),
                                        8, 0x3103, int(round(self.ends[ii][1] * multiplier))))
            if any(self.offsets[:, ii] != 0):
                points = numpy.zeros(self.points.shape)
                sa = self.points[:-1] + un * self.offsets[:-1, ii:ii + 1]
                sb = self.points[1:] + un * self.offsets[1:, ii:ii + 1]
                vn = sb - sa
                den = vn[1:, 0] * vn[:-1, 1] - vn[1:, 1] * vn[:-1, 0]
                idx = numpy.nonzero(den**2 < 1e-12 * (vn[1:, 0]**2 + vn[1:, 1]**2) * (vn[:-1, 0]**2 + vn[:-1, 1]**2))[0]
                if len(idx) > 0:
                    den[idx] = 1
                u0 = (vn[1:, 1] * (sb[:-1, 0] - sa[1:, 0])
                      - vn[1:, 0] * (sb[:-1, 1] - sa[1:, 1])) / den
                points[1:-1] = sb[:-1] + u0.reshape((u0.shape[0], 1)) * vn[:-1]
                if len(idx) > 0:
                    points[idx + 1] = 0.5 * (sa[idx + 1] + sb[idx])
                points[0] = sa[0]
                points[-1] = sb[-1]
            else:
                points = self.points
            if self.corners[ii] == 'circular bend':
                r = self.bend_radius[ii]
                p0 = points[0]
                p1 = points[1]
                v0 = p1 - p0
                bends = [p0]
                for jj in range(1, points.shape[0] - 1):
                    p2 = points[jj + 1]
                    v1 = p2 - p1
                    direction = v0[0] * v1[1] - v0[1] * v1[0]
                    if direction == 0:
                        bends.append(p1)
                    else:
                        if direction > 0:
                            a0 = numpy.arctan2(-v0[0], v0[1])
                            a1 = numpy.arctan2(-v1[0], v1[1])
                        elif direction < 0:
                            a0 = numpy.arctan2(v0[0], -v0[1])
                            a1 = numpy.arctan2(v1[0], -v1[1])
                        if abs(a1 - a0) > numpy.pi:
                            if a1 > a0:
                                a0 += 2 * numpy.pi
                            else:
                                a1 += 2 * numpy.pi
                        da = 0.5 * abs(a1 - a0)
                        d = r * numpy.tan(da) / (v0[0]**2 + v0[1]**2)**0.5
                        np = max(2, 1 + int(da / numpy.arccos(1 - self.tolerance / (r + 0.5 * self.widths[0, ii])) + 0.5))
                        angles = numpy.linspace(a0, a1, np)
                        bpts = r * numpy.vstack((numpy.cos(angles), numpy.sin(angles))).T
                        bends.extend(bpts - bpts[0] + p1 - d * v0)
                    p0 = p1
                    p1 = p2
                    v0 = v1
                bends.append(p1)
                points = numpy.array(bends)
            points = numpy.round(points * multiplier).astype('>i4')
            if points.shape[0] > 8191:
                warnings.warn("[GDSPY] Paths with more than 8191 are not supported by the official GDSII specification.  This GDSII file might not be compatible with all readers.", stacklevel=4)
                i0 = 0
                while i0 < points.shape[0]:
                    i1 = min(i0 + 8191, points.shape[0])
                    data.append(struct.pack('>2H', 4 + 8 * (i1 - i0), 0x1003))
                    data.append(points[i0:i1].tostring())
                    i0 = i1
            else:
                data.append(struct.pack('>2H', 4 + 8 * points.shape[0], 0x1003))
                data.append(points.tostring())
            data.append(struct.pack('>2H', 4, 0x1100))
        return b''.join(data)

    def area(self, by_spec=False):
        """
        Calculate the total area of this object.

        This functions creates a `PolgonSet` from this object and
        calculates its area, which means it is computationally
        expensive.

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
        return self.to_polygonset().area(by_spec)

    def translate(self, dx, dy):
        """
        Translate this path.

        Parameters
        ----------
        dx : number
            Distance to move in the x-direction
        dy : number
            Distance to move in the y-direction

        Returns
        -------
        out : `FlexPath`
            This object.
        """
        self._polygon_dict = None
        self.points = self.points + numpy.array((dx, dy))
        return self

    def rotate(self, angle, center=(0, 0)):
        """
        Rotate this path.

        Parameters
        ----------
        angle : number
            The angle of rotation (in *radians*).
        center : array-like[2]
            Center point for the rotation.

        Returns
        -------
        out : `FlexPath`
            This object.
        """
        self._polygon_dict = None
        ca = numpy.cos(angle)
        sa = numpy.sin(angle) * _mpone
        c0 = numpy.array(center)
        x = self.points - c0
        self.points = x * ca + x[:, ::-1] * sa + c0
        return self

    def scale(self, scale, center=(0, 0)):
        """
        Scale this path.

        Parameters
        ----------
        scale : number
            Scaling factor.
        center : array-like[2]
            Center point for the scaling operation.

        Returns
        -------
        out : `FlexPath`
            This object.
        """
        self._polygon_dict = None
        c0 = numpy.array(center) * (1 - scale)
        self.points = self.points * scale + c0
        self.widths = self.widths * scale
        self.offsets = self.offsets * scale
        return self

    def transform(self, translation, rotation, scale, x_reflection, array_trans=None):
        """
        Apply a transform to this path.

        Parameters
        ----------
        translation : Numpy array[2]
            Translation vector.
        rotation : number
            Rotation angle.
        scale : number
            Scaling factor.
        x_reflection : bool
            Reflection around the first axis.
        array_trans : Numpy aray[2]
            Translation vector before rotation and reflection.

        Returns
        -------
        out : `FlexPath`
            This object.

        Notes
        -----
        Applies the transformations in the same order as a
        `CellReference` or a `CellArray`.  If `width_transform` is
        False, the widths are not scaled.
        """
        self._polygon_dict = None
        if translation is None:
            translation = _zero
        if array_trans is None:
            array_trans = _zero
        if rotation is None:
            cos = _one
            sin = _zero
        else:
            cos = numpy.cos(rotation) * _one
            sin = numpy.sin(rotation) * _mpone
        if scale is not None:
            cos = cos * scale
            sin = sin * scale
            array_trans = array_trans / scale
            if self.width_transform or not self.gdsii_path:
                self.widths = self.widths * scale
            self.offsets = self.offsets * scale
        if x_reflection:
            cos[1] = -cos[1]
            sin[0] = -sin[0]
        pts = self.points + array_trans
        self.points = pts * cos + pts[:,::-1] * sin + translation
        return self

    def segment(self, end_point, width=None, offset=None, relative=False):
        """
        Add a straight section to the path.

        Parameters
        ----------
        end_point : array-like[2]
            End position of the straight segment.
        width : number, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  A list can be used where each
            element defines the width for one of the parallel paths in
            this object.  This argument has no effect if the path was
            created with `gdsii_path` True.
        offset : number, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  A list
            can be used where each element defines the *absolute* offset
            (not offset increase) for one of the parallel paths in this
            object.
        relative : bool
            If True, `end_point` is used as an offset from the current
            path position, i.e., if the path is at (1, -2) and the
            `end_point` is (10, 25), the segment will be constructed
            from (1, -2) to (1 + 10, -2 + 25) = (11, 23).  Otherwise,
            `end_point` is used as an absolute coordinate.

        Returns
        -------
        out : `FlexPath`
            This object.
        """
        self._polygon_dict = None
        self.points = numpy.vstack((self.points, (self.points[-1] + numpy.array(end_point))
                                                 if relative else end_point))
        if self.gdsii_path or width is None:
            self.widths = numpy.vstack((self.widths, self.widths[-1]))
        elif hasattr(width, '__iter__'):
            self.widths = numpy.vstack((self.widths, width))
        else:
            self.widths = numpy.vstack((self.widths, numpy.repeat(width, self.n)))
        if offset is None:
            self.offsets = numpy.vstack((self.offsets, self.offsets[-1]))
        elif hasattr(offset, '__iter__'):
            self.offsets = numpy.vstack((self.offsets, offset))
        else:
            self.offsets = numpy.vstack((self.offsets, self.offsets[-1] + offset))
        return self

    def arc(self, radius, initial_angle, final_angle, width=None, offset=None):
        """
        Add a circular arc section to the path.

        Parameters
        ----------
        radius : number
            Radius of the circular arc.
        initial_angle : number
            Initial angle of the arc.
        final_angle : number
            Final angle of the arc.
        width : number, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  A list can be used where each
            element defines the width for one of the parallel paths in
            this object.  This argument has no effect if the path was
            created with `gdsii_path` True.
        offset : number, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  A list
            can be used where each element defines the *absolute* offset
            (not offset increase) for one of the parallel paths in this
            object.

        Returns
        -------
        out : `FlexPath`
            This object.
        """
        self._polygon_dict = None
        if self.gdsii_path:
            width = None
        if width is None:
            wid = self.widths[-1]
        elif hasattr(width, '__iter__'):
            wid = numpy.array(width)
        else:
            wid = numpy.full(self.n, width)
        if offset is None:
            off = self.offsets[-1]
        elif hasattr(offset, '__iter__'):
            off = numpy.array(offset)
        else:
            off = self.offsets[-1] + offset
        rmax = radius + max((self.offsets[-1] + self.widths[-1]).max(), (off + wid).max())
        np = max(3, 1 + int(0.5 * abs(final_angle - initial_angle) / numpy.arccos(1 - self.tolerance / rmax) + 0.5))
        ang = numpy.linspace(initial_angle, final_angle, np)
        pts = radius * numpy.vstack((numpy.cos(ang), numpy.sin(ang))).T
        self.points = numpy.vstack((self.points, pts[1:] + (self.points[-1] - pts[0])))
        if width is None:
            self.widths = numpy.vstack((self.widths, numpy.tile(wid, (np - 1, 1))))
        else:
            u = numpy.linspace(0, 1, np)[1:]
            self.widths = numpy.vstack((self.widths, numpy.outer(1 - u, self.widths[-1]) + numpy.outer(u, wid)))
        if offset is None:
            self.offsets = numpy.vstack((self.offsets, numpy.tile(off, (np - 1, 1))))
        else:
            u = numpy.linspace(0, 1, np)[1:]
            self.offsets = numpy.vstack((self.offsets, numpy.outer(1 - u, self.offsets[-1]) + numpy.outer(u, off)))
        return self

    def turn(self, radius, angle, width=None, offset=None):
        """
        Add a circular turn to the path.

        The initial angle of the arc is calculated from the last path
        segment.

        Parameters
        ----------
        radius : number
            Radius of the circular arc.
        angle : 'r', 'l', 'rr', 'll' or number
            Angle (in *radians*) of rotation of the path.  The values
            'r' and 'l' represent 90-degree turns cw and ccw,
            respectively; the values 'rr' and 'll' represent analogous
            180-degree turns.
        width : number, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  A list can be used where each
            element defines the width for one of the parallel paths in
            this object.  This argument has no effect if the path was
            created with `gdsii_path` True.
        offset : number, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  A list
            can be used where each element defines the *absolute* offset
            (not offset increase) for one of the parallel paths in this
            object.

        Returns
        -------
        out : `FlexPath`
            This object.
        """
        self._polygon_dict = None
        if self.points.shape[0] < 2:
            raise ValueError("[GDSPY] Cannot define initial angle for turn on a FlexPath withouth previous segments.")
        v = self.points[-1] - self.points[-2]
        angle = _angle_dic.get(angle, angle)
        initial_angle = numpy.arctan2(v[1], v[0]) + (_halfpi if angle < 0 else -_halfpi)
        self.arc(radius, initial_angle, initial_angle + angle, width, offset)
        return self

    def parametric(self, curve_function, width=None, offset=None, relative=True):
        """
        Add a parametric curve to the path.

        Parameters
        ----------
        curve_function : callable
            Function that defines the curve.  Must be a function of one
            argument (that varies from 0 to 1) that returns a 2-element
            Numpy array with the coordinates of the curve.
        width : number, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  A list can be used where each
            element defines the width for one of the parallel paths in
            this object.  This argument has no effect if the path was
            created with `gdsii_path` True.
        offset : number, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  A list
            can be used where each element defines the *absolute* offset
            (not offset increase) for one of the parallel paths in this
            object.
        relative : bool
            If True, the return values of `curve_function` are used as
            offsets from the current path position, i.e., to ensure a
            continuous path, ``curve_function(0)`` must be (0, 0).
            Otherwise, they are used as absolute coordinates.

        Returns
        -------
        out : `FlexPath`
            This object.
        """
        self._polygon_dict = None
        if self.gdsii_path:
            width = None
        if width is None:
            wid = self.widths[-1]
        elif hasattr(width, '__iter__'):
            wid = numpy.array(width)
        else:
            wid = numpy.full(self.n, width)
        if offset is None:
            off = self.offsets[-1]
        elif hasattr(offset, '__iter__'):
            off = numpy.array(offset)
        else:
            off = self.offsets[-1] + offset
        tol = self.tolerance**2
        u = [0, 1]
        pts = [numpy.array(curve_function(0)), numpy.array(curve_function(1))]
        i = 1
        while i < len(pts):
            f = 0.2
            while f < 1:
                test_u = u[i - 1] * (1 - f) +  u[i] * f
                test_pt = numpy.array(curve_function(test_u))
                test_err = pts[i - 1] * (1 - f) +  pts[i] * f - test_pt
                if test_err[0]**2 + test_err[1]**2 > tol:
                    u.insert(i, test_u)
                    pts.insert(i, test_pt)
                    f = 1
                    i -=1
                else:
                    f += 0.3
            i += 1
        pts = numpy.array(pts[1:])
        np = pts.shape[0] + 1
        self.points = numpy.vstack((self.points, (pts + self.points[-1]) if relative else pts))
        if width is None:
            self.widths = numpy.vstack((self.widths, numpy.tile(wid, (np - 1, 1))))
        else:
            u = numpy.linspace(0, 1, np)[1:]
            self.widths = numpy.vstack((self.widths, numpy.outer(1 - u, self.widths[-1]) + numpy.outer(u, wid)))
        if offset is None:
            self.offsets = numpy.vstack((self.offsets, numpy.tile(off, (np - 1, 1))))
        else:
            u = numpy.linspace(0, 1, np)[1:]
            self.offsets = numpy.vstack((self.offsets, numpy.outer(1 - u, self.offsets[-1]) + numpy.outer(u, off)))
        return self


    def bezier(self, points, width=None, offset=None, relative=True):
        """
        Add a Bezier curve to the path.

        A Bezier curve is added to the path starting from its current
        position and finishing at the last point in the `points` array.

        Parameters
        ----------
        points : array-like[N][2]
            Control points defining the Bezier curve.
        width : number, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  A list can be used where each
            element defines the width for one of the parallel paths in
            this object.  This argument has no effect if the path was
            created with `gdsii_path` True.
        offset : number, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  A list
            can be used where each element defines the *absolute* offset
            (not offset increase) for one of the parallel paths in this
            object.
        relative : bool
            If True, all coordinates in the `points` array are used as
            offsets from the current path position, i.e., if the path is
            at (1, -2) and the last point in the array is (10, 25), the
            constructed Bezier will end at (1 + 10, -2 + 25) = (11, 23).
            Otherwise, the points are used as absolute coordinates.

        Returns
        -------
        out : `FlexPath`
            This object.
        """
        self._polygon_dict = None
        if relative:
            ctrl = self.points[-1] + numpy.vstack(([(0, 0)], points))
        else:
            ctrl = numpy.vstack((self.points[-1:], points))
        self.parametric(_func_bezier(ctrl), width, offset, False)
        return self


    def smooth(self, points, angles=None, curl_start=1, curl_end=1, t_in=1, t_out=1, cycle=False,
               width=None, offset=None, relative=True):
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
        width : number, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  A list can be used where each
            element defines the width for one of the parallel paths in
            this object.  This argument has no effect if the path was
            created with `gdsii_path` True.
        offset : number, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  A list
            can be used where each element defines the *absolute* offset
            (not offset increase) for one of the parallel paths in this
            object.
        relative : bool
            If True, all coordinates in the `points` array are used as
            offsets from the current path position, i.e., if the path is
            at (1, -2) and the last point in the array is (10, 25), the
            constructed curve will end at (1 + 10, -2 + 25) = (11, 23).
            Otherwise, the points are used as absolute coordinates.

        Returns
        -------
        out : `FlexPath`
            This object.

        Notes
        -----
        Arguments `width` and `offset` are repeated for *each* cubic
        Bezier that composes this path element.

        References
        ----------
        .. [1] Hobby, J.D.  *Discrete Comput. Geom.* (1986) 1: 123.
           `DOI: 10.1007/BF02187690
           <https://doi.org/10.1007/BF02187690>`_
        """
        if relative:
            points = self.points[-1] + numpy.vstack(([(0, 0)], points))
        else:
            points = numpy.vstack((self.points[-1:], points))
        cta, ctb = _hobby(points, angles, curl_start, curl_end, t_in, t_out, cycle)
        for i in range(points.shape[0] - 1):
            self.bezier((cta[i], ctb[i], points[i + 1]), width, offset, False)
        if cycle:
            self.bezier((cta[-1], ctb[-1], points[0]), width, offset, False)
        return self


class RobustPath(object):
    """
    Path object with lazy evaluation.

    This class keeps information about the constructive parameters of
    the path and calculates its boundaries only upon request.  The
    benefits are that joins and path components can be calculated
    automatically to ensure continuity (except in extreme cases).

    It can be stored as a proper path element in the GDSII format,
    unlike `Path`.  In this case, the width must be constant along the
    whole path.

    The downside of `RobustPath` is that it is more computationally
    expensive than the other path classes.

    Parameters
    ----------
    initial_point : array-like[2]
        Starting position of the path.
    width : number, list
        Width of each parallel path being created.  The number of
        parallel paths being created is defined by the length of this
        list.
    offset : number, list
        Offsets of each parallel path from the center.  If `width` is
        not a list, the length of this list is used to determine the
        number of parallel paths being created.  Otherwise, offset must
        be a list with the same length as width, or a number, which is
        used as distance between adjacent paths.
    ends : 'flush', 'extended', 'round', 'smooth', 2-tuple, list
        Type of end caps for the paths.  A 2-element tuple represents
        the start and end extensions to the paths.  A list can be used
        to define the end type for each parallel path.
    tolerance : number
        Tolerance used to draw the paths and calculate joins.
    precision : number
        Precision for rounding the coordinates of vertices when
        fracturing the final polygonal boundary.
    max_points : integer
        If the number of points in the polygonal path boundary is
        greater than `max_points`, it will be fractured in smaller
        polygons with at most `max_points` each.  If `max_points` is
        zero no fracture will occur.
    max_evals : integer
        Limit to the maximal number of evaluations when calculating each
        path component.
    gdsii_path : bool
        If True, treat this object as a GDSII path element.
        Otherwise, it will be converted into polygonal boundaries when
        required.
    width_transform : bool
        If `gdsii_path` is True, this flag indicates whether the width
        of the path should transform when scaling this object.  It has
        no effect when `gdsii_path` is False.
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
    The value of `tolerance` should not be smaller than `precision`,
    otherwise there would be wasted computational effort in calculating
    the paths.
    """
    __slots__ = ('n', 'ends', 'x', 'offsets', 'widths', 'paths', 'layers', 'datatypes', 'tolerance',
                 'precision', 'max_points', 'max_evals', 'gdsii_path', 'width_transform', '_polygon_dict')

    _pathtype_dict = {'flush': 0, 'round': 1, 'extended': 2, 'smooth': 1}

    def __init__(self, initial_point, width, offset=0, ends='flush', tolerance=0.01, precision=1e-3,
                 max_points=199, max_evals=1000, gdsii_path=False, width_transform=True, layer=0, datatype=0):
        self._polygon_dict = None
        if isinstance(width, list):
            self.n = len(width)
            self.widths = width
            if isinstance(offset, list):
                self.offsets = offset
            else:
                self.offsets = [(i - 0.5 * (self.n - 1)) * offset for i in range(self.n)]
        else:
            if isinstance(offset, list):
                self.n = len(offset)
                self.offsets = offset
            else:
                self.n = 1
                self.offsets = [offset]
            self.widths = [width] * self.n
        self.x = numpy.array(initial_point)
        self.paths = [[] for _ in range(self.n)]
        if isinstance(ends, list):
            self.ends = [ends[i % len(ends)] for i in range(self.n)]
        else:
            self.ends = [ends for _ in range(self.n)]
        if isinstance(layer, list):
            self.layers = [layer[i % len(layer)] for i in range(self.n)]
        else:
            self.layers = [layer] * self.n
        if isinstance(datatype, list):
            self.datatypes = [datatype[i % len(datatype)] for i in range(self.n)]
        else:
            self.datatypes = [datatype] * self.n
        self.tolerance = tolerance
        self.precision = precision
        self.max_points = max_points
        self.max_evals = max_evals
        self.gdsii_path = gdsii_path
        self.width_transform = width_transform
        if self.gdsii_path and any(end == 'smooth' for end in self.ends):
            warnings.warn("[GDSPY] Smooth end caps not supported in `RobustPath` with `gdsii_path == True`.", stacklevel=3)

    def __str__(self):
        if self.n > 1:
            return "RobustPath (x{}, end at ({}, {}), length {}, layers {}, datatypes {})".format(self.n, self.x[0], self.x[1], len(self), self.layers, self.datatypes)
        else:
            return "RobustPath (end at ({}, {}), length {}, layer {}, datatype {})".format(self.x[0], self.x[1], len(self), self.layers[0], self.datatypes[0])

    def __len__(self):
        """
        Number of path components.
        """
        return len(self.paths[0])

    def __call__(self, u, arm=0):
        """
        Calculate the positions of each parallel path.

        Parameters
        ----------
        u : number
            Position along the `RobustPath` to compute.  This argument
            can range from 0 (start of the path) to ``len(self)`` (end
            of the path).
        arm : -1, 0, 1
            Wether to calculate one of the path boundaries (-1 or 1) or
            its central spine (0).

        Returns
        -------
        out : Numpy array[N, 2]
            Coordinates for each of the N parallel paths in this object.
        """
        i = int(u)
        u -= i
        if i == len(self.paths[0]):
            i -= 1
            u = 1
        return numpy.array([p[i](u, arm) for p in self.paths])

    def grad(self, u, arm=0, side='-'):
        """
        Calculate the direction vector of each parallel path.

        Parameters
        ----------
        u : number
            Position along the `RobustPath` to compute.  This argument
            can range from 0 (start of the path) to `len(self)` (end of
            the path).
        arm : -1, 0, 1
            Wether to calculate one of the path boundaries (-1 or 1) or
            its central spine (0).
        side : '-' or '+'
            At path joins, whether to calculate the direction using the
            component before or after the join.

        Returns
        -------
        out : Numpy array[N, 2]
            Direction vectors for each of the N parallel paths in this
            object.
        """
        i = int(u)
        u -= i
        if u == 0 and ((i > 0 and side == '-') or i == len(self.paths[0])):
            i -= 1
            u = 1
        return numpy.array([p[i].grad(u, arm) for p in self.paths])

    def width(self, u):
        """
        Calculate the width of each parallel path.

        Parameters
        ----------
        u : number
            Position along the `RobustPath` to compute.  This argument
            can range from 0 (start of the path) to `len(self)` (end of
            the path).

        Returns
        -------
        out : Numpy array[N]
            Width for each of the N parallel paths in this object.
        """
        i = int(u)
        u -= i
        if u == 0 and i == len(self.paths[0]):
            i -= 1
            u = 1
        return numpy.array([p[i].wid(u) for p in self.paths])

    def get_polygons(self, by_spec=False):
        """
        Calculate the polygonal boundaries described by this path.

        Parameters
        ----------
        by_spec : bool
            If True, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype).

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons (if
            `by_spec` is True).
        """
        if self._polygon_dict is None:
            self._polygon_dict = {}
            for path, end, layer, datatype in zip(self.paths, self.ends, self.layers, self.datatypes):
                poly = []
                for arm in [-1, 1]:
                    if not (end == 'flush'):
                        i = 0 if arm == 1 else -1
                        u = abs(i)
                        if end == 'smooth':
                            v1 = -arm * path[i].grad(u, -arm)
                            v2 = arm * path[i].grad(u, arm)
                            angles = [numpy.arctan2(v1[1], v1[0]), numpy.arctan2(v2[1], v2[0])]
                            points = numpy.array([path[i](u, -arm), path[i](u, arm)])
                            cta, ctb = _hobby(points, angles)
                            f = _func_bezier(numpy.array([points[0], cta[0], ctb[0], points[1]]))
                            tol = self.tolerance ** 2
                            uu = [0, 1]
                            fu = [f(0), f(1)]
                            iu = 1
                            while iu < len(fu) < self.max_evals:
                                test_u = 0.5 * (uu[iu - 1] +  uu[iu])
                                test_pt = f(test_u)
                                test_err = 0.5 * (fu[iu - 1] +  fu[iu]) - test_pt
                                if test_err[0]**2 + test_err[1]**2 > tol:
                                    uu.insert(iu, test_u)
                                    fu.insert(iu, test_pt)
                                else:
                                    iu += 1
                            poly.extend(fu[1:-1])
                        else:
                            p = path[i](u, 0)
                            v = -arm * path[i].grad(u, 0)
                            r = 0.5 * path[i].wid(u)
                            if end == 'round':
                                np = max(5, 1 + int(_halfpi / numpy.arccos(1 - self.tolerance / r) + 0.5))
                                ang = numpy.linspace(-_halfpi, _halfpi, np)[1:-1] + numpy.arctan2(v[1], v[0])
                                endpts = p + r * numpy.vstack((numpy.cos(ang), numpy.sin(ang))).T
                                poly.extend(endpts)
                            else:
                                v /= (v[0]**2 + v[1]**2)**0.5
                                w = v[::-1] * _pmone
                                d = r if end == 'extended' else end[u]
                                poly.append(p + d * v + r * w)
                                poly.append(p + d * v - r * w)
                    path_arm = []
                    start = 0
                    tol = self.tolerance**2
                    for sub0, sub1 in zip(path[:-1], path[1:]):
                        p0 = sub0(1, arm)
                        v0 = sub0.grad(1, arm)
                        p1 = sub1(0, arm)
                        v1 = sub1.grad(0, arm)
                        den = v1[1] * v0[0] - v1[0] * v0[1]
                        lim = 1e-12 * (v0[0]**2 + v0[1]**2) * (v1[0]**2 + v1[1]**2)
                        dx = p1[0] - p0[0]
                        dy = p1[1] - p0[1]
                        if den**2 < lim or dx**2 + dy**2 <= tol:
                            u0 = u1 = 0
                            px = 0.5 * (p0 + p1)
                        else:
                            u0 = (v1[1] * dx - v1[0] * dy) / den
                            u1 = (v0[1] * dx - v0[0] * dy) / den
                            px = 0.5 * (p0 + v0 * u0 + p1 + v1 * u1)
                        u0 = 1 + u0
                        if u0 < 1 and u1 > 0:
                            delta = sub0(u0, arm) - sub1(u1, arm)
                            err = delta[0]**2 + delta[1]**2
                            iters = 0
                            step = 0.5
                            while err > tol:
                                iters += 1
                                if iters > self.max_evals:
                                    warnings.warn('[GDSPY] Intersection not found.')
                                    break
                                du = delta * sub0.grad(u0, arm)
                                new_u0 = min(1, max(0, u0 - step * (du[0] + du[1])))
                                du = delta * sub1.grad(u1, arm)
                                new_u1 = min(1, max(0, u1 + step * (du[0] + du[1])))
                                new_delta = sub0(new_u0, arm) - sub1(new_u1, arm)
                                new_err = new_delta[0]**2 + new_delta[1]**2
                                if new_err >= err:
                                    step /= 2
                                    continue
                                u0 = new_u0
                                u1 = new_u1
                                delta = new_delta
                                err = new_err
                            px = 0.5 * (sub0(u0, arm) + sub1(u1, arm))
                        if u1 >= 0:
                            if u0 <= 1:
                                path_arm.extend(sub0.points(start, u0, arm)[:-1])
                            else:
                                path_arm.extend(sub0.points(start, 1, arm))
                                warnings.warn("[GDSPY] RobustPath join at ({}, {}) cannot be ensured.  Please check the resulting polygon.".format(path_arm[-1][0], path_arm[-1][1]), stacklevel=3)
                            start = u1
                        else:
                            if u0 <= 1:
                                path_arm.extend(sub0.points(start, u0, arm))
                                warnings.warn("[GDSPY] RobustPath join at ({}, {}) cannot be ensured.  Please check the resulting polygon.".format(path_arm[-1][0], path_arm[-1][1]), stacklevel=2)
                            else:
                                path_arm.extend(sub0.points(start, 1, arm))
                                path_arm.append(px)
                            start = 0
                    path_arm.extend(path[-1].points(start, 1, arm))
                    poly.extend(path_arm[::arm])
                polygons = [numpy.array(poly)]
                if self.max_points > 4 and polygons[0].shape[0] > self.max_points:
                    ii = 0
                    while ii < len(polygons):
                        if len(polygons[ii]) > self.max_points:
                            pts0 = sorted(polygons[ii][:, 0])
                            pts1 = sorted(polygons[ii][:, 1])
                            ncuts = len(pts0) // self.max_points
                            if pts0[-1] - pts0[0] > pts1[-1] - pts1[0]:
                                # Vertical cuts
                                cuts = [pts0[int(i * len(pts0) / (ncuts + 1.0) + 0.5)] for i in range(1, ncuts + 1)]
                                chopped = clipper._chop(polygons[ii], cuts, 0, 1 / self.precision)
                            else:
                                # Horizontal cuts
                                cuts = [pts1[int(i * len(pts1) / (ncuts + 1.0) + 0.5)] for i in range(1, ncuts + 1)]
                                chopped = clipper._chop(polygons[ii], cuts, 1, 1 / self.precision)
                            polygons.pop(ii)
                            polygons.extend(numpy.array(x) for x in itertools.chain.from_iterable(chopped))
                        else:
                            ii += 1
                key = (layer, datatype)
                if key in self._polygon_dict:
                    self._polygon_dict[key].extend(polygons)
                else:
                    self._polygon_dict[key] = polygons
        if by_spec:
            return libcopy.deepcopy(self._polygon_dict)
        else:
            return list(itertools.chain.from_iterable(self._polygon_dict.values()))

    def to_polygonset(self):
        """
        Create a `PolygonSet` representation of this object.

        The resulting object will be fractured according to the
        parameter `max_points` used when instantiating this object.

        Returns
        -------
        out : `PolygonSet` or None
            A `PolygonSet` that contains all boundaries for this path.
            If the path is empty, returns None.
        """
        if len(self.paths[0]) == 0:
            return None
        polygons = self.get_polygons(True)
        pol = PolygonSet([])
        for k, v in polygons.items():
            pol.layers.extend([k[0]] * len(v))
            pol.datatypes.extend([k[1]] * len(v))
            pol.polygons.extend(v)
        return pol.fracture(self.max_points, self.precision)

    def to_gds(self, multiplier):
        """
        Convert this object to a series of GDSII elements.

        If `RobustPath.gdsii_path` is True, GDSII path elements are
        created instead of boundaries.  Such paths do not support
        variable widths, but their memeory footprint is smaller than
        full polygonal boundaries.

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
        if len(self.paths[0]) == 0:
            return b''
        if self.gdsii_path:
            sign = 1 if self.width_transform else -1
        else:
            return self.to_polygonset().to_gds(multiplier)
        data = []
        for ii in range(self.n):
            pathtype = RobustPath._pathtype_dict.get(self.ends[ii], 4)
            data.append(struct.pack('>4Hh2Hh2Hh2Hl', 4, 0x0900, 6, 0x0D02, self.layers[ii],
                                    6, 0x0E02, self.datatypes[ii], 6, 0x2102, pathtype,
                                    8, 0x0F03, sign * int(round(self.widths[ii] * multiplier))))
            if pathtype == 4:
                data.append(struct.pack('>2Hl2Hl', 8, 0x3003, int(round(self.ends[ii][0] * multiplier)),
                                        8, 0x3103, int(round(self.ends[ii][1] * multiplier))))
            points = []
            for path in self.paths[ii]:
                new_points = numpy.round(numpy.array(path.points(0, 1, 0)) * multiplier)
                if len(points) > 0 and new_points[0, 0] == points[-1][-1, 0] and new_points[0, 1] == points[-1][-1, 1]:
                    points.append(new_points[1:])
                else:
                    points.append(new_points)
            points = numpy.vstack(points).astype('>i4')
            if points.shape[0] > 8191:
                warnings.warn("[GDSPY] Paths with more than 8191 are not supported by the official GDSII specification.  This GDSII file might not be compatible with all readers.", stacklevel=4)
                i0 = 0
                while i0 < points.shape[0]:
                    i1 = min(i0 + 8191, points.shape[0])
                    data.append(struct.pack('>2H', 4 + 8 * (i1 - i0), 0x1003))
                    data.append(points[i0:i1].tostring())
                    i0 = i1
            else:
                data.append(struct.pack('>2H', 4 + 8 * points.shape[0], 0x1003))
                data.append(points.tostring())
            data.append(struct.pack('>2H', 4, 0x1100))
        return b''.join(data)

    def area(self, by_spec=False):
        """
        Calculate the total area of this object.

        This functions creates a `PolgonSet` from this object and
        calculates its area, which means it is computationally
        expensive.

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
        return self.to_polygonset().area(by_spec)

    def translate(self, dx, dy):
        """
        Translate this path.

        Parameters
        ----------
        dx : number
            Distance to move in the x-direction
        dy : number
            Distance to move in the y-direction

        Returns
        -------
        out : `RobustPath`
            This object.
        """
        self._polygon_dict = None
        offset = numpy.array((dx, dy))
        self.x = self.x + offset
        for path in self.paths:
            for sub in path:
                sub.x = _func_multadd(sub.x, None, offset)
        return self

    def rotate(self, angle, center=(0, 0)):
        """
        Rotate this path.

        Parameters
        ----------
        angle : number
            The angle of rotation (in *radians*).
        center : array-like[2]
            Center point for the rotation.

        Returns
        -------
        out : `RobustPath`
            This object.
        """
        self._polygon_dict = None
        ca = numpy.cos(angle)
        sa = numpy.sin(angle) * _mpone
        c0 = numpy.array(center)
        x = self.x - c0
        self.x = x * ca + x[::-1] * sa + c0
        for path in self.paths:
            for sub in path:
                sub.x = _func_rotate(sub.x, ca, sa, c0)
                sub.dx = _func_rotate(sub.dx, ca, sa, nargs=2)
        return self

    def scale(self, scale, center=(0, 0)):
        """
        Scale this path.

        Parameters
        ----------
        scale : number
            Scaling factor.
        center : array-like[2]
            Center point for the scaling operation.

        Returns
        -------
        out : `RobustPath`
            This object.
        """
        self._polygon_dict = None
        c0 = numpy.array(center) * (1 - scale)
        self.x = self.x * scale + c0
        self.widths = [wid * scale for wid in self.widths]
        self.offsets = [off * scale for off in self.offsets]
        for path in self.paths:
            for sub in path:
                sub.x = _func_multadd(sub.x, scale, c0)
                sub.dx = _func_multadd(sub.dx, scale, None, nargs=2)
                sub.wid = _func_multadd(sub.wid, abs(scale), None)
                sub.off = _func_multadd(sub.off, scale, None)
        return self

    def transform(self, translation, rotation, scale, x_reflection, array_trans=None):
        """
        Apply a transform to this path.

        Parameters
        ----------
        translation : Numpy array[2]
            Translation vector.
        rotation : number
            Rotation angle.
        scale : number
            Scaling factor.
        x_reflection : bool
            Reflection around the first axis.
        array_trans : Numpy aray[2]
            Translation vector before rotation and reflection.

        Returns
        -------
        out : `RobustPath`
            This object.

        Notes
        -----
        Applies the transformations in the same order as a
        `CellReference` or a `CellArray`.  If `width_transform` is
        False, the widths are not scaled.
        """
        self._polygon_dict = None
        for ii in range(self.n):
            for sub in self.paths[ii]:
                sub.x = _func_trafo(sub.x, translation, rotation, scale, x_reflection, array_trans)
                sub.dx = _func_trafo(sub.dx, None, rotation, scale, x_reflection, None, nargs=2)
                if self.width_transform or not self.gdsii_path:
                    sub.wid = _func_multadd(sub.wid, scale, None)
                sub.off = _func_multadd(sub.off, scale, None)
            self.x[ii] = self.paths[-1].x(1)
            self.widths[ii] = self.paths[-1].wid(1)
            self.offsets[ii] = self.offsets[-1].off(1)
        return self

    def _parse_offset(self, arg, idx):
        if arg is None:
            return _func_const(self.offsets[idx])
        elif hasattr(arg, '__getitem__'):
            if callable(arg[idx]):
                return arg[idx]
            return _func_linear(self.offsets[idx], arg[idx])
        elif callable(arg):
            return _func_multadd(arg, None, self.offsets[idx])
        return _func_linear(self.offsets[idx], self.offsets[idx] + arg)

    def _parse_width(self, arg, idx):
        if arg is None or self.gdsii_path:
            if arg is not None:
                warnings.warn("[GDSPY] Argument `width` ignored in RobustPath with `gdsii_path == True`.", stacklevel=3)
            return _func_const(self.widths[idx])
        elif hasattr(arg, '__getitem__'):
            if callable(arg[idx]):
                return arg[idx]
            return _func_linear(self.widths[idx], arg[idx])
        elif callable(arg):
            return arg
        return _func_linear(self.widths[idx], arg)

    def segment(self, end_point, width=None, offset=None, relative=False):
        """
        Add a straight section to the path.

        Parameters
        ----------
        end_point : array-like[2]
            End position of the straight segment.
        width : number, callable, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  If this is callable, it must be a
            function of one argument (that varies from 0 to 1) that
            returns the width of the path.  A list can be used where
            each element (number or callable) defines the width for one
            of the parallel paths in this object.
        offset : number, callable, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  If this
            is callable, it must be a function of one argument (that
            varies from 0 to 1) that returns the offset *increase*.  A
            list can be used where each element (number or callable)
            defines the *absolute* offset (not offset increase) for one
            of the parallel paths in this object.
        relative : bool
            If True, `end_point` is used as an offset from the current
            path position, i.e., if the path is at (1, -2) and the
            `end_point` is (10, 25), the segment will be constructed
            from (1, -2) to (1 + 10, -2 + 25) = (11, 23).  Otherwise,
            `end_point` is used as an absolute coordinate.

        Returns
        -------
        out : `RobustPath`
            This object.
        """
        self._polygon_dict = None
        x = (numpy.array(end_point) + self.x) if relative else numpy.array(end_point)
        f = _func_linear(self.x, x)
        df = _func_const(x - self.x, 2)
        self.x = x
        for i in range(self.n):
            off = self._parse_offset(offset, i)
            wid = self._parse_width(width, i)
            self.paths[i].append(_SubPath(f, df, off, wid, self.tolerance, self.max_evals))
            self.widths[i] = wid(1)
            self.offsets[i] = off(1)
        return self

    def arc(self, radius, initial_angle, final_angle, width=None, offset=None):
        """
        Add a circular arc section to the path.

        Parameters
        ----------
        radius : number
            Radius of the circular arc.
        initial_angle : number
            Initial angle of the arc.
        final_angle : number
            Final angle of the arc.
        width : number, callable, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  If this is callable, it must be a
            function of one argument (that varies from 0 to 1) that
            returns the width of the path.  A list can be used where
            each element (number or callable) defines the width for one
            of the parallel paths in this object.
        offset : number, callable, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  If this
            is callable, it must be a function of one argument (that
            varies from 0 to 1) that returns the offset *increase*.  A
            list can be used where each element (number or callable)
            defines the *absolute* offset (not offset increase) for one
            of the parallel paths in this object.

        Returns
        -------
        out : `RobustPath`
            This object.
        """
        self._polygon_dict = None
        x0 = self.x - numpy.array((radius * numpy.cos(initial_angle), radius * numpy.sin(initial_angle)))
        def f(u):
            angle = initial_angle * (1 - u) + final_angle * u
            return x0 + numpy.array((radius * numpy.cos(angle), radius * numpy.sin(angle)))
        def df(u, h):
            angle = initial_angle * (1 - u) + final_angle * u
            r = radius * (final_angle - initial_angle)
            return numpy.array((-r * numpy.sin(angle), r * numpy.cos(angle)))
        self.x = f(1)
        for i in range(self.n):
            off = self._parse_offset(offset, i)
            wid = self._parse_width(width, i)
            self.paths[i].append(_SubPath(f, df, off, wid, self.tolerance, self.max_evals))
            self.widths[i] = wid(1)
            self.offsets[i] = off(1)
        return self

    def turn(self, radius, angle, width=None, offset=None):
        """
        Add a circular turn to the path.

        The initial angle of the arc is calculated from an average of
        the current directions of all parallel paths in this object.

        Parameters
        ----------
        radius : number
            Radius of the circular arc.
        angle : 'r', 'l', 'rr', 'll' or number
            Angle (in *radians*) of rotation of the path.  The values
            'r' and 'l' represent 90-degree turns cw and ccw,
            respectively; the values 'rr' and 'll' represent analogous
            180-degree turns.
        width : number, callable, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  If this is callable, it must be a
            function of one argument (that varies from 0 to 1) that
            returns the width of the path.  A list can be used where
            each element (number or callable) defines the width for one
            of the parallel paths in this object.
        offset : number, callable, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  If this
            is callable, it must be a function of one argument (that
            varies from 0 to 1) that returns the offset *increase*.  A
            list can be used where each element (number or callable)
            defines the *absolute* offset (not offset increase) for one
            of the parallel paths in this object.

        Returns
        -------
        out : `RobustPath`
            This object.
        """
        self._polygon_dict = None
        i = len(self.paths[0]) - 1
        if i < 0:
            raise ValueError("[GDSPY] Cannot define initial angle for turn on an empty RobustPath.")
        angle = _angle_dic.get(angle, angle)
        initial_angle = 0
        for p in self.paths:
            v = p[i].grad(1, 0)
            initial_angle += numpy.arctan2(v[1], v[0])
        initial_angle = initial_angle / len(self.paths) + (_halfpi if angle < 0 else -_halfpi)
        self.arc(radius, initial_angle, initial_angle + angle, width, offset)
        return self

    def parametric(self, curve_function, curve_derivative=None, width=None, offset=None, relative=True):
        """
        Add a parametric curve to the path.

        Parameters
        ----------
        curve_function : callable
            Function that defines the curve.  Must be a function of one
            argument (that varies from 0 to 1) that returns a 2-element
            Numpy array with the coordinates of the curve.
        curve_derivative : callable
            If set, it should be the derivative of the curve function.
            Must be a function of one argument (that varies from 0 to 1)
            that returns a 2-element Numpy array.  If None, the
            derivative will be calculated numerically.
        width : number, callable, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  If this is callable, it must be a
            function of one argument (that varies from 0 to 1) that
            returns the width of the path.  A list can be used where
            each element (number or callable) defines the width for one
            of the parallel paths in this object.
        offset : number, callable, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  If this
            is callable, it must be a function of one argument (that
            varies from 0 to 1) that returns the offset *increase*.  A
            list can be used where each element (number or callable)
            defines the *absolute* offset (not offset increase) for one
            of the parallel paths in this object.
        relative : bool
            If True, the return values of `curve_function` are used as
            offsets from the current path position, i.e., to ensure a
            continuous path, ``curve_function(0)`` must be (0, 0).
            Otherwise, they are used as absolute coordinates.

        Returns
        -------
        out : `RobustPath`
            This object.
        """
        self._polygon_dict = None
        f = _func_multadd(curve_function, None, self.x) if relative else curve_function
        if curve_derivative is None:
            def df(u, h):
                u0 = max(0, u - h)
                u1 = min(1, u + h)
                return (curve_function(u1) - curve_function(u0)) / (u1 - u0)
        else:
            def df(u, h):
                return curve_derivative(u)
        self.x = f(1)
        for i in range(self.n):
            off = self._parse_offset(offset, i)
            wid = self._parse_width(width, i)
            self.paths[i].append(_SubPath(f, df, off, wid, self.tolerance, self.max_evals))
            self.widths[i] = wid(1)
            self.offsets[i] = off(1)
        return self


    def bezier(self, points, width=None, offset=None, relative=True):
        """
        Add a Bezier curve to the path.

        A Bezier curve is added to the path starting from its current
        position and finishing at the last point in the `points` array.

        Parameters
        ----------
        points : array-like[N][2]
            Control points defining the Bezier curve.
        width : number, callable, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  If this is callable, it must be a
            function of one argument (that varies from 0 to 1) that
            returns the width of the path.  A list can be used where
            each element (number or callable) defines the width for one
            of the parallel paths in this object.
        offset : number, callable, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  If this
            is callable, it must be a function of one argument (that
            varies from 0 to 1) that returns the offset *increase*.  A
            list can be used where each element (number or callable)
            defines the *absolute* offset (not offset increase) for one
            of the parallel paths in this object.
        relative : bool
            If True, all coordinates in the `points` array are used as
            offsets from the current path position, i.e., if the path is
            at (1, -2) and the last point in the array is (10, 25), the
            constructed Bezier will end at (1 + 10, -2 + 25) = (11, 23).
            Otherwise, the points are used as absolute coordinates.

        Returns
        -------
        out : `RobustPath`
            This object.
        """
        self._polygon_dict = None
        if relative:
            ctrl = self.x + numpy.vstack(([(0, 0)], points))
        else:
            ctrl = numpy.vstack(([self.x], points))
        dctrl = (ctrl.shape[0] - 1) * (ctrl[1:] - ctrl[:-1])
        self.x = ctrl[-1]
        f = _func_bezier(ctrl)
        df = _func_bezier(dctrl, 2)
        for i in range(self.n):
            off = self._parse_offset(offset, i)
            wid = self._parse_width(width, i)
            self.paths[i].append(_SubPath(f, df, off, wid, self.tolerance, self.max_evals))
            self.widths[i] = wid(1)
            self.offsets[i] = off(1)
        return self


    def smooth(self, points, angles=None, curl_start=1, curl_end=1, t_in=1, t_out=1, cycle=False,
               width=None, offset=None, relative=True):
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
        width : number, callable, list
            If a number, all parallel paths are linearly tapered to this
            width along the segment.  If this is callable, it must be a
            function of one argument (that varies from 0 to 1) that
            returns the width of the path.  A list can be used where
            each element (number or callable) defines the width for one
            of the parallel paths in this object.
        offset : number, callable, list
            If a number, all parallel paths offsets are linearly
            *increased* by this amount (which can be negative).  If this
            is callable, it must be a function of one argument (that
            varies from 0 to 1) that returns the offset *increase*.  A
            list can be used where each element (number or callable)
            defines the *absolute* offset (not offset increase) for one
            of the parallel paths in this object.
        relative : bool
            If True, all coordinates in the `points` array are used as
            offsets from the current path position, i.e., if the path is
            at (1, -2) and the last point in the array is (10, 25), the
            constructed curve will end at (1 + 10, -2 + 25) = (11, 23).
            Otherwise, the points are used as absolute coordinates.

        Returns
        -------
        out : `RobustPath`
            This object.

        Notes
        -----
        Arguments `width` and `offset` are repeated for *each* cubic
        Bezier that composes this path element.

        References
        ----------
        .. [1] Hobby, J.D.  *Discrete Comput. Geom.* (1986) 1: 123.
           `DOI: 10.1007/BF02187690
           <https://doi.org/10.1007/BF02187690>`_
        """
        if relative:
            points = self.x + numpy.vstack(([(0, 0)], points))
        else:
            points = numpy.vstack(([self.x], points))
        cta, ctb = _hobby(points, angles, curl_start, curl_end, t_in, t_out, cycle)
        for i in range(points.shape[0] - 1):
            self.bezier((cta[i], ctb[i], points[i + 1]), width, offset, False)
        if cycle:
            self.bezier((cta[-1], ctb[-1], points[0]), width, offset, False)
        return self


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
        If True, the label is reflected parallel to the x direction
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
        If True, the label is reflected parallel to the x direction
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
        Translate this label.

        Parameters
        ----------
        dx : number
            Distance to move in the x-direction
        dy : number
            Distance to move in the y-direction

        Returns
        -------
        out : `Label`
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
    Collection of polygons, paths, labels and raferences to other cells.

    Parameters
    ----------
    name : string
        The name of the cell.
    exclude_from_current : bool
        If True, the cell will not be automatically included in the
        current library.

    Attributes
    ----------
    name : string
        The name of this cell.
    polygons : list of `PolygonSet`
        List of cell polygons.
    paths : list of `RobustPath` or `FlexPath`
        List of cell paths.
    labels : list of `Label`
        List of cell labels.
    references : list of `CellReference` or `CellArray`
        List of cell references.
    """
    __slots__ = 'name', 'polygons', 'paths', 'labels', 'references', '_bb_valid'

    def __init__(self, name, exclude_from_current=False):
        self.name = name
        self.polygons = []
        self.paths = []
        self.labels = []
        self.references = []
        self._bb_valid = False
        if not exclude_from_current:
            current_library.add(self)

    def __str__(self):
        return "Cell (\"{}\", {} polygons, {} paths, {} labels, {} references)".format(self.name, len(self.polygons), len(self.paths), len(self.labels), len(self.references))

    def to_gds(self, multiplier, timestamp=None):
        """
        Convert this cell to a GDSII structure.

        Parameters
        ----------
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            structure.
        timestamp : datetime object
            Sets the GDSII timestamp.  If None, the current time is
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
        data = [struct.pack('>2H12h2H', 28, 0x0502, now.year, now.month, now.day, now.hour, now.minute,
                            now.second, now.year, now.month, now.day, now.hour, now.minute, now.second,
                            4 + len(name), 0x0606), name.encode('ascii')]
        data.extend(polygon.to_gds(multiplier) for polygon in self.polygons)
        data.extend(path.to_gds(multiplier) for path in self.paths)
        data.extend(label.to_gds(multiplier) for label in self.labels)
        data.extend(reference.to_gds(multiplier) for reference in self.references)
        data.append(struct.pack('>2H', 4, 0x0700))
        return b''.join(data)

    def copy(self, name, exclude_from_current=False, deep_copy=False):
        """
        Creates a copy of this cell.

        Parameters
        ----------
        name : string
            The name of the cell.
        exclude_from_current : bool
            If True, the cell will not be included in the global list of
            cells maintained by `gdspy`.
        deep_copy : bool
            If False, the new cell will contain only references to the
            existing elements.  If True, copies of all elements are also
            created.

        Returns
        -------
        out : `Cell`
            The new copy of this cell.
        """
        new_cell = Cell(name, exclude_from_current)
        if deep_copy:
            new_cell.polygons = libcopy.deepcopy(self.polygons)
            new_cell.paths = libcopy.deepcopy(self.paths)
            new_cell.labels = libcopy.deepcopy(self.labels)
            new_cell.references = libcopy.deepcopy(self.references)
            for ref in new_cell.get_dependencies(True):
                if ref._bb_valid:
                    ref._bb_valid = False
        else:
            new_cell.polygons = list(self.polygons)
            new_cell.paths = list(self.paths)
            new_cell.labels = list(self.labels)
            new_cell.references = list(self.references)
        return new_cell

    def add(self, element):
        """
        Add a new element or list of elements to this cell.

        Parameters
        ----------
        element : `PolygonSet`, `CellReference`, `CellArray` or iterable
            The element or iterable of elements to be inserted in this
            cell.

        Returns
        -------
        out : `Cell`
            This cell.
        """
        if isinstance(element, PolygonSet):
            self.polygons.append(element)
        elif isinstance(element, RobustPath) or isinstance(element, FlexPath):
            self.paths.append(element)
        elif isinstance(element, Label):
            self.labels.append(element)
        elif isinstance(element, CellReference) or isinstance(element, CellArray):
            self.references.append(element)
        else:
            for e in element:
                if isinstance(e, PolygonSet):
                    self.polygons.append(e)
                elif isinstance(e, RobustPath) or isinstance(e, FlexPath):
                    self.paths.append(e)
                elif isinstance(e, Label):
                    self.labels.append(e)
                elif isinstance(e, CellReference) or isinstance(e, CellArray):
                    self.references.append(e)
                else:
                    raise ValueError("[GDSPY] Only instances of `PolygonSet`, `FlexPath`, `RobustPath`, `Label`, `CellReference`, and `CellArray` can be added to `Cell`.")
        self._bb_valid = False
        return self

    def remove_polygons(self, test):
        """
        Remove polygons from this cell.

        The function or callable `test` is called for each polygon in
        the cell.  If its return value evaluates to True, the
        corresponding polygon is removed from the cell.

        Parameters
        ----------
        test : callable
            Test function to query whether a polygon should be removed.
            The function is called with arguments:
            ``(points, layer, datatype)``

        Returns
        -------
        out : `Cell`
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
        for element in self.polygons:
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
            self.polygons.remove(element)
        return self

    def remove_paths(self, test):
        """
        Remove paths from this cell.

        The function or callable `test` is called for each `FlexPath`
        or `RobustPath` in the cell.  If its return value evaluates to
        True, the corresponding label is removed from the cell.

        Parameters
        ----------
        test : callable
            Test function to query whether a path should be removed.
            The function is called with the path as the only argument.

        Returns
        -------
        out : `Cell`
            This cell.
        """
        ii = 0
        while ii < len(self.paths):
            if test(self.paths[ii]):
                self.paths.pop(ii)
            else:
                ii += 1
        return self

    def remove_labels(self, test):
        """
        Remove labels from this cell.

        The function or callable `test` is called for each label in
        the cell.  If its return value evaluates to True, the
        corresponding label is removed from the cell.

        Parameters
        ----------
        test : callable
            Test function to query whether a label should be removed.
            The function is called with the label as the only argument.

        Returns
        -------
        out : `Cell`
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
            If True, the return value is a dictionary with the areas
            of each individual pair (layer, datatype).

        Returns
        -------
        out : number, dictionary
            Area of this cell.
        """
        if by_spec:
            cell_area = {}
            for element in itertools.chain(self.polygons, self.paths, self.references):
                element_area = element.area(True)
                for ll in element_area.keys():
                    if ll in cell_area:
                        cell_area[ll] += element_area[ll]
                    else:
                        cell_area[ll] = element_area[ll]
        else:
            cell_area = 0
            for element in itertools.chain(self.polygons, self.paths, self.references):
                cell_area += element.area()
        return cell_area

    def get_layers(self):
        """
        Return the set of layers in this cell.

        Returns
        -------
        out : set
            Set of the layers used in this cell.
        """
        layers = set()
        for element in itertools.chain(self.polygons, self.paths):
            layers.update(element.layers)
        for reference in self.references:
            layers.update(reference.ref_cell.get_layers())
        for label in self.labels:
            layers.add(label.layer)
        return layers

    def get_datatypes(self):
        """
        Return the set of datatypes in this cell.

        Returns
        -------
        out : set
            Set of the datatypes used in this cell.
        """
        datatypes = set()
        for element in itertools.chain(self.polygons, self.paths):
            datatypes.update(element.datatypes)
        for reference in self.references:
            datatypes.update(reference.ref_cell.get_datatypes())
        return datatypes

    def get_bounding_box(self):
        """
        Calculate the bounding box for this cell.

        Returns
        -------
        out : Numpy array[2, 2] or None
            Bounding box of this cell [[x_min, y_min], [x_max, y_max]],
            or None if the cell is empty.
        """
        if len(self.polygons) == 0 and len(self.paths) == 0 and len(self.references) == 0:
            return None
        if not (self._bb_valid and all(ref._bb_valid for ref in self.get_dependencies(True))):
            bb = numpy.array(((1e300, 1e300), (-1e300, -1e300)))
            all_polygons = []
            for polygon in self.polygons:
                all_polygons.extend(polygon.polygons)
            for path in self.paths:
                all_polygons.extend(path.to_polygonset().polygons)
            for reference in self.references:
                reference_bb = reference.get_bounding_box()
                if reference_bb is not None:
                    bb[0, 0] = min(bb[0, 0], reference_bb[0, 0])
                    bb[0, 1] = min(bb[0, 1], reference_bb[0, 1])
                    bb[1, 0] = max(bb[1, 0], reference_bb[1, 0])
                    bb[1, 1] = max(bb[1, 1], reference_bb[1, 1])
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
        Return a list of polygons in this cell.

        Parameters
        ----------
        by_spec : bool
            If True, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype).
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve polygons.  References below this level will result
            in a bounding box.  If `by_spec` is True the key will be the
            name of this cell.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons (if
            `by_spec` is True).

        Note
        ----
        Instances of `FlexPath` and `RobustPath` are also included in
        the result by computing their polygonal boundary.
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
                for polyset in self.polygons:
                    for ii in range(len(polyset.polygons)):
                        key = (polyset.layers[ii], polyset.datatypes[ii])
                        if key in polygons:
                            polygons[key].append(numpy.array(polyset.polygons[ii]))
                        else:
                            polygons[key] = [numpy.array(polyset.polygons[ii])]
                for path in self.paths:
                    path_polygons = path.get_polygons(True)
                    for kk in path_polygons.keys():
                        if kk in polygons:
                            polygons[kk].extend(path_polygons[kk])
                        else:
                            polygons[kk] = path_polygons[kk]
                for reference in self.references:
                    if depth is None:
                        next_depth = None
                    else:
                        next_depth = depth - 1
                    cell_polygons = reference.get_polygons(True, next_depth)
                    for kk in cell_polygons.keys():
                        if kk in polygons:
                            polygons[kk].extend(cell_polygons[kk])
                        else:
                            polygons[kk] = cell_polygons[kk]
            else:
                polygons = []
                for polyset in self.polygons:
                    for points in polyset.polygons:
                        polygons.append(numpy.array(points))
                for path in self.paths:
                    polygons.extend(path.get_polygons())
                for reference in self.references:
                    if depth is None:
                        next_depth = None
                    else:
                        next_depth = depth - 1
                    polygons.extend(reference.get_polygons(depth=next_depth))
        return polygons

    def get_polygonsets(self, depth=None):
        """
        Return a list with a copy of the polygons in this cell.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve polygons from.

        Returns
        -------
        out : list of `PolygonSet`
            List containing the polygons in this cell and its
            references.
        """
        polys = libcopy.deepcopy(self.polygons)
        if depth is None or depth > 0:
            for reference in self.references:
                if depth is None:
                    next_depth = None
                else:
                    next_depth = depth - 1
                polys.extend(reference.get_polygonsets(next_depth))
        return polys

    def get_paths(self, depth=None):
        """
        Return a list with a copy of the paths in this cell.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve paths from.

        Returns
        -------
        out : list of `FlexPath` or `RobustPath`
            List containing the paths in this cell and its references.
        """
        paths = libcopy.deepcopy(self.paths)
        if depth is None or depth > 0:
            for reference in self.references:
                if depth is None:
                    next_depth = None
                else:
                    next_depth = depth - 1
                paths.extend(reference.get_paths(next_depth))
        return paths

    def get_labels(self, depth=None):
        """
        Return a list with a copy of the labels in this cell.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve labels from.

        Returns
        -------
        out : list of `Label`
            List containing the labels in this cell and its references.
        """
        labels = libcopy.deepcopy(self.labels)
        if depth is None or depth > 0:
            for reference in self.references:
                if depth is None:
                    next_depth = None
                else:
                    next_depth = depth - 1
                labels.extend(reference.get_labels(next_depth))
        return labels

    def get_dependencies(self, recursive=False):
        """
        Return a list of the cells included in this cell as references.

        Parameters
        ----------
        recursive : bool
            If True returns cascading dependencies.

        Returns
        -------
        out : set of `Cell`
            List of the cells referenced by this cell.
        """
        dependencies = set()
        for reference in self.references:
            if recursive:
                dependencies.update(reference.ref_cell.get_dependencies(True))
            dependencies.add(reference.ref_cell)
        return dependencies

    def flatten(self, single_layer=None, single_datatype=None, single_texttype=None):
        """
        Convert all references into polygons, paths and labels.

        Parameters
        ----------
        single_layer : integer or None
            If not None, all polygons will be transfered to the
            layer indicated by this number.
        single_datatype : integer or None
            If not None, all polygons will be transfered to the
            datatype indicated by this number.
        single_datatype : integer or None
            If not None, all labels will be transfered to the
            texttype indicated by this number.

        Returns
        -------
        out : `Cell`
            This cell.
        """
        self.labels = self.get_labels()
        if single_layer is not None and single_datatype is not None:
            for lbl in self.labels:
                lbl.layer = single_layer
                lbl.texttype = single_texttype
        elif single_layer is not None:
            for lbl in self.labels:
                lbl.layer = single_layer
        elif single_datatype is not None:
            for lbl in self.labels:
                lbl.texttype = single_texttype
        self.polygons = self.get_polygonsets()
        self.paths = self.get_paths()
        if single_layer is not None and single_datatype is not None:
            for poly in self.polygons:
                poly.layers = [single_layer] * len(poly.polygons)
                poly.datatypes = [single_datatype] * len(poly.polygons)
            for path in self.paths:
                path.layers = [single_layer] * path.n
                path.datatypes = [single_datatype] * path.n
        elif single_layer is not None:
            for poly in self.polygons:
                poly.layers = [single_layer] * len(poly.polygons)
            for path in self.paths:
                path.layers = [single_layer] * path.n
        elif single_datatype is not None:
            for poly in self.polygons:
                poly.datatypes = [single_datatype] * len(poly.polygons)
            for path in self.paths:
                path.datatypes = [single_datatype] * path.n
        self.references = []
        return self


class CellReference(object):
    """
    Simple reference to an existing cell.

    Parameters
    ----------
    ref_cell : `Cell` or string
        The referenced cell or its name.
    origin : array-like[2]
        Position where the reference is inserted.
    rotation : number
        Angle of rotation of the reference (in *degrees*).
    magnification : number
        Magnification factor for the reference.
    x_reflection : bool
        If True the reference is reflected parallel to the x
        direction before being rotated.
    ignore_missing : bool
        If False a warning is issued when the referenced cell is not
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
            If True, the return value is a dictionary with the areas
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
        Return the list of polygons created by this reference.

        Parameters
        ----------
        by_spec : bool
            If True, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype).
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve polygons.  References below this level will result
            in a bounding box.  If `by_spec` is True the key will be the
            name of the referenced cell.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons (if
            `by_spec` is True).

        Note
        ----
        Instances of `FlexPath` and `RobustPath` are also included in
        the result by computing their polygonal boundary.
        """
        if not isinstance(self.ref_cell, Cell):
            return dict() if by_spec else []
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0) * _mpone
        if self.x_reflection:
            xrefl = _pmone_int
        if self.magnification is not None:
            mag = self.magnification * _one
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

    def get_polygonsets(self, depth=None):
        """
        Return the list of polygons created by this reference.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve polygons from.

        Returns
        -------
        out : list of `PolygonSet`
            List containing the polygons in this cell and its
            references.
        """
        if not isinstance(self.ref_cell, Cell):
            return []
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0) * _mpone
        if self.x_reflection:
            xrefl = _pmone_int
        if self.magnification is not None:
            mag = self.magnification * _one
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        polygonsets = self.ref_cell.get_polygonsets(depth=depth)
        for ps in polygonsets:
            for ii in range(len(ps.polygons)):
                if self.x_reflection:
                    ps.polygons[ii] = ps.polygons[ii] * xrefl
                if self.magnification is not None:
                    ps.polygons[ii] = ps.polygons[ii] * mag
                if self.rotation is not None:
                    ps.polygons[ii] = (ps.polygons[ii] * ct + ps.polygons[ii][:, ::-1] * st)
                if self.origin is not None:
                    ps.polygons[ii] = ps.polygons[ii] + orgn
        return polygonsets

    def get_paths(self, depth=None):
        """
        Return the list of paths created by this reference.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve paths from.

        Returns
        -------
        out : list of `FlexPath` or `RobustPath`
            List containing the paths in this cell and its references.
        """
        if not isinstance(self.ref_cell, Cell):
            return []
        if self.origin is not None:
            trans = numpy.array(self.origin)
        else:
            trans = None
        if self.rotation is not None:
            rot = self.rotation * numpy.pi / 180.0
        else:
            rot = None
        return [p.transform(trans, rot, self.magnification, self.x_reflection)
                for p in self.ref_cell.get_paths(depth=depth)]

    def get_labels(self, depth=None):
        """
        Return the list of labels created by this reference.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve labels from.

        Returns
        -------
        out : list of `Label`
            List containing the labels in this cell and its references.
        """
        if not isinstance(self.ref_cell, Cell):
            return []
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0) * _mpone
        if self.x_reflection:
            xrefl = _pmone_int
        if self.magnification is not None:
            mag = self.magnification * _one
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
        Calculate the bounding box for this reference.

        Returns
        -------
        out : Numpy array[2, 2] or None
            Bounding box of this cell [[x_min, y_min], [x_max, y_max]],
            or None if the cell is empty.
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
        Translate this reference.

        Parameters
        ----------
        dx : number
            Distance to move in the x-direction.
        dy : number
            Distance to move in the y-direction.

        Returns
        -------
        out : `CellReference`
            This object.
        """
        self.origin = (self.origin[0] + dx, self.origin[1] + dy)
        return self


class CellArray(object):
    """
    Multiple references to an existing cell in an array format.

    Parameters
    ----------
    ref_cell : `Cell` or string
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
        If True, the reference is reflected parallel to the x
        direction before being rotated.
    ignore_missing : bool
        If False a warning is issued when the referenced cell is not
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
            If True, the return value is a dictionary with the areas
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
        Return the list of polygons created by this reference.

        Parameters
        ----------
        by_spec : bool
            If True, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype).
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve polygons.  References below this level will result
            in a bounding box.  If `by_spec` is True the key will be
            name of the referenced cell.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with the list of polygons (if
            `by_spec` is True).

        Note
        ----
        Instances of `FlexPath` and `RobustPath` are also included in
        the result by computing their polygonal boundary.
        """
        if not isinstance(self.ref_cell, Cell):
            return dict() if by_spec else []
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0) * _mpone
        if self.magnification is not None:
            mag = self.magnification * _one
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        if self.x_reflection:
            xrefl = _pmone_int
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
                        if self.magnification is not None:
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

    def get_polygonsets(self, depth=None):
        """
        Return the list of polygons created by this reference.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve polygons from.

        Returns
        -------
        out : list of `PolygonSet`
            List containing the polygons in this cell and its
            references.
        """
        if not isinstance(self.ref_cell, Cell):
            return []
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0) * _mpone
        if self.x_reflection:
            xrefl = _pmone_int
        if self.magnification is not None:
            mag = self.magnification * _one
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        polygonsets = self.ref_cell.get_polygonsets(depth=depth)
        array = []
        for i in range(self.columns):
            for j in range(self.rows):
                spc = numpy.array([self.spacing[0] * i, self.spacing[1] * j])
                for polygonset in polygonsets:
                    ps = libcopy.deepcopy(polygonset)
                    for ii in range(len(ps.polygons)):
                        if self.magnification is not None:
                            ps.polygons[ii] = ps.polygons[ii] * mag + spc
                        else:
                            ps.polygons[ii] = ps.polygons[ii] + spc
                        if self.x_reflection:
                            ps.polygons[ii] = ps.polygons[ii] * xrefl
                        if self.rotation is not None:
                            ps.polygons[ii] = (ps.polygons[ii] * ct + ps.polygons[ii][:, ::-1] * st)
                        if self.origin is not None:
                            ps.polygons[ii] = ps.polygons[ii] + orgn
                    array.append(ps)
        return array

    def get_paths(self, depth=None):
        """
        Return the list of paths created by this reference.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve paths from.

        Returns
        -------
        out : list of `FlexPath` or `RobustPath`
            List containing the paths in this cell and its references.
        """
        if not isinstance(self.ref_cell, Cell):
            return []
        if self.origin is not None:
            trans = numpy.array(self.origin)
        else:
            trans = None
        if self.rotation is not None:
            rot = self.rotation * numpy.pi / 180.0
        else:
            rot = None
        paths = self.ref_cell.get_paths(depth=depth)
        array = []
        for i in range(self.columns):
            for j in range(self.rows):
                spc = numpy.array([self.spacing[0] * i, self.spacing[1] * j])
                for path in paths:
                    array.append(libcopy.deepcopy(path).transform(trans, rot, self.magnification,
                                                                  self.x_reflection, spc))
        return array

    def get_labels(self, depth=None):
        """
        Return the list of labels created by this reference.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve labels from.

        Returns
        -------
        out : list of `Label`
            List containing the labels in this cell and its references.
        """
        if not isinstance(self.ref_cell, Cell):
            return []
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0) * _mpone
        if self.magnification is not None:
            mag = self.magnification * _one
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        if self.x_reflection:
            xrefl = _pmone_int
        cell_labels = self.ref_cell.get_labels(depth=depth)
        labels = []
        for ii in range(self.columns):
            for jj in range(self.rows):
                spc = numpy.array([self.spacing[0] * ii, self.spacing[1] * jj])
                for clbl in cell_labels:
                    lbl = libcopy.deepcopy(clbl)
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
        Calculate the bounding box for this reference.

        Returns
        -------
        out : Numpy array[2, 2] or None
            Bounding box of this cell [[x_min, y_min], [x_max, y_max]],
            or None if the cell is empty.
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
        Translate this reference.

        Parameters
        ----------
        dx : number
            Distance to move in the x-direction.
        dy : number
            Distance to move in the y-direction.

        Returns
        -------
        out : `CellArray`
            This object.
        """
        self.origin = (self.origin[0] + dx, self.origin[1] + dy)
        return self


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
    _pathtype_dict = {0: 'flush', 1: 'round', 2: 'extended'}

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
        cell : `Cell` or iterable
            Cells to be included in the library.
        overwrite_duplicate : bool
            If True an existing cell with the same name in the library
            will be overwritten.

        Returns
        -------
        out : `GdsLibrary`
            This object.

        Notes
        -----
        `CellReference` or `CellArray` instances that referred to an
        overwritten cell are not automatically updated.
        """
        if isinstance(cell, Cell):
            if not overwrite_duplicate and cell.name in self.cell_dict and self.cell_dict[cell.name] is not cell:
                raise ValueError("[GDSPY] Cell named {0} already present in library.".format(cell.name))
            self.cell_dict[cell.name] = cell
        else:
            for c in cell:
                if not overwrite_duplicate and c.name in self.cell_dict and self.cell_dict[c.name] is not c:
                    raise ValueError("[GDSPY] Cell named {0} already present in library.".format(c.name))
                self.cell_dict[c.name] = c
        return self

    def write_gds(self, outfile, cells=None, timestamp=None, binary_cells=None):
        """
        Write the GDSII library to a file.

        The dimensions actually written on the GDSII file will be the
        dimensions of the objects created times the ratio
        unit/precision.  For example, if a circle with radius 1.5 is
        created and we set `GdsLibrary.unit` to 1.0e-6 (1 um) and
        `GdsLibrary.precision` to 1.0e-9` (1 nm), the radius of the
        circle will be 1.5 um and the GDSII file will contain the
        dimension 1500 nm.

        Parameters
        ----------
        outfile : file or string
            The file (or path) where the GDSII stream will be written.
            It must be opened for writing operations in binary format.
        cells : iterable
            The cells or cell names to be included in the library.  If
            None, all cells are used.
        timestamp : datetime object
            Sets the GDSII timestamp.  If None, the current time is
            used.
        binary_cells : iterable of bytes
            Iterable with binary data for GDSII cells (from
            `get_binary_cells`, for example).

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
            file.  'convert': the imported geometry is scaled to
            this library units. 'import': the unit and precision in
            this library are replaced by those from the imported file.
            'skip': the imported geometry is not scaled and units
            are not replaced; the geometry is imported in the *user
            units* of the file.
        rename : dictionary
            Dictionary used to rename the imported cells.  Keys and
            values must be strings.
        rename_template : string
            Template string used to rename the imported cells. Appiled
            only if the cell name is not in the `rename` dictionary.
            Examples: 'prefix-{name}', '{name}-suffix'
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
        out : `GdsLibrary`
            This object.

        Notes
        -----
        Not all features from the GDSII specification are currently
        supported.  A warning will be produced if any unsupported
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
                if record[1][0] < 0:
                    kwargs['width_transform'] = False
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
                kwargs['ends'] = GdsLibrary._pathtype_dict.get(record[1][0], 'extended')
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
            kwargs['ends'] = (kwargs.pop('bgnextn', 0), kwargs.pop('endextn', 0))
        kwargs['points'] = xy.reshape((xy.size // 2, 2))
        kwargs['gdsii_path'] = True
        return FlexPath(**kwargs)

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
        cell : `Cell` or string
            Cell or name of the cell to be extracted from the imported
            file.  Referenced cells will be automatically extracted as
            well.
        overwrite_duplicate : bool
            If True an existing cell with the same name in the current
            global library will be overwritten.

        Returns
        -------
        out : `Cell`
            The extracted cell.

        Notes
        -----
        `CellReference` or `CellArray` instances that referred to an
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
    dimensions of the objects created times the ratio unit/precision.
    For example, if a circle with radius 1.5 is created and we set
    `unit` to 1.0e-6 (1 um) and `precision` to 1.0e-9 (1 nm), the radius
    of the circle will be 1.5 um and the GDSII file will contain the
    dimension 1500 nm.

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
        Sets the GDSII timestamp.  If None, the current time is
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
        cell : `Cell`
            Cell to be written.
        timestamp : datetime object
            Sets the GDSII timestamp.  If None, the current time is
            used.

        Notes
        -----
        Only the specified cell is written.  Dependencies must be
        manually included.

        Returns
        -------
        out : `GdsWriter`
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
            `get_binary_cells`, for example).

        Returns
        -------
        out : `GdsWriter`
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
            result[i] = PolygonSet(result[i], layer[i % len(layer)], datatype[i % len(datatype)])
    return result


def offset(polygons, distance, join='miter', tolerance=2, precision=0.001, join_first=False,
           max_points=199, layer=0, datatype=0):
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
        if operation in ['not', 'xor']:
            if len(poly1) == 0:
                return None
            return PolygonSet(poly1, layer, datatype).fracture(max_points, precision)
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
        pts = (points, )
        sc = 0
    else:
        pts = points
        sc = 1 if short_circuit == 'any' else -1
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


def write_gds(outfile, cells=None, name='library', unit=1.0e-6, precision=1.0e-9, timestamp=None,
              binary_cells=None):
    """
    Write the current GDSII library to a file.

    The dimensions actually written on the GDSII file will be the
    dimensions of the objects created times the ratio unit/precision.
    For example, if a circle with radius 1.5 is created and we set
    `unit` to 1.0e-6 (1 um) and `precision` to 1.0e-9 (1 nm), the radius
    of the circle will be 1.5 um and the GDSII file will contain the
    dimension 1500 nm.

    Parameters
    ----------
    outfile : file or string
        The file (or path) where the GDSII stream will be written.  It
        must be opened for writing operations in binary format.
    cells : array-like
        The sequence of cells or cell names to be included in the
        library.  If None, all cells are used.
    name : string
        Name of the GDSII library.
    unit : number
        Unit size for the objects in the library (in *meters*).
    precision : number
        Precision for the dimensions of the objects in the library (in
        *meters*).
    timestamp : datetime object
        Sets the GDSII timestamp.  If None, the current time is
        used.
    binary_cells : iterable of bytes
        Iterable with binary data for GDSII cells (from
        `get_binary_cells`, for example).
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
        the methods `update` and `hexdigest` as defined in the hashlib
        module.  If None, the dafault `hashlib.sha1()` is used.

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
Current `GdsLibrary` instance for automatic creation of GDSII files.

This variable can be freely overwritten by the user with a new instance
of `GdsLibrary`.

Examples
--------
>>> gdspy.Cell('MAIN')
>>> gdspy.current_library = GdsLibrary()  # Reset current library
>>> gdspy.Cell('MAIN')  # A new MAIN cell is created without error
"""

from gdspy.curve import Curve
from gdspy import clipper
try:
    from gdspy.viewer import LayoutViewer
except ImportError as e:
    warnings.warn("[GDSPY] LayoutViewer not available: " + repr(e), category=ImportWarning, stacklevel=2)
