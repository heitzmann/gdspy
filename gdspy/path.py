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
import warnings
import itertools
import struct

from gdspy import clipper
from gdspy.hobby import _hobby
from gdspy.polygon import PolygonSet

_halfpi = 0.5 * numpy.pi
_angle_dict = {"l": _halfpi, "r": -_halfpi, "ll": numpy.pi, "rr": -numpy.pi}

_one = numpy.array((1.0, 1.0))
_mpone = numpy.array((-1.0, 1.0))
_pmone = numpy.array((1.0, -1.0))


def _func_const(c, nargs=1):
    if nargs == 1:
        return lambda u: c
    elif nargs == 2:
        return lambda u, h: c
    return lambda *args: c


def _func_linear(c0, c1):
    return lambda u: c0 * (1 - u) + c1 * u


def _func_offset(f, a):
    return lambda u: f(u) + a


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


class _SubPath(object):
    """
    Single path component.
    """

    __slots__ = "x", "dx", "off", "wid", "h", "err", "max_evals", "transform"

    def __init__(self, x, dx, off, wid, tolerance, max_evals):
        self.x = x
        self.dx = dx
        self.off = off
        self.wid = wid
        self.err = tolerance ** 2
        self.h = 0.5 / max_evals
        self.max_evals = max_evals
        self.transform = numpy.eye(3)

    def __str__(self):
        return "SubPath ({} - {})".format(self(0, 1e-6, 0), self(1, 1e-6, 0))

    def __call__(self, u, arm, transform=True):
        v = self.dx(u, self.h)[::-1] * _pmone
        v /= (v[0] ** 2 + v[1] ** 2) ** 0.5
        x = self.x(u) + self.off(u) * v
        if arm != 0:
            u0 = max(0, u - self.h)
            v = self.dx(u0, self.h)[::-1] * _pmone
            v /= (v[0] ** 2 + v[1] ** 2) ** 0.5
            x0 = self.x(u0) + self.off(u0) * v
            u1 = min(1, u + self.h)
            v = self.dx(u1, self.h)[::-1] * _pmone
            v /= (v[0] ** 2 + v[1] ** 2) ** 0.5
            x1 = self.x(u1) + self.off(u1) * v
            w = (x1 - x0)[::-1] * _pmone
            w /= (w[0] ** 2 + w[1] ** 2) ** 0.5
            if arm < 0:
                x -= 0.5 * self.wid(u) * w
            else:
                x += 0.5 * self.wid(u) * w
        if transform:
            cx = self.transform[0]
            cy = self.transform[1]
            return numpy.array(
                (
                    cx[0] * x[0] + cx[1] * x[1] + cx[2],
                    cy[0] * x[0] + cy[1] * x[1] + cy[2],
                )
            )
        return x

    def grad(self, u, arm, transform=True):
        u0 = max(0, u - self.h)
        u1 = min(1, u + self.h)
        return (self(u1, arm, transform) - self(u0, arm, transform)) / (u1 - u0)

    def width(self, u):
        d = self(u, -1) - self(u, 1)
        return (d[0] ** 2 + d[1] ** 2) ** 0.5

    def points(self, u0, u1, arm):
        u = [u0, u1]
        pts = [numpy.array(self(u[0], arm)), numpy.array(self(u[1], arm))]
        i = 1
        while i < len(pts) < self.max_evals:
            f = 0.2
            while f < 1:
                test_u = u[i - 1] * (1 - f) + u[i] * f
                test_pt = numpy.array(self(test_u, arm))
                test_err = pts[i - 1] * (1 - f) + pts[i] * f - test_pt
                if test_err[0] ** 2 + test_err[1] ** 2 > self.err:
                    u.insert(i, test_u)
                    pts.insert(i, test_pt)
                    f = 1
                    i -= 1
                else:
                    f += 0.3
            i += 1
        return pts

    def translate(self, offset):
        self.transform = numpy.matmul(
            numpy.array(
                ((1.0, 0.0, offset[0]), (0.0, 1.0, offset[1]), (0.0, 0.0, 1.0))
            ),
            self.transform,
        )
        return self

    def scale(self, sx, sy=None):
        if sy is None:
            sy = sx
        self.transform = numpy.matmul(
            numpy.array(((sx, 0.0, 0.0), (0.0, sy, 0.0), (0.0, 0.0, 1.0))),
            self.transform,
        )
        return self

    def rotate(self, angle):
        c = numpy.cos(angle)
        s = numpy.sin(angle)
        self.transform = numpy.matmul(
            numpy.array(((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0))),
            self.transform,
        )
        return self


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

    Attributes
    ----------
    properties : {integer: string} dictionary
        Properties for these elements.

    Notes
    -----
    The value of `tolerance` should not be smaller than `precision`,
    otherwise there would be wasted computational effort in calculating
    the paths.
    """

    __slots__ = (
        "n",
        "ends",
        "corners",
        "points",
        "offsets",
        "widths",
        "layers",
        "datatypes",
        "tolerance",
        "precision",
        "max_points",
        "gdsii_path",
        "width_transform",
        "bend_radius",
        "properties",
        "_polygon_dict",
    )

    _pathtype_dict = {"flush": 0, "round": 1, "extended": 2, "smooth": 1}

    def __init__(
        self,
        points,
        width,
        offset=0,
        corners="natural",
        ends="flush",
        bend_radius=None,
        tolerance=0.01,
        precision=1e-3,
        max_points=199,
        gdsii_path=False,
        width_transform=True,
        layer=0,
        datatype=0,
    ):
        self._polygon_dict = None
        if isinstance(width, list):
            self.n = len(width)
            self.widths = width
            if isinstance(offset, list):
                self.offsets = offset
            else:
                self.offsets = [
                    (i - 0.5 * (self.n - 1)) * offset for i in range(self.n)
                ]
        else:
            if isinstance(offset, list):
                self.n = len(offset)
                self.offsets = offset
            else:
                self.n = 1
                self.offsets = [offset]
            self.widths = [width] * self.n
        self.points = numpy.array(points)
        if self.points.shape == (2,):
            self.points.resize((1, 2))
        self.widths = numpy.tile(self.widths, (len(points), 1))
        self.offsets = numpy.tile(self.offsets, (len(points), 1))
        if isinstance(ends, list):
            self.ends = [ends[i % len(ends)] for i in range(self.n)]
        else:
            self.ends = [ends for _ in range(self.n)]
        if isinstance(corners, list):
            self.corners = [corners[i % len(corners)] for i in range(self.n)]
        else:
            self.corners = [corners for _ in range(self.n)]
        if isinstance(bend_radius, list):
            self.bend_radius = [
                bend_radius[i % len(bend_radius)] for i in range(self.n)
            ]
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
        self.properties = {}
        if self.gdsii_path:
            if any(end == "smooth" or callable(end) for end in self.ends):
                warnings.warn(
                    "[GDSPY] Smooth and custom end caps are not supported "
                    "in `FlexPath` with `gdsii_path == True`.",
                    stacklevel=3,
                )
            if any(
                corner != "natural" and corner != "circular bend"
                for corner in self.corners
            ):
                warnings.warn(
                    "[GDSPY] Corner specification not supported in "
                    "`FlexPath` with `gdsii_path == True`.",
                    stacklevel=3,
                )

    def __str__(self):
        if self.n > 1:
            return "FlexPath (x{}, {} segments, layers {}, datatypes {})".format(
                self.n, self.points.shape[0], self.layers, self.datatypes
            )
        else:
            return "FlexPath ({} segments, layer {}, datatype {})".format(
                self.points.shape[0], self.layers[0], self.datatypes[0]
            )

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
                if un[0] == 0 and un[1] == 0:
                    return {} if by_spec else []
                un = un[::-1] * _mpone / (un[0] ** 2 + un[1] ** 2) ** 0.5
                for kk in range(self.n):
                    end = self.ends[kk]
                    pts = numpy.array(
                        (
                            self.points[0] + un * self.offsets[0, kk],
                            self.points[1] + un * self.offsets[1, kk],
                        )
                    )
                    vn = pts[1] - pts[0]
                    vn = vn[::-1] * _mpone / (vn[0] ** 2 + vn[1] ** 2) ** 0.5
                    v = (
                        vn * (0.5 * self.widths[0, kk]),
                        vn * (0.5 * self.widths[1, kk]),
                    )
                    poly = numpy.array(
                        (pts[0] - v[0], pts[0] + v[0], pts[1] + v[1], pts[1] - v[1])
                    )
                    if end != "flush":
                        v0 = poly[3] - poly[0]
                        v1 = poly[2] - poly[1]
                        if callable(end):
                            cap0 = end(poly[0], -v0, poly[1], v1)
                            cap1 = end(poly[3], v0, poly[2], -v1)
                            poly = numpy.array(cap0[::-1] + cap1)
                        elif end == "smooth":
                            angles = [
                                numpy.arctan2(-v0[1], -v0[0]),
                                numpy.arctan2(v1[1], v1[0]),
                            ]
                            cta, ctb = _hobby(poly[:2], angles)
                            f = _func_bezier(
                                numpy.array([poly[0], cta[0], ctb[0], poly[1]])
                            )
                            tol = self.tolerance ** 2
                            uu = [0, 1]
                            fu = [f(0), f(1)]
                            iu = 1
                            while iu < len(fu):
                                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                                test_pt = f(test_u)
                                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                                if test_err[0] ** 2 + test_err[1] ** 2 > tol:
                                    uu.insert(iu, test_u)
                                    fu.insert(iu, test_pt)
                                else:
                                    iu += 1
                            cap = fu
                            cta, ctb = _hobby(poly[2:], [angles[1], angles[0]])
                            f = _func_bezier(
                                numpy.array([poly[2], cta[0], ctb[0], poly[3]])
                            )
                            tol = self.tolerance ** 2
                            uu = [0, 1]
                            fu = [f(0), f(1)]
                            iu = 1
                            while iu < len(fu):
                                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                                test_pt = f(test_u)
                                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                                if test_err[0] ** 2 + test_err[1] ** 2 > tol:
                                    uu.insert(iu, test_u)
                                    fu.insert(iu, test_pt)
                                else:
                                    iu += 1
                            cap.extend(fu)
                            poly = numpy.array(cap)
                        elif end == "round":
                            v = pts[1] - pts[0]
                            r = 0.5 * self.widths[0, kk]
                            np = max(
                                5,
                                1
                                + int(
                                    _halfpi / numpy.arccos(1 - self.tolerance / r) + 0.5
                                ),
                            )
                            ang = numpy.linspace(_halfpi, -_halfpi, np) + numpy.arctan2(
                                -v[1], -v[0]
                            )
                            poly = (
                                pts[0]
                                + r * numpy.vstack((numpy.cos(ang), numpy.sin(ang))).T
                            )
                            r = 0.5 * self.widths[1, kk]
                            np = max(
                                5,
                                1
                                + int(
                                    _halfpi / numpy.arccos(1 - self.tolerance / r) + 0.5
                                ),
                            )
                            ang = numpy.linspace(_halfpi, -_halfpi, np) + numpy.arctan2(
                                v[1], v[0]
                            )
                            poly = numpy.vstack(
                                (
                                    poly,
                                    pts[1]
                                    + r
                                    * numpy.vstack((numpy.cos(ang), numpy.sin(ang))).T,
                                )
                            )
                        else:  # 'extended'/list
                            v = pts[1] - pts[0]
                            v /= (v[0] ** 2 + v[1] ** 2) ** 0.5
                            if end == "extended":
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
                                poly = numpy.array(
                                    (
                                        poly[0],
                                        poly[0] - v0,
                                        poly[1] - v0,
                                        poly[1],
                                        poly[2],
                                        poly[2] + v1,
                                        poly[3] + v1,
                                        poly[3],
                                    )
                                )
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
                                    cuts = [
                                        pts0[int(i * len(pts0) / (ncuts + 1.0) + 0.5)]
                                        for i in range(1, ncuts + 1)
                                    ]
                                    chopped = clipper._chop(
                                        polygons[ii], cuts, 0, 1 / self.precision
                                    )
                                else:
                                    # Horizontal cuts
                                    cuts = [
                                        pts1[int(i * len(pts1) / (ncuts + 1.0) + 0.5)]
                                        for i in range(1, ncuts + 1)
                                    ]
                                    chopped = clipper._chop(
                                        polygons[ii], cuts, 1, 1 / self.precision
                                    )
                                polygons.pop(ii)
                                polygons.extend(
                                    numpy.array(x)
                                    for x in itertools.chain.from_iterable(chopped)
                                )
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
                un2 = un[:, 0] ** 2 + un[:, 1] ** 2
                if not un2.all():
                    nz = [0]
                    nz.extend(un2.nonzero()[0] + 1)
                    self.points = self.points[nz, :]
                    self.widths = self.widths[nz, :]
                    self.offsets = self.offsets[nz, :]
                    un = self.points[1:] - self.points[:-1]
                    un2 = un[:, 0] ** 2 + un[:, 1] ** 2
                un = un[:, ::-1] * _mpone / (un2 ** 0.5).reshape((un.shape[0], 1))
                for kk in range(self.n):
                    corner = self.corners[kk]
                    end = self.ends[kk]
                    if any(self.offsets[:, kk] != 0):
                        pts = numpy.empty(self.points.shape)
                        sa = self.points[:-1] + un * self.offsets[:-1, kk : kk + 1]
                        sb = self.points[1:] + un * self.offsets[1:, kk : kk + 1]
                        vn = sb - sa
                        den = vn[1:, 0] * vn[:-1, 1] - vn[1:, 1] * vn[:-1, 0]
                        idx = numpy.nonzero(
                            den ** 2
                            < 1e-12
                            * (vn[1:, 0] ** 2 + vn[1:, 1] ** 2)
                            * (vn[:-1, 0] ** 2 + vn[:-1, 1] ** 2)
                        )[0]
                        if len(idx) > 0:
                            den[idx] = 1
                        ds = sb[:-1] - sa[1:]
                        u0 = (vn[1:, 1] * ds[:, 0] - vn[1:, 0] * ds[:, 1]) / den
                        u1 = (vn[:-1, 1] * ds[:, 0] - vn[:-1, 0] * ds[:, 1]) / den
                        if any(u0 < -1) or any(u1 > 1):
                            warnings.warn(
                                "[GDSPY] Possible inconsistency found in "
                                "`FlexPath` due to sharp corner."
                            )
                        pts[1:-1] = sb[:-1] + u0.reshape((u0.shape[0], 1)) * vn[:-1]
                        if len(idx) > 0:
                            pts[idx + 1] = 0.5 * (sa[idx + 1] + sb[idx])
                        pts[0] = sa[0]
                        pts[-1] = sb[-1]
                    else:
                        pts = self.points
                    vn = pts[1:] - pts[:-1]
                    vn = (
                        vn[:, ::-1]
                        * _mpone
                        / ((vn[:, 0] ** 2 + vn[:, 1] ** 2) ** 0.5).reshape(
                            (vn.shape[0], 1)
                        )
                    )
                    arms = [[], []]
                    caps = [[], []]
                    for ii in (0, 1):
                        sign = -1 if ii == 0 else 1
                        pa = pts[:-1] + vn * (
                            sign * 0.5 * self.widths[:-1, kk : kk + 1]
                        )
                        pb = pts[1:] + vn * (sign * 0.5 * self.widths[1:, kk : kk + 1])
                        vec = pb - pa
                        caps[0].append(pa[0])
                        caps[1].append(pb[-1])
                        for jj in range(1, self.points.shape[0] - 1):
                            p0 = pb[jj - 1]
                            v0 = vec[jj - 1]
                            p1 = pa[jj]
                            v1 = vec[jj]
                            half_w = 0.5 * self.widths[jj, kk]
                            if corner == "natural":
                                v0 = v0 * (half_w / (v0[0] ** 2 + v0[1] ** 2) ** 0.5)
                                v1 = v1 * (half_w / (v1[0] ** 2 + v1[1] ** 2) ** 0.5)
                                den = v1[1] * v0[0] - v1[0] * v0[1]
                                if den ** 2 < 1e-12 * half_w ** 4:
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
                                    arms[ii].append(
                                        0.5
                                        * (p0 + min(1, u0) * v0 + p1 + max(-1, u1) * v1)
                                    )
                                else:
                                    arms[ii].append(p0 + min(1, u0) * v0)
                                    arms[ii].append(p1 + max(-1, u1) * v1)
                            elif corner == "circular bend":
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
                                    d = (
                                        self.bend_radius[kk]
                                        * numpy.tan(da)
                                        / (v0[0] ** 2 + v0[1] ** 2) ** 0.5
                                    )
                                    np = max(
                                        2,
                                        1
                                        + int(
                                            da / numpy.arccos(1 - self.tolerance / r)
                                            + 0.5
                                        ),
                                    )
                                    angles = numpy.linspace(a0, a1, np)
                                    points = (
                                        r
                                        * numpy.vstack(
                                            (numpy.cos(angles), numpy.sin(angles))
                                        ).T
                                    )
                                    arms[ii].extend(points - points[0] + p0 - d * v0)
                            elif callable(corner):
                                arms[ii].extend(
                                    corner(p0, v0, p1, v1, pts[jj], self.widths[jj, kk])
                                )
                            else:
                                den = v1[1] * v0[0] - v1[0] * v0[1]
                                lim = (
                                    1e-12
                                    * (v0[0] ** 2 + v0[1] ** 2)
                                    * (v1[0] ** 2 + v1[1] ** 2)
                                )
                                if den ** 2 < lim:
                                    u0 = u1 = 0
                                    p = 0.5 * (p0 + p1)
                                else:
                                    dx = p1[0] - p0[0]
                                    dy = p1[1] - p0[1]
                                    u0 = (v1[1] * dx - v1[0] * dy) / den
                                    u1 = (v0[1] * dx - v0[0] * dy) / den
                                    p = 0.5 * (p0 + v0 * u0 + p1 + v1 * u1)
                                if corner == "miter":
                                    arms[ii].append(p)
                                elif u0 <= 0 and u1 >= 0:
                                    arms[ii].append(p)
                                elif corner == "bevel":
                                    arms[ii].append(p0)
                                    arms[ii].append(p1)
                                elif corner == "round":
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
                                    np = max(
                                        4,
                                        1
                                        + int(
                                            0.5
                                            * abs(a1 - a0)
                                            / numpy.arccos(1 - self.tolerance / half_w)
                                            + 0.5
                                        ),
                                    )
                                    angles = numpy.linspace(a0, a1, np)
                                    arms[ii].extend(
                                        pts[jj]
                                        + half_w
                                        * numpy.vstack(
                                            (numpy.cos(angles), numpy.sin(angles))
                                        ).T
                                    )
                                elif corner == "smooth":
                                    angles = [
                                        numpy.arctan2(v0[1], v0[0]),
                                        numpy.arctan2(v1[1], v1[0]),
                                    ]
                                    bezpts = numpy.vstack((p0, p1))
                                    cta, ctb = _hobby(bezpts, angles)
                                    f = _func_bezier(
                                        numpy.array(
                                            [bezpts[0], cta[0], ctb[0], bezpts[1]]
                                        )
                                    )
                                    tol = self.tolerance ** 2
                                    uu = [0, 1]
                                    fu = [f(0), f(1)]
                                    iu = 1
                                    while iu < len(fu):
                                        test_u = 0.5 * (uu[iu - 1] + uu[iu])
                                        test_pt = f(test_u)
                                        test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                                        if test_err[0] ** 2 + test_err[1] ** 2 > tol:
                                            uu.insert(iu, test_u)
                                            fu.insert(iu, test_pt)
                                        else:
                                            iu += 1
                                    arms[ii].extend(fu)
                    if end != "flush":
                        for ii in (0, 1):
                            if callable(end):
                                vecs = [
                                    caps[ii][0] - arms[0][-ii],
                                    arms[1][-ii] - caps[ii][1],
                                ]
                                caps[ii] = end(
                                    caps[ii][0], vecs[0], caps[ii][1], vecs[1]
                                )
                            elif end == "smooth":
                                points = numpy.array(caps[ii])
                                vecs = [
                                    caps[ii][0] - arms[0][-ii],
                                    arms[1][-ii] - caps[ii][1],
                                ]
                                angles = [
                                    numpy.arctan2(vecs[0][1], vecs[0][0]),
                                    numpy.arctan2(vecs[1][1], vecs[1][0]),
                                ]
                                cta, ctb = _hobby(points, angles)
                                f = _func_bezier(
                                    numpy.array([points[0], cta[0], ctb[0], points[1]])
                                )
                                tol = self.tolerance ** 2
                                uu = [0, 1]
                                fu = [f(0), f(1)]
                                iu = 1
                                while iu < len(fu):
                                    test_u = 0.5 * (uu[iu - 1] + uu[iu])
                                    test_pt = f(test_u)
                                    test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                                    if test_err[0] ** 2 + test_err[1] ** 2 > tol:
                                        uu.insert(iu, test_u)
                                        fu.insert(iu, test_pt)
                                    else:
                                        iu += 1
                                caps[ii] = fu
                            elif end == "round":
                                v = pts[0] - pts[1] if ii == 0 else pts[-1] - pts[-2]
                                r = 0.5 * self.widths[-ii, kk]
                                np = max(
                                    5,
                                    1
                                    + int(
                                        _halfpi / numpy.arccos(1 - self.tolerance / r)
                                        + 0.5
                                    ),
                                )
                                ang = (2 * ii - 1) * numpy.linspace(
                                    -_halfpi, _halfpi, np
                                ) + numpy.arctan2(v[1], v[0])
                                caps[ii] = list(
                                    pts[-ii]
                                    + r
                                    * numpy.vstack((numpy.cos(ang), numpy.sin(ang))).T
                                )
                            else:  # 'extended'/list
                                v = pts[0] - pts[1] if ii == 0 else pts[-1] - pts[-2]
                                v = v / (v[0] ** 2 + v[1] ** 2) ** 0.5
                                w = (2 * ii - 1) * v[::-1] * _pmone
                                r = 0.5 * self.widths[-ii, kk]
                                d = r if end == "extended" else end[ii]
                                caps[ii] = [
                                    pts[-ii] + r * w,
                                    pts[-ii] + r * w + d * v,
                                    pts[-ii] - r * w + d * v,
                                    pts[-ii] - r * w,
                                ]
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
                                    cuts = [
                                        pts0[int(i * len(pts0) / (ncuts + 1.0) + 0.5)]
                                        for i in range(1, ncuts + 1)
                                    ]
                                    chopped = clipper._chop(
                                        polygons[ii], cuts, 0, 1 / self.precision
                                    )
                                else:
                                    # Horizontal cuts
                                    cuts = [
                                        pts1[int(i * len(pts1) / (ncuts + 1.0) + 0.5)]
                                        for i in range(1, ncuts + 1)
                                    ]
                                    chopped = clipper._chop(
                                        polygons[ii], cuts, 1, 1 / self.precision
                                    )
                                polygons.pop(ii)
                                polygons.extend(
                                    numpy.array(x)
                                    for x in itertools.chain.from_iterable(chopped)
                                )
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
        if self.properties is not None and len(self.properties) > 0:
            pol.properties.update(self.properties)
        return pol.fracture(self.max_points, self.precision)

    def to_gds(self, outfile, multiplier):
        """
        Convert this object to a series of GDSII elements.

        If `FlexPath.gdsii_path` is True, GDSII path elements are
        created instead of boundaries.  Such paths do not support
        variable widths, but their memeory footprint is smaller than
        full polygonal boundaries.

        Parameters
        ----------
        outfile : open file
            Output to write the GDSII.
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            elements.
        """
        if len(self.points) == 0:
            return
        if self.gdsii_path:
            sign = 1 if self.width_transform else -1
        else:
            self.to_polygonset().to_gds(outfile, multiplier)
            return
        un = self.points[1:] - self.points[:-1]
        un2 = un[:, 0] ** 2 + un[:, 1] ** 2
        if not un2.all():
            nz = [0]
            nz.extend(un2.nonzero()[0] + 1)
            self.points = self.points[nz, :]
            self.widths = self.widths[nz, :]
            self.offsets = self.offsets[nz, :]
            un = self.points[1:] - self.points[:-1]
            un2 = un[:, 0] ** 2 + un[:, 1] ** 2
        un = un[:, ::-1] * _mpone / (un2 ** 0.5).reshape((un.shape[0], 1))
        for ii in range(self.n):
            pathtype = (
                0
                if callable(self.ends[ii])
                else FlexPath._pathtype_dict.get(self.ends[ii], 4)
            )
            outfile.write(
                struct.pack(
                    ">4Hh2Hh2Hh2Hl",
                    4,
                    0x0900,
                    6,
                    0x0D02,
                    self.layers[ii],
                    6,
                    0x0E02,
                    self.datatypes[ii],
                    6,
                    0x2102,
                    pathtype,
                    8,
                    0x0F03,
                    sign * int(round(self.widths[0, ii] * multiplier)),
                )
            )
            if pathtype == 4:
                outfile.write(
                    struct.pack(
                        ">2Hl2Hl",
                        8,
                        0x3003,
                        int(round(self.ends[ii][0] * multiplier)),
                        8,
                        0x3103,
                        int(round(self.ends[ii][1] * multiplier)),
                    )
                )
            if any(self.offsets[:, ii] != 0):
                points = numpy.zeros(self.points.shape)
                sa = self.points[:-1] + un * self.offsets[:-1, ii : ii + 1]
                sb = self.points[1:] + un * self.offsets[1:, ii : ii + 1]
                vn = sb - sa
                den = vn[1:, 0] * vn[:-1, 1] - vn[1:, 1] * vn[:-1, 0]
                idx = numpy.nonzero(
                    den ** 2
                    < 1e-12
                    * (vn[1:, 0] ** 2 + vn[1:, 1] ** 2)
                    * (vn[:-1, 0] ** 2 + vn[:-1, 1] ** 2)
                )[0]
                if len(idx) > 0:
                    den[idx] = 1
                u0 = (
                    vn[1:, 1] * (sb[:-1, 0] - sa[1:, 0])
                    - vn[1:, 0] * (sb[:-1, 1] - sa[1:, 1])
                ) / den
                points[1:-1] = sb[:-1] + u0.reshape((u0.shape[0], 1)) * vn[:-1]
                if len(idx) > 0:
                    points[idx + 1] = 0.5 * (sa[idx + 1] + sb[idx])
                points[0] = sa[0]
                points[-1] = sb[-1]
            else:
                points = self.points
            if self.corners[ii] == "circular bend":
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
                        d = r * numpy.tan(da) / (v0[0] ** 2 + v0[1] ** 2) ** 0.5
                        np = max(
                            2,
                            1
                            + int(
                                da
                                / numpy.arccos(
                                    1 - self.tolerance / (r + 0.5 * self.widths[0, ii])
                                )
                                + 0.5
                            ),
                        )
                        angles = numpy.linspace(a0, a1, np)
                        bpts = (
                            r * numpy.vstack((numpy.cos(angles), numpy.sin(angles))).T
                        )
                        bends.extend(bpts - bpts[0] + p1 - d * v0)
                    p0 = p1
                    p1 = p2
                    v0 = v1
                bends.append(p1)
                points = numpy.array(bends)
            points = numpy.round(points * multiplier).astype(">i4")
            if points.shape[0] > 8191:
                warnings.warn(
                    "[GDSPY] Paths with more than 8191 are not supported "
                    "by the official GDSII specification.  This GDSII "
                    "file might not be compatible with all readers.",
                    stacklevel=4,
                )
                i0 = 0
                while i0 < points.shape[0]:
                    i1 = min(i0 + 8191, points.shape[0])
                    outfile.write(struct.pack(">2H", 4 + 8 * (i1 - i0), 0x1003))
                    outfile.write(points[i0:i1].tobytes())
                    i0 = i1
            else:
                outfile.write(struct.pack(">2H", 4 + 8 * points.shape[0], 0x1003))
                outfile.write(points.tobytes())
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
        for (l, d), polygons in self.get_polygons(True).items():
            for p in polygons:
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

        Notes
        -----
        If `width_transform` is False, the widths are not scaled.
        """
        self._polygon_dict = None
        c0 = numpy.array(center) * (1 - scale)
        self.points = self.points * scale + c0
        if self.width_transform or not self.gdsii_path:
            self.widths = self.widths * scale
        self.offsets = self.offsets * scale
        for i, end in enumerate(self.ends):
            # CustomPlus created by bgnextn and endextn
            if isinstance(end, tuple):
                self.ends[i] = tuple([e * scale for e in end])
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
        array_trans : Numpy array[2]
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
            translation = numpy.array((0.0, 0.0))
        if array_trans is None:
            array_trans = numpy.array((0.0, 0.0))
        if rotation is None:
            cos = numpy.array((1.0, 1.0))
            sin = numpy.array((0.0, 0.0))
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
        self.points = pts * cos + pts[:, ::-1] * sin + translation
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
        self.points = numpy.vstack(
            (
                self.points,
                (self.points[-1] + numpy.array(end_point)) if relative else end_point,
            )
        )
        if self.gdsii_path or width is None:
            self.widths = numpy.vstack((self.widths, self.widths[-1]))
        elif hasattr(width, "__iter__"):
            self.widths = numpy.vstack((self.widths, width))
        else:
            self.widths = numpy.vstack((self.widths, numpy.repeat(width, self.n)))
        if offset is None:
            self.offsets = numpy.vstack((self.offsets, self.offsets[-1]))
        elif hasattr(offset, "__iter__"):
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
        elif hasattr(width, "__iter__"):
            wid = numpy.array(width)
        else:
            wid = numpy.full(self.n, width)
        if offset is None:
            off = self.offsets[-1]
        elif hasattr(offset, "__iter__"):
            off = numpy.array(offset)
        else:
            off = self.offsets[-1] + offset
        rmax = radius + max(
            (self.offsets[-1] + self.widths[-1]).max(), (off + wid).max()
        )
        np = max(
            3,
            1
            + int(
                0.5
                * abs(final_angle - initial_angle)
                / numpy.arccos(1 - self.tolerance / rmax)
                + 0.5
            ),
        )
        ang = numpy.linspace(initial_angle, final_angle, np)
        pts = radius * numpy.vstack((numpy.cos(ang), numpy.sin(ang))).T
        self.points = numpy.vstack((self.points, pts[1:] + (self.points[-1] - pts[0])))
        if width is None:
            self.widths = numpy.vstack((self.widths, numpy.tile(wid, (np - 1, 1))))
        else:
            u = numpy.linspace(0, 1, np)[1:]
            self.widths = numpy.vstack(
                (self.widths, numpy.outer(1 - u, self.widths[-1]) + numpy.outer(u, wid))
            )
        if offset is None:
            self.offsets = numpy.vstack((self.offsets, numpy.tile(off, (np - 1, 1))))
        else:
            u = numpy.linspace(0, 1, np)[1:]
            self.offsets = numpy.vstack(
                (
                    self.offsets,
                    numpy.outer(1 - u, self.offsets[-1]) + numpy.outer(u, off),
                )
            )
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
            raise ValueError(
                "[GDSPY] Cannot define initial angle for turn on a "
                "FlexPath withouth previous segments."
            )
        v = self.points[-1] - self.points[-2]
        angle = _angle_dict.get(angle, angle)
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
        elif hasattr(width, "__iter__"):
            wid = numpy.array(width)
        else:
            wid = numpy.full(self.n, width)
        if offset is None:
            off = self.offsets[-1]
        elif hasattr(offset, "__iter__"):
            off = numpy.array(offset)
        else:
            off = self.offsets[-1] + offset
        tol = self.tolerance ** 2
        u = [0, 1]
        pts = [numpy.array(curve_function(0)), numpy.array(curve_function(1))]
        i = 1
        while i < len(pts):
            f = 0.2
            while f < 1:
                test_u = u[i - 1] * (1 - f) + u[i] * f
                test_pt = numpy.array(curve_function(test_u))
                test_err = pts[i - 1] * (1 - f) + pts[i] * f - test_pt
                if test_err[0] ** 2 + test_err[1] ** 2 > tol:
                    u.insert(i, test_u)
                    pts.insert(i, test_pt)
                    f = 1
                    i -= 1
                else:
                    f += 0.3
            i += 1
        pts = numpy.array(pts[1:])
        np = pts.shape[0] + 1
        self.points = numpy.vstack(
            (self.points, (pts + self.points[-1]) if relative else pts)
        )
        if width is None:
            self.widths = numpy.vstack((self.widths, numpy.tile(wid, (np - 1, 1))))
        else:
            u = numpy.linspace(0, 1, np)[1:]
            self.widths = numpy.vstack(
                (self.widths, numpy.outer(1 - u, self.widths[-1]) + numpy.outer(u, wid))
            )
        if offset is None:
            self.offsets = numpy.vstack((self.offsets, numpy.tile(off, (np - 1, 1))))
        else:
            u = numpy.linspace(0, 1, np)[1:]
            self.offsets = numpy.vstack(
                (
                    self.offsets,
                    numpy.outer(1 - u, self.offsets[-1]) + numpy.outer(u, off),
                )
            )
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

    def smooth(
        self,
        points,
        angles=None,
        curl_start=1,
        curl_end=1,
        t_in=1,
        t_out=1,
        cycle=False,
        width=None,
        offset=None,
        relative=True,
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

    Attributes
    ----------
    properties : {integer: string} dictionary
        Properties for these elements.

    Notes
    -----
    The value of `tolerance` should not be smaller than `precision`,
    otherwise there would be wasted computational effort in calculating
    the paths.
    """

    __slots__ = (
        "n",
        "ends",
        "x",
        "offsets",
        "widths",
        "paths",
        "layers",
        "datatypes",
        "tolerance",
        "precision",
        "max_points",
        "max_evals",
        "gdsii_path",
        "width_transform",
        "properties",
        "_polygon_dict",
    )

    _pathtype_dict = {"flush": 0, "round": 1, "extended": 2, "smooth": 1}

    def __init__(
        self,
        initial_point,
        width,
        offset=0,
        ends="flush",
        tolerance=0.01,
        precision=1e-3,
        max_points=199,
        max_evals=1000,
        gdsii_path=False,
        width_transform=True,
        layer=0,
        datatype=0,
    ):
        self._polygon_dict = None
        if isinstance(width, list):
            self.n = len(width)
            self.widths = width
            if isinstance(offset, list):
                self.offsets = offset
            else:
                self.offsets = [
                    (i - 0.5 * (self.n - 1)) * offset for i in range(self.n)
                ]
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
        self.properties = {}
        if self.gdsii_path and any(end == "smooth" for end in self.ends):
            warnings.warn(
                "[GDSPY] Smooth end caps not supported in `RobustPath` "
                "with `gdsii_path == True`.",
                stacklevel=3,
            )

    def __str__(self):
        if self.n > 1:
            return "RobustPath (x{}, end at ({}, {}), length {}, layers {}, datatypes {})".format(
                self.n, self.x[0], self.x[1], len(self), self.layers, self.datatypes
            )
        else:
            return (
                "RobustPath (end at ({}, {}), length {}, layer {}, datatype {})".format(
                    self.x[0], self.x[1], len(self), self.layers[0], self.datatypes[0]
                )
            )

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

    def grad(self, u, arm=0, side="-"):
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
        if u == 0 and ((i > 0 and side == "-") or i == len(self.paths[0])):
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
        return numpy.array([p[i].width(u) for p in self.paths])

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
            for path, end, layer, datatype in zip(
                self.paths, self.ends, self.layers, self.datatypes
            ):
                poly = []
                for arm in [-1, 1]:
                    if not (end == "flush"):
                        i = 0 if arm == 1 else -1
                        u = abs(i)
                        if end == "smooth":
                            v1 = -arm * path[i].grad(u, -arm)
                            v2 = arm * path[i].grad(u, arm)
                            angles = [
                                numpy.arctan2(v1[1], v1[0]),
                                numpy.arctan2(v2[1], v2[0]),
                            ]
                            points = numpy.array([path[i](u, -arm), path[i](u, arm)])
                            cta, ctb = _hobby(points, angles)
                            f = _func_bezier(
                                numpy.array([points[0], cta[0], ctb[0], points[1]])
                            )
                            tol = self.tolerance ** 2
                            uu = [0, 1]
                            fu = [f(0), f(1)]
                            iu = 1
                            while iu < len(fu) < self.max_evals:
                                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                                test_pt = f(test_u)
                                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                                if test_err[0] ** 2 + test_err[1] ** 2 > tol:
                                    uu.insert(iu, test_u)
                                    fu.insert(iu, test_pt)
                                else:
                                    iu += 1
                            poly.extend(fu[1:-1])
                        else:
                            p = path[i](u, 0)
                            v = -arm * path[i].grad(u, 0)
                            r = 0.5 * path[i].width(u)
                            if end == "round":
                                np = max(
                                    5,
                                    1
                                    + int(
                                        _halfpi / numpy.arccos(1 - self.tolerance / r)
                                        + 0.5
                                    ),
                                )
                                ang = numpy.linspace(-_halfpi, _halfpi, np)[
                                    1:-1
                                ] + numpy.arctan2(v[1], v[0])
                                endpts = (
                                    p
                                    + r
                                    * numpy.vstack((numpy.cos(ang), numpy.sin(ang))).T
                                )
                                poly.extend(endpts)
                            else:
                                v /= (v[0] ** 2 + v[1] ** 2) ** 0.5
                                w = v[::-1] * _pmone
                                d = r if end == "extended" else end[u]
                                poly.append(p + d * v + r * w)
                                poly.append(p + d * v - r * w)
                    path_arm = []
                    start = 0
                    tol = self.tolerance ** 2
                    for sub0, sub1 in zip(path[:-1], path[1:]):
                        p0 = sub0(1, arm)
                        v0 = sub0.grad(1, arm)
                        p1 = sub1(0, arm)
                        v1 = sub1.grad(0, arm)
                        den = v1[1] * v0[0] - v1[0] * v0[1]
                        lim = (
                            1e-12
                            * (v0[0] ** 2 + v0[1] ** 2)
                            * (v1[0] ** 2 + v1[1] ** 2)
                        )
                        dx = p1[0] - p0[0]
                        dy = p1[1] - p0[1]
                        if den ** 2 < lim or dx ** 2 + dy ** 2 <= tol:
                            u0 = u1 = 0
                            px = 0.5 * (p0 + p1)
                        else:
                            u0 = (v1[1] * dx - v1[0] * dy) / den
                            u1 = (v0[1] * dx - v0[0] * dy) / den
                            px = 0.5 * (p0 + v0 * u0 + p1 + v1 * u1)
                        u0 = 1 + u0
                        if u0 < 1 and u1 > 0:
                            delta = sub0(u0, arm) - sub1(u1, arm)
                            err = delta[0] ** 2 + delta[1] ** 2
                            iters = 0
                            step = 0.5
                            while err > tol:
                                iters += 1
                                if iters > self.max_evals:
                                    warnings.warn("[GDSPY] Intersection not found.")
                                    break
                                du = delta * sub0.grad(u0, arm)
                                new_u0 = min(1, max(0, u0 - step * (du[0] + du[1])))
                                du = delta * sub1.grad(u1, arm)
                                new_u1 = min(1, max(0, u1 + step * (du[0] + du[1])))
                                new_delta = sub0(new_u0, arm) - sub1(new_u1, arm)
                                new_err = new_delta[0] ** 2 + new_delta[1] ** 2
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
                                warnings.warn(
                                    "[GDSPY] RobustPath join at ({}, {}) cannot be ensured.  "
                                    "Please check the resulting polygon.".format(
                                        path_arm[-1][0], path_arm[-1][1]
                                    ),
                                    stacklevel=3,
                                )
                            start = u1
                        else:
                            if u0 <= 1:
                                path_arm.extend(sub0.points(start, u0, arm))
                                warnings.warn(
                                    "[GDSPY] RobustPath join at ({}, {}) cannot be ensured.  "
                                    "Please check the resulting polygon.".format(
                                        path_arm[-1][0], path_arm[-1][1]
                                    ),
                                    stacklevel=2,
                                )
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
                                cuts = [
                                    pts0[int(i * len(pts0) / (ncuts + 1.0) + 0.5)]
                                    for i in range(1, ncuts + 1)
                                ]
                                chopped = clipper._chop(
                                    polygons[ii], cuts, 0, 1 / self.precision
                                )
                            else:
                                # Horizontal cuts
                                cuts = [
                                    pts1[int(i * len(pts1) / (ncuts + 1.0) + 0.5)]
                                    for i in range(1, ncuts + 1)
                                ]
                                chopped = clipper._chop(
                                    polygons[ii], cuts, 1, 1 / self.precision
                                )
                            polygons.pop(ii)
                            polygons.extend(
                                numpy.array(x)
                                for x in itertools.chain.from_iterable(chopped)
                            )
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
        if self.properties is not None and len(self.properties) > 0:
            pol.properties.update(self.properties)
        return pol.fracture(self.max_points, self.precision)

    def to_gds(self, outfile, multiplier):
        """
        Convert this object to a series of GDSII elements.

        If `RobustPath.gdsii_path` is True, GDSII path elements are
        created instead of boundaries.  Such paths do not support
        variable widths, but their memeory footprint is smaller than
        full polygonal boundaries.

        Parameters
        ----------
        outfile : open file
            Output to write the GDSII.
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            elements.
        """
        if len(self.paths[0]) == 0:
            return
        if self.gdsii_path:
            sign = 1 if self.width_transform else -1
        else:
            self.to_polygonset().to_gds(outfile, multiplier)
            return
        for ii in range(self.n):
            pathtype = RobustPath._pathtype_dict.get(self.ends[ii], 4)
            outfile.write(
                struct.pack(
                    ">4Hh2Hh2Hh2Hl",
                    4,
                    0x0900,
                    6,
                    0x0D02,
                    self.layers[ii],
                    6,
                    0x0E02,
                    self.datatypes[ii],
                    6,
                    0x2102,
                    pathtype,
                    8,
                    0x0F03,
                    sign * int(round(self.widths[ii] * multiplier)),
                )
            )
            if pathtype == 4:
                outfile.write(
                    struct.pack(
                        ">2Hl2Hl",
                        8,
                        0x3003,
                        int(round(self.ends[ii][0] * multiplier)),
                        8,
                        0x3103,
                        int(round(self.ends[ii][1] * multiplier)),
                    )
                )
            points = []
            start = 0
            tol = self.tolerance ** 2
            for sub0, sub1 in zip(self.paths[ii][:-1], self.paths[ii][1:]):
                p0 = sub0(1, 0)
                v0 = sub0.grad(1, 0)
                p1 = sub1(0, 0)
                v1 = sub1.grad(0, 0)
                den = v1[1] * v0[0] - v1[0] * v0[1]
                lim = 1e-12 * (v0[0] ** 2 + v0[1] ** 2) * (v1[0] ** 2 + v1[1] ** 2)
                dx = p1[0] - p0[0]
                dy = p1[1] - p0[1]
                if den ** 2 < lim or dx ** 2 + dy ** 2 <= tol:
                    u0 = u1 = 0
                    px = 0.5 * (p0 + p1)
                else:
                    u0 = (v1[1] * dx - v1[0] * dy) / den
                    u1 = (v0[1] * dx - v0[0] * dy) / den
                    px = 0.5 * (p0 + v0 * u0 + p1 + v1 * u1)
                u0 = 1 + u0
                if u0 < 1 and u1 > 0:
                    delta = sub0(u0, 0) - sub1(u1, 0)
                    err = delta[0] ** 2 + delta[1] ** 2
                    iters = 0
                    step = 0.5
                    while err > tol:
                        iters += 1
                        if iters > self.max_evals:
                            warnings.warn("[GDSPY] Intersection not found.")
                            break
                        du = delta * sub0.grad(u0, 0)
                        new_u0 = min(1, max(0, u0 - step * (du[0] + du[1])))
                        du = delta * sub1.grad(u1, 0)
                        new_u1 = min(1, max(0, u1 + step * (du[0] + du[1])))
                        new_delta = sub0(new_u0, 0) - sub1(new_u1, 0)
                        new_err = new_delta[0] ** 2 + new_delta[1] ** 2
                        if new_err >= err:
                            step /= 2
                            continue
                        u0 = new_u0
                        u1 = new_u1
                        delta = new_delta
                        err = new_err
                    px = 0.5 * (sub0(u0, 0) + sub1(u1, 0))
                if u1 >= 0:
                    if u0 <= 1:
                        points.extend(sub0.points(start, u0, 0)[:-1])
                    else:
                        points.extend(sub0.points(start, 1, 0))
                        warnings.warn(
                            "[GDSPY] RobustPath join at ({}, {}) cannot be ensured.  "
                            "Please check the resulting polygon.".format(
                                points[-1][0], points[-1][1]
                            ),
                            stacklevel=3,
                        )
                    start = u1
                else:
                    if u0 <= 1:
                        points.extend(sub0.points(start, u0, 0))
                        warnings.warn(
                            "[GDSPY] RobustPath join at ({}, {}) cannot be ensured.  "
                            "Please check the resulting polygon.".format(
                                points[-1][0], points[-1][1]
                            ),
                            stacklevel=2,
                        )
                    else:
                        points.extend(sub0.points(start, 1, 0))
                        points.append(px)
                    start = 0
            points.extend(self.paths[ii][-1].points(start, 1, 0))
            points = (numpy.array(points) * multiplier).astype(">i4")
            if points.shape[0] > 8191:
                warnings.warn(
                    "[GDSPY] Paths with more than 8191 are not supported "
                    "by the official GDSII specification.  This GDSII "
                    "file might not be compatible with all readers.",
                    stacklevel=4,
                )
                i0 = 0
                while i0 < points.shape[0]:
                    i1 = min(i0 + 8191, points.shape[0])
                    outfile.write(struct.pack(">2H", 4 + 8 * (i1 - i0), 0x1003))
                    outfile.write(points[i0:i1].tobytes())
                    i0 = i1
            else:
                outfile.write(struct.pack(">2H", 4 + 8 * points.shape[0], 0x1003))
                outfile.write(points.tobytes())
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
        for (l, d), polygons in self.get_polygons(True).items():
            for p in polygons:
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
                sub.translate(offset)
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
                sub.translate(-c0)
                sub.rotate(angle)
                sub.translate(c0)
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

        Notes
        -----
        If `width_transform` is False, the widths are not scaled.
        """
        self._polygon_dict = None
        c0 = numpy.array(center) * (1 - scale)
        self.x = self.x * scale + c0
        if self.width_transform or not self.gdsii_path:
            self.widths = [wid * scale for wid in self.widths]
        self.offsets = [off * scale for off in self.offsets]
        for path in self.paths:
            for sub in path:
                sub.scale(scale)
                sub.translate(c0)
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
        array_trans : Numpy array[2]
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
        if array_trans is not None:
            self.x = self.x + array_trans
            for path in self.paths:
                for sub in path:
                    sub.translate(array_trans)
        if x_reflection:
            self.x = numpy.array((self.x[0], -self.x[1]))
            self.widths = [wid * -1 for wid in self.widths]
            self.offsets = [off * -1 for off in self.offsets]
            for path in self.paths:
                for sub in path:
                    sub.scale(1, -1)
        if scale is not None:
            self.x = self.x * scale
            if self.width_transform or not self.gdsii_path:
                self.widths = [wid * scale for wid in self.widths]
            self.offsets = [off * scale for off in self.offsets]
            for path in self.paths:
                for sub in path:
                    sub.scale(scale)
        if rotation is not None:
            ca = numpy.cos(rotation)
            sa = numpy.sin(rotation) * _mpone
            self.x = self.x * ca + self.x[::-1] * sa
            for path in self.paths:
                for sub in path:
                    sub.rotate(rotation)
        if translation is not None:
            self.x = self.x + translation
            for path in self.paths:
                for sub in path:
                    sub.translate(translation)
        return self

    def _parse_offset(self, arg, idx):
        if arg is None:
            return _func_const(self.offsets[idx])
        elif hasattr(arg, "__getitem__"):
            if callable(arg[idx]):
                return arg[idx]
            return _func_linear(self.offsets[idx], arg[idx])
        elif callable(arg):
            return _func_offset(arg, self.offsets[idx])
        return _func_linear(self.offsets[idx], self.offsets[idx] + arg)

    def _parse_width(self, arg, idx):
        if arg is None or self.gdsii_path:
            if arg is not None:
                warnings.warn(
                    "[GDSPY] Argument `width` ignored in RobustPath with "
                    "`gdsii_path == True`.",
                    stacklevel=3,
                )
            return _func_const(self.widths[idx])
        elif hasattr(arg, "__getitem__"):
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
        self.x = numpy.array(x)
        for i in range(self.n):
            off = self._parse_offset(offset, i)
            wid = self._parse_width(width, i)
            self.paths[i].append(
                _SubPath(f, df, off, wid, self.tolerance, self.max_evals)
            )
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
        x0 = self.x - numpy.array(
            (radius * numpy.cos(initial_angle), radius * numpy.sin(initial_angle))
        )

        def f(u):
            angle = initial_angle * (1 - u) + final_angle * u
            return x0 + numpy.array(
                (radius * numpy.cos(angle), radius * numpy.sin(angle))
            )

        def df(u, h):
            angle = initial_angle * (1 - u) + final_angle * u
            r = radius * (final_angle - initial_angle)
            return numpy.array((-r * numpy.sin(angle), r * numpy.cos(angle)))

        self.x = f(1)
        for i in range(self.n):
            off = self._parse_offset(offset, i)
            wid = self._parse_width(width, i)
            self.paths[i].append(
                _SubPath(f, df, off, wid, self.tolerance, self.max_evals)
            )
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
            raise ValueError(
                "[GDSPY] Cannot define initial angle for turn on an "
                "empty RobustPath."
            )
        angle = _angle_dict.get(angle, angle)
        initial_angle = 0
        for p in self.paths:
            v = p[i].grad(1, 0, False)
            initial_angle += numpy.arctan2(v[1], v[0])
        initial_angle = initial_angle / len(self.paths) + (
            _halfpi if angle < 0 else -_halfpi
        )
        self.arc(radius, initial_angle, initial_angle + angle, width, offset)
        return self

    def parametric(
        self,
        curve_function,
        curve_derivative=None,
        width=None,
        offset=None,
        relative=True,
    ):
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
        f = _func_offset(curve_function, self.x) if relative else curve_function
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
            self.paths[i].append(
                _SubPath(f, df, off, wid, self.tolerance, self.max_evals)
            )
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
            self.paths[i].append(
                _SubPath(f, df, off, wid, self.tolerance, self.max_evals)
            )
            self.widths[i] = wid(1)
            self.offsets[i] = off(1)
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
        width=None,
        offset=None,
        relative=True,
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
