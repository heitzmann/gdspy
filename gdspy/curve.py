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

import numpy
from gdspy.path import _func_bezier, _hobby


class Curve(object):
    """
    Generation of curves loosely based on SVG paths.

    Short summary of available methods:

    ====== =============================
    Method Primitive
    ====== =============================
    L/l    Line segments
    H/h    Horizontal line segments
    V/v    Vertical line segments
    C/c    Cubic Bezier curve
    S/s    Smooth cubic Bezier curve
    Q/q    Quadratic Bezier curve
    T/t    Smooth quadratic Bezier curve
    B/b    General degree Bezier curve
    I/i    Smooth interpolating curve
    arc    Elliptical arc
    ====== =============================

    The uppercase version of the methods considers that all coordinates
    are absolute, whereas the lowercase considers that they are relative
    to the current end point of the curve.

    Parameters
    ----------
    x : number
        X-coordinate of the starting point of the curve.  If this is a
        complex number, the value of `y` is ignored and the starting
        point becomes ``(x.real, x.imag)``.
    y : number
        Y-coordinate of the starting point of the curve.
    tolerance : number
        Tolerance used to calculate a polygonal approximation to the
        curve.

    Notes
    -----

    In all methods of this class that accept coordinate pairs, a single
    complex number can be passed to be split into its real and imaginary
    parts.
    This feature can be useful in expressing coordinates in polar form.

    All commands follow the SVG 2 specification, except for elliptical
    arcs and smooth interpolating curves, which are inspired by the
    Metapost syntax.

    Examples
    --------
    >>> curve = gdspy.Curve(3, 4).H(1).q(0.5, 1, 2j).L(2 + 3j, 2, 2)
    >>> pol = gdspy.Polygon(curve.get_points())
    """

    __slots__ = "points", "tol", "last_c", "last_q"

    def __init__(self, x, y=0, tolerance=0.01):
        self.last_c = self.last_q = None
        self.tol = tolerance ** 2
        if isinstance(x, complex):
            self.points = [numpy.array((x.real, x.imag))]
        else:
            self.points = [numpy.array((x, y))]

    def get_points(self):
        """
        Get the polygonal points that approximate this curve.

        Returns
        -------
        out : Numpy array[N, 2]
            Vertices of the polygon.
        """
        delta = (self.points[-1] - self.points[0]) ** 2
        if delta[0] + delta[1] < self.tol:
            return numpy.array(self.points[:-1])
        return numpy.array(self.points)

    def L(self, *xy):
        """
        Add straight line segments to the curve.

        Parameters
        ----------
        xy : numbers
            Endpoint coordinates of the line segments.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = self.last_q = None
        i = 0
        while i < len(xy):
            if isinstance(xy[i], complex):
                self.points.append(numpy.array((xy[i].real, xy[i].imag)))
                i += 1
            else:
                self.points.append(numpy.array((xy[i], xy[i + 1])))
                i += 2
        return self

    def l(self, *xy):
        """
        Add straight line segments to the curve.

        Parameters
        ----------
        xy : numbers
            Endpoint coordinates of the line segments relative to the
            current end point.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = self.last_q = None
        o = self.points[-1]
        i = 0
        while i < len(xy):
            if isinstance(xy[i], complex):
                self.points.append(o + numpy.array((xy[i].real, xy[i].imag)))
                i += 1
            else:
                self.points.append(o + numpy.array((xy[i], xy[i + 1])))
                i += 2
        return self

    def H(self, *x):
        """
        Add horizontal line segments to the curve.

        Parameters
        ----------
        x : numbers
            Endpoint x-coordinates of the line segments.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = self.last_q = None
        y0 = self.points[-1][1]
        self.points.extend(numpy.array((xx, y0)) for xx in x)
        return self

    def h(self, *x):
        """
        Add horizontal line segments to the curve.

        Parameters
        ----------
        x : numbers
            Endpoint x-coordinates of the line segments relative to the
            current end point.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = self.last_q = None
        x0, y0 = self.points[-1]
        self.points.extend(numpy.array((x0 + xx, y0)) for xx in x)
        return self

    def V(self, *y):
        """
        Add vertical line segments to the curve.

        Parameters
        ----------
        y : numbers
            Endpoint y-coordinates of the line segments.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = self.last_q = None
        x0 = self.points[-1][0]
        self.points.extend(numpy.array((x0, yy)) for yy in y)
        return self

    def v(self, *y):
        """
        Add vertical line segments to the curve.

        Parameters
        ----------
        y : numbers
            Endpoint y-coordinates of the line segments relative to the
            current end point.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = self.last_q = None
        x0, y0 = self.points[-1]
        self.points.extend(numpy.array((x0, y0 + yy)) for yy in y)
        return self

    def arc(self, radius, initial_angle, final_angle, rotation=0):
        """
        Add an elliptical arc to the curve.

        Parameters
        ----------
        radius : number, array-like[2]
            Arc radius.  An elliptical arc can be created by passing an
            array with 2 radii.
        initial_angle : number
            Initial angle of the arc (in *radians*).
        final_angle : number
            Final angle of the arc (in *radians*).
        rotation : number
            Rotation of the axis of the ellipse.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = self.last_q = None
        if hasattr(radius, "__iter__"):
            rx, ry = radius
            radius = max(radius)
        else:
            rx = ry = radius
        full_angle = abs(final_angle - initial_angle)
        number_of_points = max(
            3,
            1
            + int(0.5 * full_angle / numpy.arccos(1 - self.tol ** 0.5 / radius) + 0.5),
        )
        angles = numpy.linspace(
            initial_angle - rotation, final_angle - rotation, number_of_points
        )
        pts = numpy.vstack((rx * numpy.cos(angles), ry * numpy.sin(angles))).T
        if rotation != 0:
            rot = numpy.empty_like(pts)
            c = numpy.cos(rotation)
            s = numpy.sin(rotation)
            rot[:, 0] = pts[:, 0] * c - pts[:, 1] * s
            rot[:, 1] = pts[:, 0] * s + pts[:, 1] * c
        else:
            rot = pts
        pts = rot[1:] - rot[0] + self.points[-1]
        self.points.extend(xy for xy in pts)
        return self

    def C(self, *xy):
        """
        Add cubic Bezier curves to the curve.

        Parameters
        ----------
        xy : numbers
            Coordinate pairs. Each set of 3 pairs are interpreted as
            the control point at the beginning of the curve, the control
            point at the end of the curve and the endpoint of the curve.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_q = None
        i = 0
        while i < len(xy):
            ctrl = numpy.empty((4, 2))
            ctrl[0] = self.points[-1]
            for j in range(1, 4):
                if isinstance(xy[i], complex):
                    ctrl[j, 0] = xy[i].real
                    ctrl[j, 1] = xy[i].imag
                    i += 1
                else:
                    ctrl[j, 0] = xy[i]
                    ctrl[j, 1] = xy[i + 1]
                    i += 2
            f = _func_bezier(ctrl)
            uu = [0, 0.2, 0.5, 0.8, 1]
            fu = [f(u) for u in uu]
            iu = 1
            while iu < len(fu):
                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                test_pt = f(test_u)
                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                if test_err[0] ** 2 + test_err[1] ** 2 > self.tol:
                    uu.insert(iu, test_u)
                    fu.insert(iu, test_pt)
                else:
                    iu += 1
            self.points.extend(xy for xy in fu[1:])
        self.last_c = ctrl[2]
        return self

    def c(self, *xy):
        """
        Add cubic Bezier curves to the curve.

        Parameters
        ----------
        xy : numbers
            Coordinate pairs. Each set of 3 pairs are interpreted as
            the control point at the beginning of the curve, the control
            point at the end of the curve and the endpoint of the curve.
            All coordinates are relative to the current end point.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_q = None
        x0, y0 = self.points[-1]
        i = 0
        while i < len(xy):
            ctrl = numpy.empty((4, 2))
            ctrl[0] = self.points[-1]
            for j in range(1, 4):
                if isinstance(xy[i], complex):
                    ctrl[j, 0] = x0 + xy[i].real
                    ctrl[j, 1] = y0 + xy[i].imag
                    i += 1
                else:
                    ctrl[j, 0] = x0 + xy[i]
                    ctrl[j, 1] = y0 + xy[i + 1]
                    i += 2
            f = _func_bezier(ctrl)
            uu = [0, 0.2, 0.5, 0.8, 1]
            fu = [f(u) for u in uu]
            iu = 1
            while iu < len(fu):
                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                test_pt = f(test_u)
                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                if test_err[0] ** 2 + test_err[1] ** 2 > self.tol:
                    uu.insert(iu, test_u)
                    fu.insert(iu, test_pt)
                else:
                    iu += 1
            self.points.extend(xy for xy in fu[1:])
        self.last_c = ctrl[2]
        return self

    def S(self, *xy):
        """
        Add smooth cubic Bezier curves to the curve.

        Parameters
        ----------
        xy : numbers
            Coordinate pairs. Each set of 2 pairs are interpreted as
            the control point at the end of the curve and the endpoint
            of the curve.  The control point at the beginning of the
            curve is assumed to be the reflection of the control point
            at the end of the last curve relative to the starting point
            of the curve. If the previous curve is not a cubic Bezier,
            the control point is coincident with the starting point.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_q = None
        if self.last_c is None:
            self.last_c = self.points[-1]
        i = 0
        while i < len(xy):
            ctrl = numpy.empty((4, 2))
            ctrl[0] = self.points[-1]
            ctrl[1] = 2 * ctrl[0] - self.last_c
            for j in range(2, 4):
                if isinstance(xy[i], complex):
                    ctrl[j, 0] = xy[i].real
                    ctrl[j, 1] = xy[i].imag
                    i += 1
                else:
                    ctrl[j, 0] = xy[i]
                    ctrl[j, 1] = xy[i + 1]
                    i += 2
            f = _func_bezier(ctrl)
            uu = [0, 0.2, 0.5, 0.8, 1]
            fu = [f(u) for u in uu]
            iu = 1
            while iu < len(fu):
                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                test_pt = f(test_u)
                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                if test_err[0] ** 2 + test_err[1] ** 2 > self.tol:
                    uu.insert(iu, test_u)
                    fu.insert(iu, test_pt)
                else:
                    iu += 1
            self.points.extend(xy for xy in fu[1:])
            self.last_c = ctrl[2]
        return self

    def s(self, *xy):
        """
        Add smooth cubic Bezier curves to the curve.

        Parameters
        ----------
        xy : numbers
            Coordinate pairs. Each set of 2 pairs are interpreted as
            the control point at the end of the curve and the endpoint
            of the curve.  The control point at the beginning of the
            curve is assumed to be the reflection of the control point
            at the end of the last curve relative to the starting point
            of the curve. If the previous curve is not a cubic Bezier,
            the control point is coincident with the starting point.
            All coordinates are relative to the current end point.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_q = None
        if self.last_c is None:
            self.last_c = self.points[-1]
        x0, y0 = self.points[-1]
        i = 0
        while i < len(xy):
            ctrl = numpy.empty((4, 2))
            ctrl[0] = self.points[-1]
            ctrl[1] = 2 * ctrl[0] - self.last_c
            for j in range(2, 4):
                if isinstance(xy[i], complex):
                    ctrl[j, 0] = x0 + xy[i].real
                    ctrl[j, 1] = y0 + xy[i].imag
                    i += 1
                else:
                    ctrl[j, 0] = x0 + xy[i]
                    ctrl[j, 1] = y0 + xy[i + 1]
                    i += 2
            f = _func_bezier(ctrl)
            uu = [0, 0.2, 0.5, 0.8, 1]
            fu = [f(u) for u in uu]
            iu = 1
            while iu < len(fu):
                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                test_pt = f(test_u)
                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                if test_err[0] ** 2 + test_err[1] ** 2 > self.tol:
                    uu.insert(iu, test_u)
                    fu.insert(iu, test_pt)
                else:
                    iu += 1
            self.points.extend(xy for xy in fu[1:])
            self.last_c = ctrl[2]
        return self

    def Q(self, *xy):
        """
        Add quadratic Bezier curves to the curve.

        Parameters
        ----------
        xy : numbers
            Coordinate pairs. Each set of 2 pairs are interpreted as
            the control point and the endpoint of the curve.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = None
        i = 0
        while i < len(xy):
            ctrl = numpy.empty((3, 2))
            ctrl[0] = self.points[-1]
            for j in range(1, 3):
                if isinstance(xy[i], complex):
                    ctrl[j, 0] = xy[i].real
                    ctrl[j, 1] = xy[i].imag
                    i += 1
                else:
                    ctrl[j, 0] = xy[i]
                    ctrl[j, 1] = xy[i + 1]
                    i += 2
            f = _func_bezier(ctrl)
            uu = [0, 0.2, 0.5, 0.8, 1]
            fu = [f(u) for u in uu]
            iu = 1
            while iu < len(fu):
                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                test_pt = f(test_u)
                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                if test_err[0] ** 2 + test_err[1] ** 2 > self.tol:
                    uu.insert(iu, test_u)
                    fu.insert(iu, test_pt)
                else:
                    iu += 1
            self.points.extend(xy for xy in fu[1:])
        self.last_q = ctrl[1]
        return self

    def q(self, *xy):
        """
        Add quadratic Bezier curves to the curve.

        Parameters
        ----------
        xy : numbers
            Coordinate pairs. Each set of 2 pairs are interpreted as
            the control point and the endpoint of the curve.
            All coordinates are relative to the current end point.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = None
        x0, y0 = self.points[-1]
        i = 0
        while i < len(xy):
            ctrl = numpy.empty((3, 2))
            ctrl[0] = self.points[-1]
            for j in range(1, 3):
                if isinstance(xy[i], complex):
                    ctrl[j, 0] = x0 + xy[i].real
                    ctrl[j, 1] = y0 + xy[i].imag
                    i += 1
                else:
                    ctrl[j, 0] = x0 + xy[i]
                    ctrl[j, 1] = y0 + xy[i + 1]
                    i += 2
            f = _func_bezier(ctrl)
            uu = [0, 0.2, 0.5, 0.8, 1]
            fu = [f(u) for u in uu]
            iu = 1
            while iu < len(fu):
                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                test_pt = f(test_u)
                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                if test_err[0] ** 2 + test_err[1] ** 2 > self.tol:
                    uu.insert(iu, test_u)
                    fu.insert(iu, test_pt)
                else:
                    iu += 1
            self.points.extend(xy for xy in fu[1:])
        self.last_q = ctrl[1]
        return self

    def T(self, *xy):
        """
        Add smooth quadratic Bezier curves to the curve.

        Parameters
        ----------
        xy : numbers
            Coordinates of the endpoints of the curves.  The control
            point is assumed to be the reflection of the control point
            of the last curve relative to the starting point of the
            curve. If the previous curve is not a quadratic Bezier,
            the control point is coincident with the starting point.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = None
        if self.last_q is None:
            self.last_q = self.points[-1]
        i = 0
        while i < len(xy):
            ctrl = numpy.empty((3, 2))
            ctrl[0] = self.points[-1]
            ctrl[1] = 2 * ctrl[0] - self.last_q
            if isinstance(xy[i], complex):
                ctrl[2, 0] = xy[i].real
                ctrl[2, 1] = xy[i].imag
                i += 1
            else:
                ctrl[2, 0] = xy[i]
                ctrl[2, 1] = xy[i + 1]
                i += 2
            f = _func_bezier(ctrl)
            uu = [0, 0.2, 0.5, 0.8, 1]
            fu = [f(u) for u in uu]
            iu = 1
            while iu < len(fu):
                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                test_pt = f(test_u)
                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                if test_err[0] ** 2 + test_err[1] ** 2 > self.tol:
                    uu.insert(iu, test_u)
                    fu.insert(iu, test_pt)
                else:
                    iu += 1
            self.points.extend(xy for xy in fu[1:])
            self.last_q = ctrl[1]
        return self

    def t(self, *xy):
        """
        Add smooth quadratic Bezier curves to the curve.

        Parameters
        ----------
        xy : numbers
            Coordinates of the endpoints of the curves.  The control
            point is assumed to be the reflection of the control point
            of the last curve relative to the starting point of the
            curve. If the previous curve is not a quadratic Bezier,
            the control point is coincident with the starting point.
            All coordinates are relative to the current end point.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = None
        if self.last_q is None:
            self.last_q = self.points[-1]
        x0, y0 = self.points[-1]
        i = 0
        while i < len(xy):
            ctrl = numpy.empty((3, 2))
            ctrl[0] = self.points[-1]
            ctrl[1] = 2 * ctrl[0] - self.last_q
            if isinstance(xy[i], complex):
                ctrl[2, 0] = x0 + xy[i].real
                ctrl[2, 1] = y0 + xy[i].imag
                i += 1
            else:
                ctrl[2, 0] = x0 + xy[i]
                ctrl[2, 1] = y0 + xy[i + 1]
                i += 2
            f = _func_bezier(ctrl)
            uu = [0, 0.2, 0.5, 0.8, 1]
            fu = [f(u) for u in uu]
            iu = 1
            while iu < len(fu):
                test_u = 0.5 * (uu[iu - 1] + uu[iu])
                test_pt = f(test_u)
                test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
                if test_err[0] ** 2 + test_err[1] ** 2 > self.tol:
                    uu.insert(iu, test_u)
                    fu.insert(iu, test_pt)
                else:
                    iu += 1
            self.points.extend(xy for xy in fu[1:])
            self.last_q = ctrl[1]
        return self

    def B(self, *xy):
        """
        Add a general degree Bezier curve.

        Parameters
        ----------
        xy : numbers
            Coordinate pairs.  The last coordinate is the endpoint of
            curve and all other are control points.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = self.last_q = None
        i = 0
        ctrl = [self.points[-1]]
        while i < len(xy):
            if isinstance(xy[i], complex):
                ctrl.append((xy[i].real, xy[i].imag))
                i += 1
            else:
                ctrl.append((xy[i], xy[i + 1]))
                i += 2
        ctrl = numpy.array(ctrl)
        f = _func_bezier(ctrl)
        uu = numpy.linspace(-1, 1, ctrl.shape[0] + 1)
        uu = list(0.5 * (1 + numpy.sign(uu) * numpy.abs(uu) ** 0.8))
        fu = [f(u) for u in uu]
        iu = 1
        while iu < len(fu):
            test_u = 0.5 * (uu[iu - 1] + uu[iu])
            test_pt = f(test_u)
            test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
            if test_err[0] ** 2 + test_err[1] ** 2 > self.tol:
                uu.insert(iu, test_u)
                fu.insert(iu, test_pt)
            else:
                iu += 1
        self.points.extend(xy for xy in fu[1:])
        return self

    def b(self, *xy):
        """
        Add a general degree Bezier curve.

        Parameters
        ----------
        xy : numbers
            Coordinate pairs.  The last coordinate is the endpoint of
            curve and all other are control points.  All coordinates are
            relative to the current end point.

        Returns
        -------
        out : `Curve`
            This curve.
        """
        self.last_c = self.last_q = None
        x0, y0 = self.points[-1]
        i = 0
        ctrl = [self.points[-1]]
        while i < len(xy):
            if isinstance(xy[i], complex):
                ctrl.append((x0 + xy[i].real, y0 + xy[i].imag))
                i += 1
            else:
                ctrl.append((x0 + xy[i], y0 + xy[i + 1]))
                i += 2
        ctrl = numpy.array(ctrl)
        f = _func_bezier(ctrl)
        uu = numpy.linspace(-1, 1, ctrl.shape[0] + 1)
        uu = list(0.5 * (1 + numpy.sign(uu) * numpy.abs(uu) ** 0.8))
        fu = [f(u) for u in uu]
        iu = 1
        while iu < len(fu):
            test_u = 0.5 * (uu[iu - 1] + uu[iu])
            test_pt = f(test_u)
            test_err = 0.5 * (fu[iu - 1] + fu[iu]) - test_pt
            if test_err[0] ** 2 + test_err[1] ** 2 > self.tol:
                uu.insert(iu, test_u)
                fu.insert(iu, test_pt)
            else:
                iu += 1
        self.points.extend(xy for xy in fu[1:])
        return self

    def I(
        self,
        points,
        angles=None,
        curl_start=1,
        curl_end=1,
        t_in=1,
        t_out=1,
        cycle=False,
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

        Returns
        -------
        out : `Curve`
            This curve.

        Examples
        --------
        >>> c1 = gdspy.Curve(0, 1).I([(1, 1), (2, 1), (1, 0)])
        >>> c2 = gdspy.Curve(0, 2).I([(1, 2), (2, 2), (1, 1)],
        ...                          cycle=True)
        >>> ps = gdspy.PolygonSet([c1.get_points(), c2.get_points()])

        References
        ----------
        .. [1] Hobby, J.D.  *Discrete Comput. Geom.* (1986) 1: 123.
           `DOI: 10.1007/BF02187690
           <https://doi.org/10.1007/BF02187690>`_
        """
        pts = numpy.vstack((self.points[-1:], points))
        cta, ctb = _hobby(pts, angles, curl_start, curl_end, t_in, t_out, cycle)
        args = []
        args.extend(
            x
            for i in range(pts.shape[0] - 1)
            for x in [
                cta[i, 0],
                cta[i, 1],
                ctb[i, 0],
                ctb[i, 1],
                pts[i + 1, 0],
                pts[i + 1, 1],
            ]
        )
        if cycle:
            args.extend(
                [cta[-1, 0], cta[-1, 1], ctb[-1, 0], ctb[-1, 1], pts[0, 0], pts[0, 1]]
            )
        return self.C(*args)

    def i(
        self,
        points,
        angles=None,
        curl_start=1,
        curl_end=1,
        t_in=1,
        t_out=1,
        cycle=False,
    ):
        """
        Add a smooth interpolating curve through the given points.

        Uses the Hobby algorithm [1]_ to calculate a smooth
        interpolating curve made of cubic Bezier segments between each
        pair of points.

        Parameters
        ----------
        points : array-like[N][2]
            Vertices in the interpolating curve (relative to the current
            endpoint).
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

        Returns
        -------
        out : `Curve`
            This curve.

        Examples
        --------
        >>> c1 = gdspy.Curve(0, 1).i([(1, 0), (2, 0), (1, -1)])
        >>> c2 = gdspy.Curve(0, 2).i([(1, 0), (2, 0), (1, -1)],
        ...                          cycle=True)
        >>> ps = gdspy.PolygonSet([c1.get_points(), c2.get_points()])

        References
        ----------
        .. [1] Hobby, J.D.  *Discrete Comput. Geom.* (1986) 1: 123.
           `DOI: 10.1007/BF02187690
           <https://doi.org/10.1007/BF02187690>`_
        """
        pts = numpy.vstack((numpy.array(((0.0, 0.0),)), points)) + self.points[-1]
        cta, ctb = _hobby(pts, angles, curl_start, curl_end, t_in, t_out, cycle)
        args = []
        args.extend(
            x
            for i in range(pts.shape[0] - 1)
            for x in [
                cta[i, 0],
                cta[i, 1],
                ctb[i, 0],
                ctb[i, 1],
                pts[i + 1, 0],
                pts[i + 1, 1],
            ]
        )
        if cycle:
            args.extend(
                [cta[-1, 0], cta[-1, 1], ctb[-1, 0], ctb[-1, 1], pts[0, 0], pts[0, 1]]
            )
        return self.C(*args)
