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
        z = numpy.hstack((numpy.roll(z, -rotate), z[rotate : rotate + 1]))
        t_in = numpy.hstack((numpy.roll(t_in, -rotate), t_in[rotate : rotate + 1]))
        t_out = numpy.hstack((numpy.roll(t_out, -rotate), t_out[rotate : rotate + 1]))
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
        m[ni, i] = d[i1] * t_in[i2] * t_in[i1] ** 2
        # B_{i+1}
        m[ni, i1] = -d[i] * t_out[i] * t_out[i1] ** 2 * (1 - 3 * t_in[i2])
        # C_{i+1}
        m[ni, ni] = d[i1] * t_in[i2] * t_in[i1] ** 2 * (1 - 3 * t_out[i])
        # D_{i+2}
        m[ni, n + i1] = -d[i] * t_out[i] * t_out[i1] ** 2
        sol = numpy.linalg.solve(m, coef)
        theta = sol[:n]
        phi = sol[n:]
        w = numpy.exp(1j * (theta + delta))
        a = 2 ** 0.5
        b = 1.0 / 16
        c = (3 - 5 ** 0.5) / 2
        sintheta = numpy.sin(theta)
        costheta = numpy.cos(theta)
        sinphi = numpy.sin(phi)
        cosphi = numpy.cos(phi)
        alpha = (
            a * (sintheta - b * sinphi) * (sinphi - b * sintheta) * (costheta - cosphi)
        )
        cta = z + w * d * ((2 + alpha) / (1 + (1 - c) * costheta + c * cosphi)) / (
            3 * t_out
        )
        ctb = numpy.roll(z, -1) - numpy.roll(w, -1) * d * (
            (2 - alpha) / (1 + (1 - c) * cosphi + c * costheta)
        ) / (3 * numpy.roll(t_in, -1))
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
            coef[1:nn] = -psi[i : j - 1]
            m = numpy.zeros((2 * nn, 2 * nn))
            if nn > 1:
                ii = numpy.arange(nn - 1)  # [0 .. nn-2]
                i0 = i + ii  # [i .. j-1]
                i1 = 1 + i0  # [i+1 .. j]
                i2 = 2 + i0  # [i+2 .. j+1]
                ni = nn + ii  # [nn .. 2*nn-2]
                ii1 = 1 + ii  # [1 .. nn-1]
                m[ii1, ii1] = 1
                m[ii1, ni] = 1
                # A_ii
                m[ni, ii] = d[i1] * t_in[i2] * t_in[i1] ** 2
                # B_{ii+1}
                m[ni, ii1] = -d[i0] * t_out[i0] * t_out[i1] ** 2 * (1 - 3 * t_in[i2])
                # C_{ii+1}
                m[ni, ni] = d[i1] * t_in[i2] * t_in[i1] ** 2 * (1 - 3 * t_out[i0])
                # D_{ii+2}
                m[ni, ni + 1] = -d[i0] * t_out[i0] * t_out[i1] ** 2
            if angles[i] is None:
                to3 = t_out[0] ** 3
                cti3 = curl_start * t_in[1] ** 3
                # B_0
                m[0, 0] = to3 * (1 - 3 * t_in[1]) - cti3
                # D_1
                m[0, nn] = to3 - cti3 * (1 - 3 * t_out[0])
            else:
                coef[0] = theta[i]
                m[0, 0] = 1
                m[0, nn] = 0
            if angles[j] is None:
                ti3 = t_in[n] ** 3
                cto3 = curl_end * t_out[n - 1] ** 3
                # A_{nn-1}
                m[2 * nn - 1, nn - 1] = ti3 - cto3 * (1 - 3 * t_in[n])
                # C_nn
                m[2 * nn - 1, 2 * nn - 1] = ti3 * (1 - 3 * t_out[n - 1]) - cto3
            else:
                coef[2 * nn - 1] = phi[j - 1]
                m[2 * nn - 1, nn - 1] = 0
                m[2 * nn - 1, 2 * nn - 1] = 1
            if nn > 1 or angles[i] is None or angles[j] is None:
                # print("range:", i, j)
                # print("A =", m)
                # print("b =", coef)
                sol = numpy.linalg.solve(m, coef)
                # print("x =", sol)
                theta[i:j] = sol[:nn]
                phi[i:j] = sol[nn:]
            i = j
        w = numpy.hstack(
            (numpy.exp(1j * (delta + theta)), numpy.exp(1j * (delta[-1:] - phi[-1:])))
        )
        a = 2 ** 0.5
        b = 1.0 / 16
        c = (3 - 5 ** 0.5) / 2
        sintheta = numpy.sin(theta)
        costheta = numpy.cos(theta)
        sinphi = numpy.sin(phi)
        cosphi = numpy.cos(phi)
        alpha = (
            a * (sintheta - b * sinphi) * (sinphi - b * sintheta) * (costheta - cosphi)
        )
        cta = z[:-1] + w[:-1] * d * (
            (2 + alpha) / (1 + (1 - c) * costheta + c * cosphi)
        ) / (3 * t_out[:-1])
        ctb = z[1:] - w[1:] * d * (
            (2 - alpha) / (1 + (1 - c) * cosphi + c * costheta)
        ) / (3 * t_in[1:])
        if rotate > 0:
            cta = numpy.roll(cta, rotate)
            ctb = numpy.roll(ctb, rotate)
    return (
        numpy.vstack((cta.real, cta.imag)).transpose(),
        numpy.vstack((ctb.real, ctb.imag)).transpose(),
    )



