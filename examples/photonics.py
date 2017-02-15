######################################################################
#                                                                    #
#  Copyright 2009-2017 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import numpy
import gdspy


def waveguide(path,
              points,
              finish,
              bend_radius,
              number_of_points=0.01,
              direction=None,
              layer=0,
              datatype=0):
    '''
    Easy waveguide creation tool with absolute positioning.

    path             : `gdspy.Path` which servers as base for the waveguide
    points           : coordinates along which the waveguide will travel
    finish           : end point of the waveguide
    bend_radius      : radius of the turns in the waveguide
    number_of_points : same as in `path.turn`
    direction        : starting direction
    layer            : GDSII layer number
    datatype         : GDSII datatype number

    Return `path`.
    '''
    if direction is not None:
        path.direction = direction
    axis = 0 if path.direction[1] == 'x' else 1
    points.append(finish[(axis + len(points)) % 2])
    n = len(points)
    if points[0] > (path.x, path.y)[axis]:
        path.direction = ['+x', '+y'][axis]
    else:
        path.direction = ['-x', '-y'][axis]
    for i in range(n):
        path.segment(
            abs(points[i] - (path.x, path.y)[axis]) - bend_radius,
            layer=layer,
            datatype=datatype)
        axis = 1 - axis
        if i < n - 1:
            goto = points[i + 1]
        else:
            goto = finish[axis]
        if (goto > (path.x, path.y)[axis]) ^ ((path.direction[0] == '+') ^
                                              (path.direction[1] == 'x')):
            bend = 'l'
        else:
            bend = 'r'
        path.turn(
            bend_radius,
            bend,
            number_of_points=number_of_points,
            layer=layer,
            datatype=datatype)
    return path.segment(
        abs(finish[axis] - (path.x, path.y)[axis]),
        layer=layer,
        datatype=datatype)


def taper(path,
          length,
          final_width,
          final_distance,
          direction=None,
          layer=0,
          datatype=0):
    '''
    Linear tapers for the lazy.

    path        : `gdspy.Path` to append the taper
    length      : total length
    final_width : final width of th taper
    direction   : taper direction
    layer       : GDSII layer number (int or list)
    datatype    : GDSII datatype number (int or list)

    Parameters `layer` and `datatype` must be of the same type. If they are
    lists, they must have the same length. Their length indicate the number of
    pieces that compose the taper.

    Return `path`.
    '''
    if layer.__class__ == datatype.__class__ == [].__class__:
        assert len(layer) == len(datatype)
    elif isinstance(layer, int) and isinstance(datatype, int):
        layer = [layer]
        datatype = [datatype]
    else:
        raise ValueError('Parameters layer and datatype must have the same '
                         'type (either int or list) and length.')
    n = len(layer)
    w = numpy.linspace(2 * path.w, final_width, n + 1)[1:]
    d = numpy.linspace(path.distance, final_distance, n + 1)[1:]
    l = float(length) / n
    for i in range(n):
        path.segment(
            l, direction, w[i], d[i], layer=layer[i], datatype=datatype[i])
    return path


def grating(period,
            number_of_teeth,
            fill_frac,
            width,
            position,
            direction,
            lda=1,
            sin_theta=0,
            focus_distance=-1,
            focus_width=-1,
            evaluations=99,
            layer=0,
            datatype=0):
    '''
    Straight or focusing grating.

    period          : grating period
    number_of_teeth : number of teeth in the grating
    fill_frac       : filling fraction of the teeth with respect to the period
    width           : width of the grating
    position        : grating position (feed point)
    direction       : one of {'+x', '-x', '+y', '-y'}
    lda             : free-space wavelength
    sin_theta       : sine of incidence angle
    focus_distance  : focus distance (negative for straight grating)
    focus_width     : if non-negative, the focusing area is included in the
                      result (usually for negative resists) and this is the
                      width of the waveguide connecting to the grating
    evaluations     : number of parametric evaluations of `path.parametric`
    layer           : GDSII layer number
    datatype        : GDSII datatype number

    Return `PolygonSet`
    '''
    if focus_distance < 0:
        path = gdspy.L1Path(
            (position[0] - 0.5 * width,
             position[1] + 0.5 * (number_of_teeth - 1 + fill_frac) * period),
            '+x',
            period * fill_frac, [width], [],
            number_of_teeth,
            period,
            layer=layer,
            datatype=datatype)
    else:
        neff = lda / float(period) + sin_theta
        qmin = int(focus_distance / float(period) + 0.5)
        path = gdspy.Path(period * fill_frac, position)
        max_points = 199 if focus_width < 0 else 2 * evaluations
        c3 = neff**2 - sin_theta**2
        w = 0.5 * width
        for q in range(qmin, qmin + number_of_teeth):
            c1 = q * lda * sin_theta
            c2 = (q * lda)**2
            path.parametric(
                lambda t: (width * t - w, (c1 + neff * numpy.sqrt(c2 - c3 * (
                    width * t - w)**2)) / c3),
                number_of_evaluations=evaluations,
                max_points=max_points,
                layer=layer,
                datatype=datatype)
            path.x = position[0]
            path.y = position[1]
        if focus_width >= 0:
            path.polygons[0] = numpy.vstack(
                (path.polygons[0][:evaluations, :],
                 ([position] if focus_width == 0 else
                  [(position[0] + 0.5 * focus_width, position[1]),
                   (position[0] - 0.5 * focus_width, position[1])])))
            path.fracture()
    if direction == '-x':
        return path.rotate(0.5 * numpy.pi, position)
    elif direction == '+x':
        return path.rotate(-0.5 * numpy.pi, position)
    elif direction == '-y':
        return path.rotate(numpy.pi, position)
    else:
        return path


if __name__ == '__main__':
    # Examples

    c1 = gdspy.Cell('Example1')
    # Waveguide starts at (0, 0)...
    path = gdspy.Path(0.450, (0, 0))
    # ... then starts in the +y direction up to y=200, then through x=500,
    # and stops at (800, 400). All bends have radius 50.
    waveguide(path, [200, 500], (800, 400), 50, direction='+y')
    c1.add(path)

    # More useful example including a taper and a grating:
    c2 = gdspy.Cell('Example2')
    for i in range(3):
        path = gdspy.Path(0.120, (50 * i, 0))
        taper(
            path,
            75,
            0.450,
            0,
            '+y',
            layer=list(range(5, 1, -1)),
            datatype=list(range(1, 5)))
        waveguide(
            path, [300 - 20 * i], (500 + 50 * i, 425), 50, layer=1, datatype=1)
        c2.add(path)
        c2.add(
            grating(
                0.626,
                28,
                0.5,
                19, (path.x, path.y),
                path.direction,
                1.55,
                numpy.sin(numpy.pi * 8 / 180),
                21.5,
                2 * path.w,
                layer=1,
                datatype=1))

    # Straight grating and positive resist example
    c3 = gdspy.Cell('Example3')
    spec = {'layer': 4, 'datatype': 3}
    lda = 1.55
    gr_width = 10
    gr_per = 0.626
    gr_teeth = 20
    wg_clad = 4
    wg_width = 0.45
    tp_len = 700
    c3.add(
        grating(gr_per, gr_teeth, 0.5, gr_width, (gr_per * gr_teeth, 0), '-x',
                lda, numpy.sin(numpy.pi * 8 / 180), **spec))
    path = gdspy.Path(
        wg_clad, (0, 0), number_of_paths=2, distance=gr_width + wg_clad)
    path.segment(gr_per * gr_teeth, '+x', **spec)
    taper(path, tp_len, wg_clad, wg_width + wg_clad, **spec)
    waveguide(path, [800], (tp_len + gr_per * gr_teeth, 200), 50, 0.1, **spec)
    taper(path, tp_len, wg_clad, gr_width + wg_clad, **spec)
    c3.add(
        grating(gr_per, gr_teeth, 0.5, gr_width, (path.x, path.y),
                path.direction, lda, numpy.sin(numpy.pi * 8 / 180), **spec))
    path.segment(gr_per * gr_teeth, **spec)
    c3.add(path)

    gdspy.LayoutViewer()
