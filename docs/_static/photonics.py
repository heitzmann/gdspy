######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import numpy
import gdspy


def grating(period, number_of_teeth, fill_frac, width, position, direction, lda=1, sin_theta=0,
            focus_distance=-1, focus_width=-1, tolerance=0.001, layer=0, datatype=0):
    '''
    Straight or focusing grating.

    period          : grating period
    number_of_teeth : number of teeth in the grating
    fill_frac       : filling fraction of the teeth (w.r.t. the period)
    width           : width of the grating
    position        : grating position (feed point)
    direction       : one of {'+x', '-x', '+y', '-y'}
    lda             : free-space wavelength
    sin_theta       : sine of incidence angle
    focus_distance  : focus distance (negative for straight grating)
    focus_width     : if non-negative, the focusing area is included in
                      the result (usually for negative resists) and this
                      is the width of the waveguide connecting to the
                      grating
    tolerance       : same as in `path.parametric`
    layer           : GDSII layer number
    datatype        : GDSII datatype number

    Return `PolygonSet`
    '''
    if focus_distance < 0:
        path = gdspy.L1Path((position[0] - 0.5 * width, position[1] + 0.5 * (number_of_teeth - 1 + fill_frac) * period),
                            '+x', period * fill_frac, [width], [], number_of_teeth, period, layer=layer, datatype=datatype)
    else:
        neff = lda / float(period) + sin_theta
        qmin = int(focus_distance / float(period) + 0.5)
        path = gdspy.Path(period * fill_frac, position)
        c3 = neff**2 - sin_theta**2
        w = 0.5 * width
        for q in range(qmin, qmin + number_of_teeth):
            c1 = q * lda * sin_theta
            c2 = (q * lda)**2
            path.parametric(lambda t: (width * t - w, (c1 + neff * numpy.sqrt(c2 - c3 * (width * t - w)**2)) / c3),
                            tolerance=tolerance, max_points=0, layer=layer, datatype=datatype)
            path.x = position[0]
            path.y = position[1]
        sz = path.polygons[0].shape[0] // 2
        if focus_width == 0:
            path.polygons[0] = numpy.vstack((path.polygons[0][:sz, :], [position]))
        elif focus_width > 0:
            path.polygons[0] = numpy.vstack((path.polygons[0][:sz, :],
                                             [(position[0] + 0.5 * focus_width, position[1]),
                                              (position[0] - 0.5 * focus_width, position[1])]))
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
    w = 0.45

    # Negative resist example
    ring = gdspy.Cell('NRing')
    ring.add(gdspy.Round((20, 0), 20, 20 - w, tolerance=0.001))

    grat = gdspy.Cell('NGrat')
    grat.add(grating(0.626, 28, 0.5, 19, (0, 0), '+y', 1.55,
                     numpy.sin(numpy.pi * 8 / 180), 21.5, w,
                     tolerance=0.001))

    taper = gdspy.Cell('NTaper')
    taper.add(gdspy.Path(0.12, (0, 0)).segment(50, '+y', final_width=w))

    c = gdspy.Cell('Negative')
    for i in range(8):
        path = gdspy.SimplePath([(150 * i, 50)], width=w, corners='circular bend', bend_radius=50, gdsii_path=True)
        path.segment((0, 600 - 20 * i), relative=True)
        path.segment((500, 0), relative=True)
        path.segment((0, 300 + 20 * i), relative=True)
        c.add(path)
        c.add(gdspy.CellReference(ring, (150 * i + w / 2 + 0.06 + 0.02 * i, 300)))
    c.add(gdspy.CellArray(taper, 8, 1, (150, 0), (0, 0)))
    c.add(gdspy.CellArray(grat, 8, 1, (150, 0), (500, 950)))

    # Positive resist example
    ring_edge = gdspy.Rectangle((0, -50), (70, 50))
    ring_hole = gdspy.Round((20, 0), 20, 20 - w, tolerance=0.001)
    ring_path = gdspy.Path(5, (0, 50), number_of_paths=2, distance=5 + w).segment(400, '+y')

    grat = gdspy.Cell('PGrat')
    grat.add(grating(0.626, 28, 0.5, 19, (0, 0), '+y', 1.55,
                     numpy.sin(numpy.pi * 8 / 180), 21.5, tolerance=0.001))
    grat.add(gdspy.Path(5, (0, 0), number_of_paths=2, distance=5 + w).segment(
        21.5, '+y', final_distance=5 + 19))

    taper = gdspy.Cell('PTaper')
    taper.add(gdspy.Path(20, (0, 0), number_of_paths=2, distance=20 + 0.12).segment(
        50, '+y', final_width=5, final_distance=5 + w))

    c = gdspy.Cell('Positive')

    for i in range(8):
        path = gdspy.SimplePath([(150 * i, 450)], width=[5, 5], offset=5 + w, gdsii_path=True)
        path.segment((0, 150 - 20 * i), relative=True)
        path.turn(50, 'r')
        path.segment((400, 0), relative=True)
        path.turn(50, 'l')
        path.segment((0, 250 + 20 * i), relative=True)
        c.add(path)
        dx = w / 2 + 0.06 + 0.2 * i
        c.add(gdspy.boolean(
            gdspy.boolean(ring_path, gdspy.copy(ring_edge, dx, 300), 'or', precision=1e-4),
            gdspy.copy(ring_hole, dx, 300), 'not', precision=1e-4).translate(150 * i, 0))
    c.add(gdspy.CellArray(grat, 8, 1, (150, 0), (500, 950)))
    c.add(gdspy.CellArray(taper, 8, 1, (150, 0), (0, 0)))

    gdspy.write_gds('photonics.gds')
    gdspy.LayoutViewer()
