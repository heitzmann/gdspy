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
        p = gdspy.L1Path((position[0] - 0.5 * width,
                          position[1] + 0.5 * (number_of_teeth - 1 + fill_frac) * period),
                         '+x', period * fill_frac, [width], [], number_of_teeth, period,
                         layer=layer, datatype=datatype)
    else:
        neff = lda / float(period) + sin_theta
        qmin = int(focus_distance / float(period) + 0.5)
        p = gdspy.Path(period * fill_frac, position)
        c3 = neff**2 - sin_theta**2
        w = 0.5 * width
        for q in range(qmin, qmin + number_of_teeth):
            c1 = q * lda * sin_theta
            c2 = (q * lda)**2
            p.parametric(lambda t: (width * t - w,
                                    (c1 + neff * numpy.sqrt(c2 - c3 * (width * t - w)**2)) / c3),
                         tolerance=tolerance, max_points=0, layer=layer, datatype=datatype)
            p.x = position[0]
            p.y = position[1]
        sz = p.polygons[0].shape[0] // 2
        if focus_width == 0:
            p.polygons[0] = numpy.vstack((p.polygons[0][:sz, :], [position]))
        elif focus_width > 0:
            p.polygons[0] = numpy.vstack((p.polygons[0][:sz, :],
                                          [(position[0] + 0.5 * focus_width, position[1]),
                                           (position[0] - 0.5 * focus_width, position[1])]))
        p.fracture()
    if direction == '-x':
        return p.rotate(0.5 * numpy.pi, position)
    elif direction == '+x':
        return p.rotate(-0.5 * numpy.pi, position)
    elif direction == '-y':
        return p.rotate(numpy.pi, position)
    else:
        return p


if __name__ == '__main__':
    # Examples

    # Negative resist example
    width = 0.45
    bend_radius = 50.0
    ring_radius = 20.0
    taper_len = 50.0
    input_gap = 150.0
    io_gap = 500.0
    wg_gap = 20.0
    ring_gaps = [0.06 + 0.02 * i for i in range(8)]

    ring = gdspy.Cell('NRing')
    ring.add(gdspy.Round((ring_radius, 0), ring_radius, ring_radius - width, tolerance=0.001))

    grat = gdspy.Cell('NGrat')
    grat.add(grating(0.626, 28, 0.5, 19, (0, 0), '+y', 1.55,
                     numpy.sin(numpy.pi * 8 / 180), 21.5, width,
                     tolerance=0.001))

    taper = gdspy.Cell('NTaper')
    taper.add(gdspy.Path(0.12, (0, 0)).segment(taper_len, '+y', final_width=width))

    c = gdspy.Cell('Negative')
    for i, gap in enumerate(ring_gaps):
        path = gdspy.FlexPath([(input_gap * i, taper_len)], width=width,
                              corners='circular bend', bend_radius=bend_radius,
                              gdsii_path=True)
        path.segment((0, 600 - wg_gap * i), relative=True)
        path.segment((io_gap, 0), relative=True)
        path.segment((0, 300 + wg_gap * i), relative=True)
        c.add(path)
        c.add(gdspy.CellReference(ring, (input_gap * i + width / 2 + gap, 300)))
    c.add(gdspy.CellArray(taper, len(ring_gaps), 1, (input_gap, 0), (0, 0)))
    c.add(gdspy.CellArray(grat, len(ring_gaps), 1, (input_gap, 0), (io_gap, 900 + taper_len)))


    # Positive resist example
    width = 0.45
    ring_radius = 20.0
    big_margin = 10.0
    small_margin = 5.0
    taper_len = 50.0
    bus_len = 400.0
    input_gap = 150.0
    io_gap = 500.0
    wg_gap = 20.0
    ring_gaps = [0.06 + 0.02 * i for i in range(8)]

    ring_margin = gdspy.Rectangle((0, -ring_radius - big_margin),
                                  (2 * ring_radius + big_margin, ring_radius + big_margin))
    ring_hole = gdspy.Round((ring_radius, 0), ring_radius, ring_radius - width, tolerance=0.001)
    ring_bus = gdspy.Path(small_margin, (0, taper_len), number_of_paths=2,
                          distance=small_margin + width)
    ring_bus.segment(bus_len, '+y')

    p = gdspy.Path(small_margin, (0, 0), number_of_paths=2, distance=small_margin + width)
    p.segment(21.5, '+y', final_distance=small_margin + 19)
    grat = gdspy.Cell('PGrat').add(p)
    grat.add(grating(0.626, 28, 0.5, 19, (0, 0), '+y', 1.55,
                     numpy.sin(numpy.pi * 8 / 180), 21.5, tolerance=0.001))

    p = gdspy.Path(big_margin, (0, 0), number_of_paths=2, distance=big_margin + 0.12)
    p.segment(taper_len, '+y', final_width=small_margin, final_distance=small_margin + width)
    taper = gdspy.Cell('PTaper').add(p)

    c = gdspy.Cell('Positive')
    for i, gap in enumerate(ring_gaps):
        path = gdspy.FlexPath([(input_gap * i, taper_len + bus_len)],
                              width=[small_margin, small_margin],
                              offset=small_margin + width, gdsii_path=True)
        path.segment((0, 600 - bus_len - bend_radius - wg_gap * i), relative=True)
        path.turn(bend_radius, 'r')
        path.segment((io_gap - 2 * bend_radius, 0), relative=True)
        path.turn(bend_radius, 'l')
        path.segment((0, 300 - bend_radius + wg_gap * i), relative=True)
        c.add(path)
        dx = width / 2 + gap
        c.add(gdspy.boolean(
            gdspy.boolean(ring_bus, gdspy.copy(ring_margin, dx, 300), 'or', precision=1e-4),
            gdspy.copy(ring_hole, dx, 300), 'not', precision=1e-4).translate(input_gap * i, 0))
    c.add(gdspy.CellArray(taper, len(ring_gaps), 1, (input_gap, 0), (0, 0)))
    c.add(gdspy.CellArray(grat, len(ring_gaps), 1, (input_gap, 0), (io_gap, 900 + taper_len)))


    # Save to a gds file and check out the output
    gdspy.write_gds('photonics.gds')
    gdspy.LayoutViewer()
