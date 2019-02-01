import gdspy
import numpy
import colorsys
from PIL import Image, ImageDraw


class ColorDict(dict):
    def __missing__(self, key):
        layer, datatype = key
        rgb = tuple(int(255 * c + 0.5) for c in colorsys.hsv_to_rgb((layer % 3) / 3.0 + (layer % 6 // 3) / 6.0 + (layer // 6) / 11.0, 1 - ( (layer + datatype) % 8) / 12.0, 1 - (datatype % 3) / 4.0))
        self[key] = rgb
        return rgb

color = ColorDict()

def draw(cell, name=None, width=600, height=400, margin=20):
    global color
    (ax, ay), (bx, by) = cell.get_bounding_box()
    ax = min(0, ax)
    ay = min(0, ay)
    bx = max(1, bx)
    by = max(1, by)
    sx = 3 * width - 2 * margin
    sy = 3 * height - 2 * margin
    bx -= ax
    by -= ay
    if bx * sy > by * sx:
        scale = sx / bx
        sy = int(sx * by / bx + 0.5)
    else:
        scale = sy / by
        sx = int(sy * bx / by + 0.5)
    ox = margin
    oy = sy + margin
    width = int((sx + 2 * margin) / 3 + 0.5)
    height = int((sy + 2 * margin) / 3 + 0.5)
    img = Image.new('RGBA', (3 * width, 3 * height), (0, 0, 0, 0))
    for key, polys in cell.get_polygons(by_spec=True).items():
        lc = color[key]
        fc = color[key] + (128,)
        for p in polys:
            p[:, 0] = ox + scale * (p[:, 0] - ax)
            p[:, 1] = oy - scale * (p[:, 1] - ay)
            pts = list(p.flatten())
            tmp = Image.new('RGBA', img.size, (0, 0, 0, 0))
            dr = ImageDraw.Draw(tmp)
            dr.polygon(pts, fill=fc, outline=lc)
            img = Image.alpha_composite(img, tmp)
    z = ox - scale * ax, oy + scale * ay
    p = ox + scale * (1 - ax), oy - scale * (1 - ay)
    #n = ox + scale * (-1 - ax), oy - scale * (-1 - ay)
    dr = ImageDraw.Draw(img)
    dr.line([z[0], z[1], p[0], z[1]], fill=(0, 0, 0, 255), width=3)
    dr.line([z[0], z[1], z[0], p[1]], fill=(0, 0, 0, 255), width=3)
    img = img.resize((width, height), Image.ANTIALIAS)
    img.save('_static/' + (cell.name if name is None else name) + '.png')


def main():
    # Polygons
    # Create a polygon from a list of vertices
    points = [(0, 0), (2, 2), (2, 6), (-6, 6),
              (-6, -6), (-4, -4), (-4, 4), (0, 4)]
    poly = gdspy.Polygon(points)
    draw(gdspy.Cell('polygons').add(poly))

    # Holes
    # Manually connect the hole to the outer boundary
    cutout = gdspy.Polygon([(0, 0), (5, 0), (5, 5), (0, 5), (0, 0),
                            (2, 2), (2, 3), (3, 3), (3, 2), (2, 2)])
    draw(gdspy.Cell('holes').add(cutout))

    # Circles
    # Circle centered at (0, 0), with radius 2 and tolerance 0.1
    circle = gdspy.Round((0, 0), 2, tolerance=0.01)

    # To create an ellipse, simply pass a list with 2 radii.
    # Because the tolerance is small (resulting a large number of
    # vertices), the ellipse is fractured in 2 polygons.
    ellipse = gdspy.Round((4, 0), [1, 2], tolerance=1e-4)

    # Circular arc example
    arc = gdspy.Round((2, 4), 2, inner_radius=1,
                      initial_angle=-0.2 * numpy.pi,
                      final_angle=1.2 * numpy.pi,
                      tolerance=0.01)
    draw(gdspy.Cell('circles').add([circle, ellipse, arc]))

    # Transformations
    poly = gdspy.Rectangle((-2, -2), (2, 2))
    poly.rotate(numpy.pi / 4)
    poly.scale(1, 0.5)
    draw(gdspy.Cell('transformations').add(poly))

    # Layer and Datatype
    # Layer/datatype definitions for each step in the fabrication
    ld_fulletch = {'layer': 1, 'datatype': 3}
    ld_partetch = {'layer': 2, 'datatype': 3}
    ld_liftoff = {'layer': 0, 'datatype': 7}

    p1 = gdspy.Rectangle((-3, -3), (3, 3), **ld_fulletch)
    p2 = gdspy.Rectangle((-5, -3), (-3, 3), **ld_partetch)
    p3 = gdspy.Rectangle((5, -3), (3, 3), **ld_partetch)
    p4 = gdspy.Round((0, 0), 2.5, number_of_points=6, **ld_liftoff)
    draw(gdspy.Cell('layer_and_datatype').add([p1, p2, p3, p4]))

    # References
    # Create a cell with a component that is used repeatedly
    contact = gdspy.Cell('CONTACT')
    contact.add([p1, p2, p3, p4])

    # Create a cell with the complete device
    device = gdspy.Cell('DEVICE')
    device.add(cutout)
    # Add 2 references to the component changing size and orientation
    ref1 = gdspy.CellReference(contact, (3.5, 1), magnification=0.25)
    ref2 = gdspy.CellReference(contact, (1, 3.5), magnification=0.25,
                               rotation=90)
    device.add([ref1, ref2])

    # The final layout has several repetitions of the complete device
    main = gdspy.Cell('MAIN')
    main.add(gdspy.CellArray(device, 3, 2, (6, 7)))
    draw(main, 'references')

    # Polygonal-Only Paths
    # Start a path at (0, 0) with width 1
    path1 = gdspy.Path(1, (0, 0))

    # Add a segment to the path goin in the '+y' direction
    path1.segment(4, '+y')

    # Further segments or turns will folow the current path direction
    # to ensure continuity
    path1.turn(2, 'r')
    path1.segment(1)
    path1.turn(3, 'rr')
    draw(gdspy.Cell('polygonal-only_paths').add(path1))

    # Polygonal-Only Paths 1
    path2 = gdspy.Path(0.5, (0, 0))

    # Start the path with a smooth Bezier S-curve
    path2.bezier([(0, 5), (5, 5), (5, 10)])

    # We want to add a spiral curve to the path.  The spiral is defined
    # as a parametric curve.  We make sure spiral(0) = (0, 0) so that
    # the path is continuous.
    def spiral(u):
        r = 4 - 3 * u
        theta = 5 * u * numpy.pi
        x = r * numpy.cos(theta) - 4
        y = r * numpy.sin(theta)
        return (x, y)

    # It is recommended to also define the derivative of the parametric
    # curve, otherwise this derivative must be calculated nummerically.
    # The derivative is used to define the side boundaries of the path,
    # so, in this case, to ensure continuity with the existing S-curve,
    # we make sure the the direction at the start of the spiral is
    # pointing exactly upwards, as if is radius were constant.
    # Additionally, the exact magnitude of the derivative is not
    # important; gdspy only uses its direction.
    def dspiral_dt(u):
        theta = 5 * u * numpy.pi
        dx_dt = -numpy.sin(theta)
        dy_dt = numpy.cos(theta)
        return (dx_dt, dy_dt)

    # Add the parametric spiral to the path
    path2.parametric(spiral, dspiral_dt)
    draw(gdspy.Cell('polygonal-only_paths_1').add(path2))

    # Polygonal-Only Paths 2
    # Start 3 parallel paths with center-to-center distance of 1.5
    path3 = gdspy.Path(0.1, (-5.5, 3), number_of_paths=3, distance=1.5)

    # Add a segment tapering the widths up to 0.5
    path3.segment(2, '-y', final_width=0.5)

    # Add a bezier curve decreasing the distance between paths to 0.75
    path3.bezier([(0, -2), (1, -3), (3, -3)], final_distance=0.75)

    # Add a parametric section to modulate the width with a sinusoidal
    # shape.  Note that the algorithm that determines the number of
    # evaluations of the parametric curve does not take the width into
    # consideration, so we have to manually increase this parameter.
    path3.parametric(
        lambda u: (5 * u, 0),
        lambda u: (1, 0),
        final_width=lambda u: 0.4 + 0.1 * numpy.cos(10 * numpy.pi * u),
        number_of_evaluations=256)

    # Add a circular turn and a final tapering segment.
    path3.turn(3, 'l')
    path3.segment(2, final_width=1, final_distance=1.5)

    draw(gdspy.Cell('polygonal-only_paths_2').add(path3))

    # Simple Paths
    # Path defined by a sequence of points and stored as a GDSII path
    sp1 = gdspy.SimplePath([(0, 0), (3, 0), (3, 2),
                            (5, 3), (3, 4), (0, 4)], 1,
                            gdsii_path=True)

    # Other construction methods can still be used
    sp1.smooth([(0, 2), (2, 2), (4, 3), (5, 1)], relative=True)

    # Multiple parallel paths separated by 0.5 with different widths,
    # end caps, and joins.  Because of the join specification, they
    # cannot be stared as GDSII paths, only as polygons.
    sp2 = gdspy.SimplePath([(12, 0), (8, 0), (8, 3), (10, 2)],
                           [0.3, 0.2, 0.4], 0.5,
                           ends=['extended', 'flush', 'round'],
                           corners=['bevel', 'miter', 'round'])
    sp2.arc(2, -0.5 * numpy.pi, 0.5 * numpy.pi)
    sp2.arc(1, 0.5 * numpy.pi, 1.5 * numpy.pi)
    draw(gdspy.Cell('simple_paths').add([sp1, sp2]))

    # Simple Paths 1
    # Path corners and end caps can be custom functions.
    # This corner function creates 'broken' joins.
    def broken(p0, v0, p1, v1, p2, w):
        u0, u1, p = gdspy._intersect_lines(p0, v0, p1, v1)
        if u0 <= 0 and u1 >= 0:
            return [p]
        return [p0, p2, p1]

    # This end cap function creates pointy caps.
    def pointy(p0, v0, p1, v1):
        r = 0.5 * numpy.sqrt(numpy.sum((p0 - p1)**2))
        v0 /= numpy.sqrt(numpy.sum(v0**2))
        v1 /= numpy.sqrt(numpy.sum(v1**2))
        return [p0, 0.5 * (p0 + p1) + 0.5 * (v0 - v1) * r, p1]

    # Paths with arbitrary offsets from the center and multiple layers.
    sp3 = gdspy.SimplePath([(0, 0), (0, 1)], [0.1, 0.3, 0.5],
                           offset=[-0.2, 0, 0.4], layer=[0, 1, 2],
                           corners=broken, ends=pointy)
    sp3.segment((3, 3), offset=[-0.5, -0.1, 0.5])
    sp3.segment((4, 1), width=[0.2, 0.2, 0.2], offset=[-0.2, 0, 0.2])
    sp3.segment((0, -1), relative=True)
    draw(gdspy.Cell('simple_paths_1').add(sp3))

    # Simple Paths 2
    # Corner function to generate bends
    def bend(radius, tolerance=0.01):
        '''
        Corner function generator for bends
        '''
        def _f(p0, v0, p1, v1, p2, w):
            half_w = 0.5 * w
            v2 = p0 - p2
            direction = v0[0] * v1[1] - v0[1] * v1[0]
            if direction == 0:
                return [0.5 * (p0 + p1)]
            elif direction > 0:
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
                r = radius - half_w
                d = (r + half_w) * numpy.tan(abs(a1 - a0) / 2)
            else:
                r = radius + half_w
                d = (r - half_w) * numpy.tan(abs(a1 - a0) / 2)
            np = max(2, 1 + int(0.5 * abs(a1 - a0) / numpy.arccos(1 - tolerance / r) + 0.5))
            angles = numpy.linspace(a0, a1, np)
            points = r * numpy.vstack((numpy.cos(angles), numpy.sin(angles))).T
            return points - points[0] + p0 - d * v0 / (v0[0]**2 + v0[1]**2)**0.5
        return _f

    # Path created with automatic bends of radius 5
    points = [(0, 0), (0, 10), (20, 0), (18, 15), (8, 15)]
    sp4 = gdspy.SimplePath(points, 0.5, corners=bend(5))

    # Same path, generated with natural corners, for comparison
    sp5 = gdspy.SimplePath(points, 0.5, layer=1)
    draw(gdspy.Cell('simple_paths_2').add([sp4, sp5]))

    # Lazy Paths
    # Create 4 parallel paths in different layers
    lp = gdspy.LazyPath((50, 0), [2, 0.5, 1, 1], [0, 0, -1, 1],
                        ends=['extended', 'round', 'flush', 'flush'],
                        layer=[0, 2, 1, 1])
    lp.segment((45, 0))
    lp.segment((5, 0),
               width=[lambda u: 2 + 16 * u * (1 - u), 0.5, 1, 1],
               offset=[0,
                       lambda u: 8 * u * (1 - u) * numpy.cos(12 * numpy.pi * u),
                       lambda u: -1 - 8 * u * (1 - u),
                       lambda u: 1 + 8 * u * (1 - u)])
    lp.segment((0, 0))
    lp.smooth([(5, 10)], angles=[0.5 * numpy.pi, 0],
              width=0.5, offset=[-0.25, 0.25, -0.75, 0.75])
    lp.parametric(lambda u: numpy.array((45 * u, 4 * numpy.sin(6 * numpy.pi * u))),
                  offset=[lambda u: -0.25 * numpy.cos(24 * numpy.pi * u),
                          lambda u: 0.25 * numpy.cos(24 * numpy.pi * u),
                          -0.75, 0.75])
    draw(gdspy.Cell('lazy_paths').add(lp))
