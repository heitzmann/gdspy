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

    # Python 2 doesn't have the pathlib module.
    Path = basestring

else:
    from pathlib import Path

    # Python 3 doesn't have basestring, as unicode is type string
    # Python 2 doesn't equate unicode to string, but both are basestring
    # Now isinstance(s, basestring) will be True for any python version
    basestring = str

import numpy
import datetime
import struct
import itertools
import colorsys
import warnings
import copy as libcopy

from gdspy.polygon import PolygonSet, Polygon
from gdspy.path import FlexPath, RobustPath
from gdspy.label import Label
from gdspy.gdsiiformat import (
    _record_reader,
    _raw_record_reader,
    _eight_byte_real,
    _eight_byte_real_to_float,
)

_mpone = numpy.array((-1.0, 1.0))

use_current_library = True
"""
Globally disable add newly-created cells to the current_library.
"""


class Cell(object):
    """
    Collection of polygons, paths, labels and raferences to other cells.

    .. deprecated:: 1.5
       The parameter `exclude_from_current` has been deprecated
       alongside the use of a global library.  It will be removed in a
       future version of Gdspy.

    Parameters
    ----------
    name : string
        The name of the cell.

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

    __slots__ = (
        "name",
        "polygons",
        "paths",
        "labels",
        "references",
        "_bb_valid",
        "_bounding_box",
    )

    def __init__(self, name, exclude_from_current=False):
        self.name = name
        self.polygons = []
        self.paths = []
        self.labels = []
        self.references = []
        self._bb_valid = False
        self._bounding_box = None
        if use_current_library and not exclude_from_current:
            import gdspy

            gdspy.current_library.add(self, include_dependencies=False)

    def __str__(self):
        return 'Cell ("{}", {} polygons, {} paths, {} labels, {} references)'.format(
            self.name,
            len(self.polygons),
            len(self.paths),
            len(self.labels),
            len(self.references),
        )

    def __iter__(self):
        return itertools.chain(self.polygons, self.paths, self.labels, self.references)

    def to_gds(self, outfile, multiplier, timestamp=None):
        """
        Convert this cell to a GDSII structure.

        Parameters
        ----------
        outfile : open file
            Output to write the GDSII.
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            structure.
        timestamp : datetime object
            Sets the GDSII timestamp.  If None, the current time is
            used.
        """
        now = datetime.datetime.today() if timestamp is None else timestamp
        name = self.name
        if len(name) % 2 != 0:
            name = name + "\0"
        outfile.write(
            struct.pack(
                ">2H12h2H",
                28,
                0x0502,
                now.year,
                now.month,
                now.day,
                now.hour,
                now.minute,
                now.second,
                now.year,
                now.month,
                now.day,
                now.hour,
                now.minute,
                now.second,
                4 + len(name),
                0x0606,
            )
        )
        outfile.write(name.encode("ascii"))
        for polygon in self.polygons:
            polygon.to_gds(outfile, multiplier)
        for path in self.paths:
            path.to_gds(outfile, multiplier)
        for label in self.labels:
            label.to_gds(outfile, multiplier)
        for reference in self.references:
            reference.to_gds(outfile, multiplier)
        outfile.write(struct.pack(">2H", 4, 0x0700))

    def copy(
        self,
        name,
        deep_copy=False,
        translation=None,
        rotation=None,
        scale=None,
        x_reflection=False,
    ):
        """
        Create a copy of this cell.

        Parameters
        ----------
        name : string
            The name of the cell.
        deep_copy : bool
            If False, the new cell will contain only references to the
            existing elements.  If True, copies of all elements are also
            created.  If any transformation is performed, this argument
            is automatically set to True.
        translation : Numpy array[2]
            Amount to translate the cell contents.
        rotation : number
            Rotation angle (in *radians*).
        scale : number
            Scaling factor.
        x_reflection : bool
            Reflect the geometry accros the x axis.

        Returns
        -------
        out : `Cell`
            The new copy of this cell.
        """
        new_cell = Cell(name)

        transform = False
        if (
            x_reflection
            or scale is not None
            or rotation is not None
            or translation is not None
        ):
            transform = True
            deep_copy = True

        if not deep_copy:
            new_cell.polygons = list(self.polygons)
            new_cell.paths = list(self.paths)
            new_cell.labels = list(self.labels)
            new_cell.references = list(self.references)
            return new_cell

        new_cell.polygons = libcopy.deepcopy(self.polygons)
        new_cell.paths = libcopy.deepcopy(self.paths)
        new_cell.labels = libcopy.deepcopy(self.labels)
        new_cell.references = [libcopy.copy(ref) for ref in self.references]

        if transform:
            r = -1 if x_reflection else 1
            s = 1 if scale is None else scale
            t = 0 if rotation is None else rotation
            dx, dy = (0, 0) if translation is None else translation
            ct = numpy.cos(t)
            st = numpy.sin(t)

            for poly in new_cell.polygons:
                if x_reflection:
                    poly.scale(1, -1)
                if scale is not None:
                    poly.scale(scale)
                if rotation is not None:
                    poly.rotate(rotation)
                if translation is not None:
                    poly.translate(dx, dy)

            for path in new_cell.paths:
                path.transform(translation, rotation, scale, x_reflection)

            for lbl in new_cell.labels:
                r0 = -1 if lbl.x_reflection is None else 1
                s0 = 1 if lbl.magnification is None else lbl.magnification
                t0 = 0 if lbl.rotation is None else (lbl.rotation * numpy.pi / 180)
                dx0, dy0 = lbl.position
                lbl.position = (
                    dx + s * (dx0 * ct - r * dy0 * st),
                    dy + s * (dx0 * st + r * dy0 * ct),
                )
                lbl.rotation = (r * t0 + t) * 180 / numpy.pi
                if lbl.rotation == 0:
                    lbl.rotation = None
                lbl.magnification = s * s0
                if lbl.magnification == 1:
                    lbl.magnification = None
                lbl.x_reflection = r * r0 < 0

            for ref in new_cell.references:
                r0 = -1 if ref.x_reflection is None else 1
                s0 = 1 if ref.magnification is None else ref.magnification
                t0 = 0 if ref.rotation is None else (ref.rotation * numpy.pi / 180)
                dx0, dy0 = ref.origin
                ref.origin = (
                    dx + s * (dx0 * ct - r * dy0 * st),
                    dy + s * (dx0 * st + r * dy0 * ct),
                )
                ref.rotation = (r * t0 + t) * 180 / numpy.pi
                if ref.rotation == 0:
                    ref.rotation = None
                ref.magnification = s * s0
                if ref.magnification == 1:
                    ref.magnification = None
                ref.x_reflection = r * r0 < 0

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
                    raise ValueError(
                        "[GDSPY] Only instances of `PolygonSet`, `FlexPath`, "
                        "`RobustPath`, `Label`, `CellReference`, and "
                        "`CellArray` can be added to `Cell`."
                    )
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
        filtered_polys = []
        for element in self.polygons:
            pld = [(poly, l, dt) for poly, l, dt in zip(element.polygons, element.layers, element.datatypes)
                   if not test(poly, l, dt)]
            if len(pld) == 0:
                pass  # we don't need no empty polygons!
            else:
                polys, layers, dts = zip(*pld)
                element.polygons = polys
                element.layers = layers
                element.datatypes = dts
                filtered_polys.append(element)
        self.polygons = filtered_polys
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

    def get_texttypes(self):
        """
        Return the set of texttypes in this cell.

        Returns
        -------
        out : set
            Set of the texttypes used in this cell.
        """
        texttypes = set()
        for reference in self.references:
            texttypes.update(reference.ref_cell.get_textypes())
        for label in self.labels:
            texttypes.add(label.texttype)
        return texttypes

    def get_svg_classes(self):
        """
        Return the set of classes for the SVG representation of this
        cell.

        Returns
        -------
        out0, out1 : sets of 2-tuples
            Sets of (layer, datatype) and (layer, texttype) used in
            this cell.
        """
        ld = set()
        lt = set()
        for element in itertools.chain(self.polygons, self.paths):
            ld.update(zip(element.layers, element.datatypes))
        for label in self.labels:
            lt.add((label.layer, label.texttype))
        for reference in self.references:
            ref_cell = reference.ref_cell
            if isinstance(ref_cell, Cell):
                ref = ref_cell.get_svg_classes()
                ld.update(ref[0])
                lt.update(ref[1])
        return ld, lt

    def get_bounding_box(self):
        """
        Calculate the bounding box for this cell.

        Returns
        -------
        out : Numpy array[2, 2] or None
            Bounding box of this cell [[x_min, y_min], [x_max, y_max]],
            or None if the cell is empty.
        """
        deps_still_valid = all(ref._bb_valid for ref in self.get_dependencies(True))
        cached_bbox_still_valid = self._bb_valid and deps_still_valid
        if not cached_bbox_still_valid:
            bb = numpy.array(((1e300, 1e300), (-1e300, -1e300)))
            all_polygons = []
            for polygon in self.polygons:
                all_polygons.extend(polygon.polygons)
            for path in self.paths:
                all_polygons.extend(path.to_polygonset().polygons)
            for reference in self.references:
                reference_bb = reference.get_bounding_box()
                if reference_bb is not None:
                    all_polygons.append(reference_bb)
            if len(all_polygons) > 0:
                all_points = numpy.concatenate(all_polygons).transpose()
                bb[0, 0] = min(bb[0, 0], all_points[0].min())
                bb[0, 1] = min(bb[0, 1], all_points[1].min())
                bb[1, 0] = max(bb[1, 0], all_points[0].max())
                bb[1, 1] = max(bb[1, 1], all_points[1].max())
                self._bounding_box = bb
            else:
                self._bounding_box = None
            self._bb_valid = True

        if self._bounding_box is None:
            return None
        else:
            # return a *copy* of the cached bounding box to ensure it doesn't get inadvertently modified
            return numpy.array(self._bounding_box)

    def get_polygons(self, by_spec=False, depth=None):
        """
        Return a list of polygons in this cell.

        Parameters
        ----------
        by_spec : bool or tuple
            If True, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype), which
            are used as keys.  If set to a tuple of (layer, datatype),
            only polygons with that specification are returned.
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve polygons.  References below this level will result
            in a bounding box.  If `by_spec` is True the key will be the
            name of this cell.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary with with the list of polygons (if
            `by_spec` is True).

        Note
        ----
        Instances of `FlexPath` and `RobustPath` are also included in
        the result by computing their polygonal boundary.
        """
        if depth is not None and depth < 0:
            if not (by_spec is False or by_spec is True):
                return []
            bb = self.get_bounding_box()
            if bb is None:
                return {} if by_spec else []
            pts = [
                numpy.array(
                    [
                        (bb[0, 0], bb[0, 1]),
                        (bb[0, 0], bb[1, 1]),
                        (bb[1, 0], bb[1, 1]),
                        (bb[1, 0], bb[0, 1]),
                    ]
                )
            ]
            polygons = {self.name: pts} if by_spec else pts
        else:
            if by_spec is True:
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
            elif by_spec is False:
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
            else:
                polygons = []
                layer, datatype = by_spec
                polygons.extend(
                    numpy.array(polyset.polygons[ii])
                    for polyset in self.polygons
                    for ii in range(len(polyset.polygons))
                    if polyset.layers[ii] == layer and polyset.datatypes[ii] == datatype
                )

                for path in self.paths:
                    if any(ld == by_spec for ld in zip(path.layers, path.datatypes)):
                        path_polygons = path.get_polygons(True)
                        if by_spec in path_polygons:
                            polygons.extend(path_polygons[by_spec])
                for reference in self.references:
                    if depth is None:
                        next_depth = None
                    else:
                        next_depth = depth - 1
                    polygons.extend(reference.get_polygons(by_spec, next_depth))
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

    def get_labels(self, depth=None, set_transform=False):
        """
        Return a list with a copy of the labels in this cell.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve labels from.
        set_transform : bool
            If True, labels will include the transformations from
            the references they are from.

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
                labels.extend(reference.get_labels(next_depth, set_transform))
        return labels

    def get_dependencies(self, recursive=False):
        """
        Return a set of the cells included in this cell as references.

        Parameters
        ----------
        recursive : bool
            If True returns cascading dependencies.

        Returns
        -------
        out : set of `Cell`
            Set of the cells referenced by this cell.
        """
        dependencies = set()
        for reference in self.references:
            if isinstance(reference.ref_cell, Cell):
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

    def to_svg(self, outfile, scaling, precision, attributes):
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
        attributes : string
            Additional attributes to set for the cell group.
        """
        outfile.write('<g id="')
        outfile.write(self.name.replace("#", "_"))
        outfile.write('" ')
        outfile.write(attributes)
        outfile.write(">\n")
        for polygon in self.polygons:
            polygon.to_svg(outfile, scaling, precision)
        for path in self.paths:
            path.to_svg(outfile, scaling, precision)
        for label in self.labels:
            label.to_svg(outfile, scaling, precision)
        for reference in self.references:
            reference.to_svg(outfile, scaling, precision)
        outfile.write("</g>\n")

    def write_svg(
        self,
        outfile,
        scaling=10,
        style=None,
        fontstyle=None,
        background="#222",
        pad="5%",
        precision=None,
    ):
        """
        Export this cell to an SVG file.

        The dimensions actually written on the GDSII file will be the
        dimensions of the objects created times the ratio
        unit/precision.  For example, if a circle with radius 1.5 is
        created and we set `GdsLibrary.unit` to 1.0e-6 (1 um) and
        `GdsLibrary.precision` to 1.0e-9` (1 nm), the radius of the
        circle will be 1.5 um and the GDSII file will contain the
        dimension 1500 nm.

        Parameters
        ----------
        outfile : file, string or Path
            The file (or path) where the GDSII stream will be written.
            It must be opened for writing operations in binary format.
        scaling : number
            Scaling factor for the geometry.
        style : dict or None
            Dictionary indexed by (layer, datatype) tuples.  Entries
            must be dictionaries with CSS key-value pairs for the
            presentation attributes of the geometry in that layer and
            datatype.
        fontstyle : dict or None
            Dictionary indexed by (layer, texttype) tuples.  Entries
            must be dictionaries with CSS key-value pairs for the
            presentation attributes of the labels in that layer and
            texttype.
        background : string or None
            String specifying the background color.  If None, no
            background is inserted.
        pad : number or string
            Background margin around the cell bounding box.  It can
            be specified as a percentage of the width or height,
            whichever is the largest.
        precision : positive integer or `None`
            Maximal number of digits for coordinates after scaling.

        Examples
        --------
        >>> cell = gdspy.Cell('MAIN')
        >>> cell.add(gdspy.Rectangle((0, 0), (10, 10), layer=1))
        >>> # Define fill and stroke for layer 1 and datatype 0
        >>> mystyle = {(1, 0): {'fill': '#CC00FF',
                                'stroke': 'black'}}
        >>> cell.write_svg('main.svg', style=mystyle)
        """
        bb = self.get_bounding_box()
        if bb is None:
            return
        close = True
        if hasattr(outfile, "__fspath__"):
            outfile = open(outfile.__fspath__(), "w")
        elif isinstance(outfile, (basestring, Path)):
            outfile = open(outfile, "w")
        else:
            close = False
        if style is None:
            style = {}
        if fontstyle is None:
            fontstyle = {}
        bb *= scaling
        x = bb[0, 0]
        y = -bb[1, 1]
        w = bb[1, 0] - bb[0, 0]
        h = bb[1, 1] - bb[0, 1]
        if background is not None:
            if isinstance(pad, basestring):
                if pad[-1] == "%":
                    pad = max(w, h) * float(pad[:-1]) / 100
                else:
                    pad = float(pad)
            x -= pad
            y -= pad
            w += 2 * pad
            h += 2 * pad
        outfile.write(
            """<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
     width="{}" height="{}" viewBox="{} {} {} {}">
<defs>
<style type="text/css">
""".format(
                numpy.format_float_positional(w, trim="0", precision=precision),
                numpy.format_float_positional(h, trim="0", precision=precision),
                numpy.format_float_positional(x, trim="0", precision=precision),
                numpy.format_float_positional(y, trim="0", precision=precision),
                numpy.format_float_positional(w, trim="0", precision=precision),
                numpy.format_float_positional(h, trim="0", precision=precision),
            )
        )
        ldkeys, ltkeys = self.get_svg_classes()
        for k in ldkeys:
            l, d = k
            if k in style:
                style_dict = style[k]
            else:
                c = "rgb({}, {}, {})".format(
                    *[
                        int(255 * c + 0.5)
                        for c in colorsys.hsv_to_rgb(
                            (l % 3) / 3.0 + (l % 6 // 3) / 6.0 + (l // 6) / 11.0,
                            1 - ((l + d) % 8) / 12.0,
                            1 - (d % 3) / 4.0,
                        )
                    ]
                )
                style_dict = {"stroke": c, "fill": c, "fill-opacity": "0.5"}
            outfile.write(".l{}d{} {{".format(l, d))
            outfile.write(" ".join("{}: {};".format(*x) for x in style_dict.items()))
            outfile.write("}\n")
        for k in ltkeys:
            l, t = k
            if k in fontstyle:
                style_dict = fontstyle[k]
            else:
                c = "rgb({}, {}, {})".format(
                    *[
                        int(255 * c + 0.5)
                        for c in colorsys.hsv_to_rgb(
                            (l % 3) / 3.0 + (l % 6 // 3) / 6.0 + (l // 6) / 11.0,
                            1 - ((l + t) % 8) / 12.0,
                            1 - (t % 3) / 4.0,
                        )
                    ]
                )
                style_dict = {"stroke": "none", "fill": c}
            outfile.write(".l{}t{} {{".format(l, t))
            outfile.write(" ".join("{}: {};".format(*x) for x in style_dict.items()))
            outfile.write("}\n")
        outfile.write("</style>\n")
        for cell in self.get_dependencies(True):
            cell.to_svg(outfile, scaling, precision, "")
        outfile.write("</defs>\n")
        if background is not None:
            outfile.write(
                '<rect x="{}" y="{}" width="{}" height="{}" fill="{}" stroke="none"/>\n'.format(
                    numpy.format_float_positional(x, trim="0", precision=precision),
                    numpy.format_float_positional(y, trim="0", precision=precision),
                    numpy.format_float_positional(w, trim="0", precision=precision),
                    numpy.format_float_positional(h, trim="0", precision=precision),
                    background,
                )
            )
        self.to_svg(outfile, scaling, precision, 'transform="scale(1 -1)"')
        outfile.write("</svg>")
        if close:
            outfile.close()


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

    Attributes
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
    properties : {integer: string} dictionary
        Properties for these elements.
    """

    __slots__ = (
        "ref_cell",
        "origin",
        "rotation",
        "magnification",
        "x_reflection",
        "properties",
    )

    def __init__(
        self,
        ref_cell,
        origin=(0, 0),
        rotation=None,
        magnification=None,
        x_reflection=False,
        ignore_missing=False,
    ):
        self.origin = origin
        self.ref_cell = ref_cell
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection
        self.properties = {}
        if not isinstance(self.ref_cell, Cell) and not ignore_missing:
            warnings.warn(
                "[GDSPY] Cell {0} not found; operations on this "
                "CellReference may not work.".format(self.ref_cell),
                stacklevel=2,
            )

    def __str__(self):
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        return 'CellReference ("{0}", at ({1[0]}, {1[1]}), rotation {2}, magnification {3}, reflection {4})'.format(
            name, self.origin, self.rotation, self.magnification, self.x_reflection
        )

    def __repr__(self):
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        return 'CellReference("{0}", ({1[0]}, {1[1]}), {2}, {3}, {4})'.format(
            name, self.origin, self.rotation, self.magnification, self.x_reflection
        )

    def to_gds(self, outfile, multiplier):
        """
        Convert this object to a GDSII element.

        Parameters
        ----------
        outfile : open file
            Output to write the GDSII.
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            element.
        """
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        if len(name) % 2 != 0:
            name = name + "\0"
        outfile.write(struct.pack(">4H", 4, 0x0A00, 4 + len(name), 0x1206))
        outfile.write(name.encode("ascii"))
        if (
            (self.rotation is not None)
            or (self.magnification is not None)
            or self.x_reflection
        ):
            word = 0
            values = b""
            if self.x_reflection:
                word += 0x8000
            if not (self.magnification is None):
                # This flag indicates that the magnification is absolute, not
                # relative (not supported).
                # word += 0x0004
                values += struct.pack(">2H", 12, 0x1B05) + _eight_byte_real(
                    self.magnification
                )
            if not (self.rotation is None):
                # This flag indicates that the rotation is absolute, not
                # relative (not supported).
                # word += 0x0002
                values += struct.pack(">2H", 12, 0x1C05) + _eight_byte_real(
                    self.rotation
                )
            outfile.write(struct.pack(">3H", 6, 0x1A01, word))
            outfile.write(values)
        outfile.write(
            struct.pack(
                ">2H2l",
                12,
                0x1003,
                int(round(self.origin[0] * multiplier)),
                int(round(self.origin[1] * multiplier)),
            )
        )
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
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        transform = "translate({} {})".format(
            numpy.format_float_positional(
                scaling * self.origin[0], trim="0", precision=precision
            ),
            numpy.format_float_positional(
                scaling * self.origin[1], trim="0", precision=precision
            ),
        )
        if self.rotation is not None:
            transform += " rotate({})".format(
                numpy.format_float_positional(
                    self.rotation, trim="0", precision=precision
                )
            )
        if self.x_reflection:
            transform += " scale(1 -1)"
        if self.magnification is not None:
            transform += " scale({})".format(
                numpy.format_float_positional(
                    self.magnification, trim="0", precision=precision
                )
            )
        outfile.write('<use transform="')
        outfile.write(transform)
        outfile.write('" xlink:href="#')
        outfile.write(name.replace("#", "_"))
        outfile.write('"/>\n')

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
                factor = self.magnification ** 2
                cell_area = self.ref_cell.area(True)
                for kk in cell_area.keys():
                    cell_area[kk] *= factor
                return cell_area
            else:
                return self.ref_cell.area() * self.magnification ** 2

    def _transform_polygons(self, polygons):
        """
        Transform a set of polygons.

        This reference transformation is used to transform the given
        polygons in place.

        Parameters
        ----------
        polygons : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary of lists of polygons.

        Returns
        -------
        polygons : list of array-like[N][2] or dictionary
            Transformed polygons. Same object as `polygons` argument.
        """
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0) * _mpone
        if self.x_reflection:
            xrefl = numpy.array((1, -1))
        if self.magnification is not None:
            mag = numpy.array((self.magnification, self.magnification), dtype=float)
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        if isinstance(polygons, dict):
            for kk in polygons.keys():
                for ii in range(len(polygons[kk])):
                    if self.x_reflection:
                        polygons[kk][ii] = polygons[kk][ii] * xrefl
                    if self.magnification is not None:
                        polygons[kk][ii] = polygons[kk][ii] * mag
                    if self.rotation is not None:
                        polygons[kk][ii] = (
                            polygons[kk][ii] * ct + polygons[kk][ii][:, ::-1] * st
                        )
                    if self.origin is not None:
                        polygons[kk][ii] = polygons[kk][ii] + orgn
        else:
            for ii in range(len(polygons)):
                if self.x_reflection:
                    polygons[ii] = polygons[ii] * xrefl
                if self.magnification is not None:
                    polygons[ii] = polygons[ii] * mag
                if self.rotation is not None:
                    polygons[ii] = polygons[ii] * ct + polygons[ii][:, ::-1] * st
                if self.origin is not None:
                    polygons[ii] = polygons[ii] + orgn
        return polygons

    def get_polygons(self, by_spec=False, depth=None):
        """
        Return the list of polygons created by this reference.

        Parameters
        ----------
        by_spec : bool or tuple
            If True, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype).
            If set to a tuple of (layer, datatype), only polygons
            with that specification are returned.
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
        polygons = self.ref_cell.get_polygons(by_spec, depth)
        return self._transform_polygons(polygons)

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
            xrefl = numpy.array((1, -1))
        if self.magnification is not None:
            mag = numpy.array((self.magnification, self.magnification), dtype=float)
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
                    ps.polygons[ii] = (
                        ps.polygons[ii] * ct + ps.polygons[ii][:, ::-1] * st
                    )
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
        return [
            p.transform(trans, rot, self.magnification, self.x_reflection)
            for p in self.ref_cell.get_paths(depth=depth)
        ]

    def get_labels(self, depth=None, set_transform=False):
        """
        Return the list of labels created by this reference.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve labels from.
        set_transform : bool
            If True, labels will include the transformations from
            the reference.

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
            xrefl = numpy.array((1, -1))
        if self.magnification is not None:
            mag = numpy.array((self.magnification, self.magnification), dtype=float)
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        labels = self.ref_cell.get_labels(depth=depth, set_transform=set_transform)
        for lbl in labels:
            if self.x_reflection:
                lbl.position = lbl.position * xrefl
            if self.magnification is not None:
                lbl.position = lbl.position * mag
            if self.rotation is not None:
                lbl.position = lbl.position * ct + lbl.position[::-1] * st
            if self.origin is not None:
                lbl.position = lbl.position + orgn
            if set_transform:
                if self.magnification is not None:
                    if lbl.magnification is not None:
                        lbl.magnification *= self.magnification
                    else:
                        lbl.magnification = self.magnification
                if self.x_reflection is not None:
                    lbl.x_reflection = not lbl.x_reflection
                if self.rotation is not None:
                    if lbl.rotation is not None:
                        lbl.rotation += self.rotation
                    else:
                        lbl.rotation = self.rotation
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

        if self.rotation is None or self.rotation % 90 == 0:
            cell_bbox = self.ref_cell.get_bounding_box()
            if cell_bbox is None:
                return None
            polygons = self._transform_polygons([cell_bbox])
        else:
            # For non-cardinal rotations of a reference, we must use the
            # flattened polygons for the reference
            polygons = self.get_polygons()
        if len(polygons) == 0:
            bb = None
        else:
            all_points = numpy.concatenate(polygons).transpose()
            bb = numpy.array(
                (
                    (all_points[0].min(), all_points[1].min()),
                    (all_points[0].max(), all_points[1].max()),
                )
            )
        return bb

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

    Attributes
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
    properties : {integer: string} dictionary
        Properties for these elements.
    """

    __slots__ = (
        "ref_cell",
        "origin",
        "rotation",
        "magnification",
        "x_reflection",
        "columns",
        "rows",
        "spacing",
        "properties",
    )

    def __init__(
        self,
        ref_cell,
        columns,
        rows,
        spacing,
        origin=(0, 0),
        rotation=None,
        magnification=None,
        x_reflection=False,
        ignore_missing=False,
    ):
        self.columns = columns
        self.rows = rows
        self.spacing = spacing
        self.origin = origin
        self.ref_cell = ref_cell
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection
        self.properties = {}
        if not isinstance(self.ref_cell, Cell) and not ignore_missing:
            warnings.warn(
                "[GDSPY] Cell {0} not found; operations on this "
                "CellArray may not work.".format(self.ref_cell),
                stacklevel=2,
            )

    def __str__(self):
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        return 'CellArray ("{0}", {1} x {2}, at ({3[0]}, {3[1]}), spacing {4[0]} x {4[1]}, rotation {5}, magnification {6}, reflection {7})'.format(
            name,
            self.columns,
            self.rows,
            self.origin,
            self.spacing,
            self.rotation,
            self.magnification,
            self.x_reflection,
        )

    def __repr__(self):
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        return 'CellArray("{0}", {1}, {2}, ({4[0]}, {4[1]}), ({3[0]}, {3[1]}), {5}, {6}, {7})'.format(
            name,
            self.columns,
            self.rows,
            self.origin,
            self.spacing,
            self.rotation,
            self.magnification,
            self.x_reflection,
        )

    def to_gds(self, outfile, multiplier):
        """
        Convert this object to a GDSII element.

        Parameters
        ----------
        outfile : open file
            Output to write the GDSII.
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            element.
        """
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        if len(name) % 2 != 0:
            name = name + "\0"
        outfile.write(struct.pack(">4H", 4, 0x0B00, 4 + len(name), 0x1206))
        outfile.write(name.encode("ascii"))
        x2 = self.origin[0] + self.columns * self.spacing[0]
        y2 = self.origin[1]
        x3 = self.origin[0]
        y3 = self.origin[1] + self.rows * self.spacing[1]
        if (
            (self.rotation is not None)
            or (self.magnification is not None)
            or self.x_reflection
        ):
            word = 0
            values = b""
            if self.x_reflection:
                word += 0x8000
                y3 = 2 * self.origin[1] - y3
            if not (self.magnification is None):
                # This flag indicates that the magnification is absolute, not
                # relative (not supported).
                # word += 0x0004
                values += struct.pack(">2H", 12, 0x1B05) + _eight_byte_real(
                    self.magnification
                )
            if not (self.rotation is None):
                # This flag indicates that the rotation is absolute, not
                # relative (not supported).
                # word += 0x0002
                sa = numpy.sin(self.rotation * numpy.pi / 180.0)
                ca = numpy.cos(self.rotation * numpy.pi / 180.0)
                tmp = (
                    (x2 - self.origin[0]) * ca
                    - (y2 - self.origin[1]) * sa
                    + self.origin[0]
                )
                y2 = (
                    (x2 - self.origin[0]) * sa
                    + (y2 - self.origin[1]) * ca
                    + self.origin[1]
                )
                x2 = tmp
                tmp = (
                    (x3 - self.origin[0]) * ca
                    - (y3 - self.origin[1]) * sa
                    + self.origin[0]
                )
                y3 = (
                    (x3 - self.origin[0]) * sa
                    + (y3 - self.origin[1]) * ca
                    + self.origin[1]
                )
                x3 = tmp
                values += struct.pack(">2H", 12, 0x1C05) + _eight_byte_real(
                    self.rotation
                )
            outfile.write(struct.pack(">3H", 6, 0x1A01, word))
            outfile.write(values)
        outfile.write(
            struct.pack(
                ">2H2h2H6l",
                8,
                0x1302,
                self.columns,
                self.rows,
                28,
                0x1003,
                int(round(self.origin[0] * multiplier)),
                int(round(self.origin[1] * multiplier)),
                int(round(x2 * multiplier)),
                int(round(y2 * multiplier)),
                int(round(x3 * multiplier)),
                int(round(y3 * multiplier)),
            )
        )
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
        if isinstance(self.ref_cell, Cell):
            name = self.ref_cell.name
        else:
            name = self.ref_cell
        transform = "translate({} {})".format(
            numpy.format_float_positional(
                scaling * self.origin[0], trim="0", precision=precision
            ),
            numpy.format_float_positional(
                scaling * self.origin[1], trim="0", precision=precision
            ),
        )
        if self.rotation is not None:
            transform += " rotate({})".format(
                numpy.format_float_positional(
                    self.rotation, trim="0", precision=precision
                )
            )
        if self.x_reflection:
            transform += " scale(1 -1)"
        mag = (
            ""
            if self.magnification is None
            else " scale({})".format(
                numpy.format_float_positional(
                    self.magnification, trim="0", precision=precision
                )
            )
        )
        for ii in range(self.columns):
            dx = scaling * self.spacing[0] * ii
            for jj in range(self.rows):
                dy = scaling * self.spacing[1] * jj
                outfile.write('<use transform="')
                outfile.write(transform)
                outfile.write(
                    " translate({} {})".format(
                        numpy.format_float_positional(
                            dx, trim="0", precision=precision
                        ),
                        numpy.format_float_positional(
                            dy, trim="0", precision=precision
                        ),
                    )
                )
                outfile.write(mag)
                outfile.write('" xlink:href="#')
                outfile.write(name.replace("#", "_"))
                outfile.write('"/>\n')

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
            factor = self.columns * self.rows * self.magnification ** 2
        if by_spec:
            cell_area = self.ref_cell.area(True)
            for kk in cell_area.keys():
                cell_area[kk] *= factor
            return cell_area
        else:
            return self.ref_cell.area() * factor

    def _transform_polygons(self, polygons):
        """
        Transform a set of polygons.

        This reference transformation is used to transform the given
        polygons.

        Parameters
        ----------
        polygons : list of array-like[N][2] or dictionary
            List containing the coordinates of the vertices of each
            polygon, or dictionary of lists of polygons.

        Returns
        -------
        out : list of array-like[N][2] or dictionary
            Transformed polygons.
        """
        if self.rotation is not None:
            ct = numpy.cos(self.rotation * numpy.pi / 180.0)
            st = numpy.sin(self.rotation * numpy.pi / 180.0) * _mpone
        if self.magnification is not None:
            mag = numpy.array((self.magnification, self.magnification), dtype=float)
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        if self.x_reflection:
            xrefl = numpy.array((1, -1))
        if isinstance(polygons, dict):
            out_polygons = {}
            for kk in polygons.keys():
                out_polygons[kk] = []
                for ii in range(self.columns):
                    for jj in range(self.rows):
                        spc = numpy.array([self.spacing[0] * ii, self.spacing[1] * jj])
                        for points in polygons[kk]:
                            if self.magnification:
                                out_polygons[kk].append(points * mag + spc)
                            else:
                                out_polygons[kk].append(points + spc)
                            if self.x_reflection:
                                out_polygons[kk][-1] = out_polygons[kk][-1] * xrefl
                            if self.rotation is not None:
                                out_polygons[kk][-1] = (
                                    out_polygons[kk][-1] * ct
                                    + out_polygons[kk][-1][:, ::-1] * st
                                )
                            if self.origin is not None:
                                out_polygons[kk][-1] = out_polygons[kk][-1] + orgn
        else:
            out_polygons = []
            for ii in range(self.columns):
                for jj in range(self.rows):
                    spc = numpy.array([self.spacing[0] * ii, self.spacing[1] * jj])
                    for points in polygons:
                        if self.magnification is not None:
                            out_polygons.append(points * mag + spc)
                        else:
                            out_polygons.append(points + spc)
                        if self.x_reflection:
                            out_polygons[-1] = out_polygons[-1] * xrefl
                        if self.rotation is not None:
                            out_polygons[-1] = (
                                out_polygons[-1] * ct + out_polygons[-1][:, ::-1] * st
                            )
                        if self.origin is not None:
                            out_polygons[-1] = out_polygons[-1] + orgn
        return out_polygons

    def get_polygons(self, by_spec=False, depth=None):
        """
        Return the list of polygons created by this reference.

        Parameters
        ----------
        by_spec : bool or tuple
            If True, the return value is a dictionary with the
            polygons of each individual pair (layer, datatype).
            If set to a tuple of (layer, datatype), only polygons
            with that specification are returned.
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
        cell_polygons = self.ref_cell.get_polygons(by_spec, depth)
        return self._transform_polygons(cell_polygons)

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
            xrefl = numpy.array((1, -1))
        if self.magnification is not None:
            mag = numpy.array((self.magnification, self.magnification), dtype=float)
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
                            ps.polygons[ii] = (
                                ps.polygons[ii] * ct + ps.polygons[ii][:, ::-1] * st
                            )
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
                    array.append(
                        libcopy.deepcopy(path).transform(
                            trans, rot, self.magnification, self.x_reflection, spc
                        )
                    )
        return array

    def get_labels(self, depth=None, set_transform=False):
        """
        Return the list of labels created by this reference.

        Parameters
        ----------
        depth : integer or None
            If not None, defines from how many reference levels to
            retrieve labels from.
        set_transform : bool
            If True, labels will include the transformations from
            the reference.

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
            mag = numpy.array((self.magnification, self.magnification), dtype=float)
        if self.origin is not None:
            orgn = numpy.array(self.origin)
        if self.x_reflection:
            xrefl = numpy.array((1, -1))
        cell_labels = self.ref_cell.get_labels(depth=depth, set_transform=set_transform)
        labels = []
        for ii in range(self.columns):
            for jj in range(self.rows):
                spc = numpy.array([self.spacing[0] * ii, self.spacing[1] * jj])
                for clbl in cell_labels:
                    lbl = libcopy.deepcopy(clbl)
                    if self.magnification is not None:
                        lbl.position = lbl.position * mag + spc
                    else:
                        lbl.position = lbl.position + spc
                    if self.x_reflection:
                        lbl.position = lbl.position * xrefl
                    if self.rotation is not None:
                        lbl.position = lbl.position * ct + lbl.position[::-1] * st
                    if self.origin is not None:
                        lbl.position = lbl.position + orgn
                    if set_transform:
                        if self.magnification is not None:
                            if lbl.magnification is not None:
                                lbl.magnification *= self.magnification
                            else:
                                lbl.magnification = self.magnification
                        if self.x_reflection is not None:
                            lbl.x_reflection = not lbl.x_reflection
                        if self.rotation is not None:
                            if lbl.rotation is not None:
                                lbl.rotation += self.rotation
                            else:
                                lbl.rotation = self.rotation
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
        if self.rotation is None or self.rotation % 90 == 0:
            cell_bbox = self.ref_cell.get_bounding_box()
            if cell_bbox is None:
                return None
            polygons = self._transform_polygons([cell_bbox])
        else:
            # For non-cardinal rotations of a reference, we must use the
            # flattened polygons for the reference
            polygons = self.get_polygons()
        if len(polygons) == 0:
            bb = None
        else:
            all_points = numpy.concatenate(polygons).transpose()
            bb = numpy.array(
                (
                    (all_points[0].min(), all_points[1].min()),
                    (all_points[0].max(), all_points[1].max()),
                )
            )
        return bb

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
    cells : dictionary
        Dictionary of cells in this library, indexed by name.
    unit : number
        Unit size for the objects in the library (in *meters*).
    precision : number
        Precision for the dimensions of the objects in the library (in
        *meters*).
    """

    _record_name = (
        "HEADER",
        "BGNLIB",
        "LIBNAME",
        "UNITS",
        "ENDLIB",
        "BGNSTR",
        "STRNAME",
        "ENDSTR",
        "BOUNDARY",
        "PATH",
        "SREF",
        "AREF",
        "TEXT",
        "LAYER",
        "DATATYPE",
        "WIDTH",
        "XY",
        "ENDEL",
        "SNAME",
        "COLROW",
        "TEXTNODE",
        "NODE",
        "TEXTTYPE",
        "PRESENTATION",
        "SPACING",
        "STRING",
        "STRANS",
        "MAG",
        "ANGLE",
        "UINTEGER",
        "USTRING",
        "REFLIBS",
        "FONTS",
        "PATHTYPE",
        "GENERATIONS",
        "ATTRTABLE",
        "STYPTABLE",
        "STRTYPE",
        "ELFLAGS",
        "ELKEY",
        "LINKTYPE",
        "LINKKEYS",
        "NODETYPE",
        "PROPATTR",
        "PROPVALUE",
        "BOX",
        "BOXTYPE",
        "PLEX",
        "BGNEXTN",
        "ENDTEXTN",
        "TAPENUM",
        "TAPECODE",
        "STRCLASS",
        "RESERVED",
        "FORMAT",
        "MASK",
        "ENDMASKS",
        "LIBDIRSIZE",
        "SRFNAME",
        "LIBSECUR",
    )
    _unused_records = (0x05, 0x00, 0x01, 0x02, 0x034, 0x38)
    _import_anchors = ["nw", "n", "ne", None, "w", "o", "e", None, "sw", "s", "se"]
    _pathtype_dict = {0: "flush", 1: "round", 2: "extended"}

    __slots__ = "name", "cells", "unit", "precision", "_references"

    def __init__(
        self, name="library", infile=None, unit=1e-6, precision=1e-9, **kwargs
    ):
        self.name = name
        self.cells = {}
        self.unit = unit
        self.precision = precision
        if infile is not None:
            self.read_gds(infile, **kwargs)

    def __str__(self):
        return "GdsLibrary (" + ", ".join([c for c in self.cells]) + ")"

    def __iter__(self):
        return iter(self.cells.values())

    def new_cell(self, name, overwrite_duplicate=False, update_references=True):
        """
        Create a new cell and add it to this library.

        Parameters
        ----------
        name : string
            Name of the cell.
        overwrite_duplicate : bool
            If True, an existing cell with the same name in the library
            will be overwritten.
        update_references : bool
            If True, `CellReference` and `CellArray` instances from an
            overwritten cell are updated to the new one (used only when
            `overwrite_duplicate` is True).

        Returns
        -------
        out : `Cell`
            The created cell.

        Notes
        -----
        This is equivalent to:
        >>> cell = gdspy.Cell(name)
        >>> lib.add(cell, False, overwrite_duplicate, update_references)
        """
        cell = Cell(name)
        self.add(cell, False, overwrite_duplicate, update_references)
        return cell

    def add(
        self,
        cell,
        include_dependencies=True,
        overwrite_duplicate=False,
        update_references=True,
    ):
        """
        Add one or more cells to the library.

        Parameters
        ----------
        cell : `Cell` or iterable
            Cells to be included in the library.
        include_dependencies : bool
            If True, also add cells referenced by `cell`, recursively.
        overwrite_duplicate : bool
            If True, an existing cell with the same name in the library
            will be overwritten.
        update_references : bool
            If True, `CellReference` and `CellArray` instances from an
            overwritten cell are updated to the new one (used only when
            `overwrite_duplicate` is True).

        Returns
        -------
        out : `GdsLibrary`
            This object.
        """
        if isinstance(cell, Cell):
            cell_set = set([cell])
            if include_dependencies:
                cell_set.update(cell.get_dependencies(True))
        else:
            cell_set = set(cell)
            if include_dependencies:
                for c in cell:
                    cell_set.update(c.get_dependencies(True))
        for c in cell_set:
            if (
                not overwrite_duplicate
                and c.name in self.cells
                and self.cells[c.name] is not c
            ):
                raise ValueError(
                    "[GDSPY] Cell named {0} already present in library.".format(c.name)
                )
            if (
                overwrite_duplicate
                and update_references
                and c.name in self.cells
                and self.cells[c.name] is not c
            ):
                self.replace_references(c.name, c)
            self.cells[c.name] = c
        return self

    def remove(self, cell, remove_references=True):
        """
        Remove a cell from the library.

        Parameters
        ----------
        cell : `Cell` or string
            Cell to be removed from the library.
        remove_references : bool
            If True, `CellReference` and `CellArray` using the removed
            cell will also be removed.

        Returns
        -------
        out : integer
            Number of references removed.
        """
        if isinstance(cell, Cell):
            name = cell.name
        else:
            name = cell
        if name in self.cells:
            del self.cells[name]
        removed = 0
        if remove_references:
            for c in self.cells.values():
                removed += len(c.references)
                c.references = [
                    ref
                    for ref in c.references
                    if name
                    != (
                        ref.ref_cell.name
                        if isinstance(ref.ref_cell, Cell)
                        else ref.ref_cell
                    )
                ]
                removed -= len(c.references)
        return removed

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
        outfile : file, string or Path
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
        close = True
        if hasattr(outfile, "__fspath__"):
            outfile = open(outfile.__fspath__(), "wb")
        elif isinstance(outfile, (basestring, Path)):
            outfile = open(outfile, "wb")
        else:
            close = False
        now = datetime.datetime.today() if timestamp is None else timestamp
        name = self.name if len(self.name) % 2 == 0 else (self.name + "\0")
        outfile.write(
            struct.pack(
                ">5H12h2H",
                6,
                0x0002,
                0x0258,
                28,
                0x0102,
                now.year,
                now.month,
                now.day,
                now.hour,
                now.minute,
                now.second,
                now.year,
                now.month,
                now.day,
                now.hour,
                now.minute,
                now.second,
                4 + len(name),
                0x0206,
            )
            + name.encode("ascii")
            + struct.pack(">2H", 20, 0x0305)
            + _eight_byte_real(self.precision / self.unit)
            + _eight_byte_real(self.precision)
        )
        if cells is None:
            cells = self.cells.values()
        else:
            cells = [self.cells.get(c, c) for c in cells]
        if len(cells) == 0:
            warnings.warn("[GDSPY] Creating a GDSII file without any cells.")
        for cell in cells:
            cell.to_gds(outfile, self.unit / self.precision, timestamp=timestamp)
        if binary_cells is not None:
            for bc in binary_cells:
                outfile.write(bc)
        outfile.write(struct.pack(">2H", 4, 0x0400))
        if close:
            outfile.close()

    def read_gds(
        self,
        infile,
        units="skip",
        rename={},
        rename_template="{name}",
        layers={},
        datatypes={},
        texttypes={},
    ):
        """
        Read a GDSII file into this library.

        Parameters
        ----------
        infile : file, string or Path
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
        close = True
        if hasattr(infile, "__fspath__"):
            infile = open(infile.__fspath__(), "rb")
        elif isinstance(infile, (basestring, Path)):
            infile = open(infile, "rb")
        else:
            close = False
        emitted_warnings = []
        kwargs = {}
        create_element = None
        factor = 1
        cell = None
        properties = {}
        attr = -1
        for record in _record_reader(infile):
            # LAYER
            if record[0] == 0x0D:
                kwargs["layer"] = layers.get(record[1][0], record[1][0])
            # DATATYPE or BOXTYPE
            elif record[0] == 0x0E or record[0] == 0x2E:
                kwargs["datatype"] = datatypes.get(record[1][0], record[1][0])
            # TEXTTYPE
            elif record[0] == 0x16:
                kwargs["texttype"] = texttypes.get(record[1][0], record[1][0])
            # XY
            elif record[0] == 0x10:
                if "xy" in kwargs:
                    kwargs["xy"] = numpy.concatenate((kwargs["xy"], factor * record[1]))
                else:
                    kwargs["xy"] = factor * record[1]
            # WIDTH
            elif record[0] == 0x0F:
                kwargs["width"] = factor * abs(record[1][0])
                if record[1][0] < 0:
                    kwargs["width_transform"] = False
            # ENDEL
            elif record[0] == 0x11:
                if create_element is not None:
                    el = create_element(**kwargs)
                    if len(properties) > 0:
                        el.properties = properties
                        properties = {}
                    cell.add(el)
                    create_element = None
                kwargs = {}
            # BOUNDARY
            elif record[0] == 0x08:
                create_element = self._create_polygon
            # PATH
            elif record[0] == 0x09:
                create_element = self._create_path
            # BOX
            elif record[0] == 0x2D:
                create_element = self._create_polygon
                if record[0] not in emitted_warnings:
                    warnings.warn(
                        "[GDSPY] GDSII elements of type BOX are imported as polygons.",
                        stacklevel=2,
                    )
                    emitted_warnings.append(record[0])
            # TEXT
            elif record[0] == 0x0C:
                create_element = self._create_label
            # SNAME
            elif record[0] == 0x12:
                if record[1] in rename:
                    name = rename[record[1]]
                else:
                    name = rename_template.format(name=record[1])
                kwargs["ref_cell"] = name
            # COLROW
            elif record[0] == 0x13:
                kwargs["columns"] = record[1][0]
                kwargs["rows"] = record[1][1]
            # STRANS
            elif record[0] == 0x1A:
                kwargs["x_reflection"] = (int(record[1][0]) & 0x8000) > 0
                if (int(record[1][0]) & 0x0006) and record[0] not in emitted_warnings:
                    warnings.warn(
                        "[GDSPY] Absolute magnification or rotation of "
                        "references is not supported.  Transformations "
                        "will be interpreted as relative.",
                        stacklevel=2,
                    )
                    emitted_warnings.append(record[0])
            # MAG
            elif record[0] == 0x1B:
                kwargs["magnification"] = record[1][0]
            # ANGLE
            elif record[0] == 0x1C:
                kwargs["rotation"] = record[1][0]
            # SREF
            elif record[0] == 0x0A:
                create_element = self._create_reference
            # AREF
            elif record[0] == 0x0B:
                create_element = self._create_array
            # STRNAME
            elif record[0] == 0x06:
                if record[1] in rename:
                    name = rename[record[1]]
                else:
                    name = rename_template.format(name=record[1])
                cell = Cell(name, exclude_from_current=True)
                if name in self.cells:
                    raise ValueError("[GDSPY] Multiple cells with name: {0} in GDSII file".format(name))
                self.cells[name] = cell
            # STRING
            elif record[0] == 0x19:
                kwargs["text"] = record[1]
            # ENDSTR
            elif record[0] == 0x07:
                cell = None
            # UNITS
            elif record[0] == 0x03:
                if units == "skip":
                    factor = record[1][0]
                elif units == "import":
                    self.unit = record[1][1] / record[1][0]
                    self.precision = record[1][1]
                    factor = record[1][0]
                elif units == "convert":
                    factor = record[1][1] / self.unit
                else:
                    raise ValueError(
                        "[GDSPY] units must be one of 'convert', 'import' or 'skip'."
                    )
            # LIBNAME
            elif record[0] == 0x02:
                self.name = record[1]
            # PRESENTATION
            elif record[0] == 0x17:
                kwargs["anchor"] = GdsLibrary._import_anchors[
                    int(record[1][0]) & 0x000F
                ]
            # PATHTYPE
            elif record[0] == 0x21:
                kwargs["ends"] = GdsLibrary._pathtype_dict.get(record[1][0], "extended")
            # BGNEXTN
            elif record[0] == 0x30:
                kwargs["bgnextn"] = factor * record[1][0]
            # ENDEXTN
            elif record[0] == 0x31:
                kwargs["endextn"] = factor * record[1][0]
            # ENDLIB
            elif record[0] == 0x04:
                for ref in self._references:
                    if ref.ref_cell in self.cells:
                        ref.ref_cell = self.cells[ref.ref_cell]
            # PROPATTR
            elif record[0] == 0x2B:
                attr = record[1][0]
            # PROPVALUE
            elif record[0] == 0x2C:
                properties[attr] = record[1]
            # Not supported
            elif (
                record[0] not in emitted_warnings
                and record[0] not in GdsLibrary._unused_records
            ):
                warnings.warn(
                    "[GDSPY] Record type {0} ({1:02X}) is not supported.".format(
                        GdsLibrary._record_name[record[0]], record[0]
                    ),
                    stacklevel=2,
                )
                emitted_warnings.append(record[0])
        if close:
            infile.close()
        return self

    def _create_polygon(self, layer, datatype, xy):
        return Polygon(xy[:-2].reshape((xy.size // 2 - 1, 2)), layer, datatype)

    def _create_path(self, **kwargs):
        xy = kwargs.pop("xy")
        if "bgnextn" in kwargs or "endextn" in kwargs:
            kwargs["ends"] = (kwargs.pop("bgnextn", 0), kwargs.pop("endextn", 0))
        kwargs["points"] = xy.reshape((xy.size // 2, 2))
        kwargs["gdsii_path"] = True
        return FlexPath(**kwargs)

    def _create_label(self, xy, width=None, width_transform=None, ends=None, **kwargs):
        kwargs["position"] = xy
        return Label(**kwargs)

    def _create_reference(self, **kwargs):
        kwargs["origin"] = kwargs.pop("xy")
        kwargs["ignore_missing"] = True
        ref = CellReference(**kwargs)
        ref.ref_cell = kwargs["ref_cell"]
        self._references.append(ref)
        return ref

    def _create_array(self, **kwargs):
        xy = kwargs.pop("xy")
        kwargs["origin"] = xy[0:2]
        if "x_reflection" in kwargs:
            if "rotation" in kwargs:
                sa = -numpy.sin(kwargs["rotation"] * numpy.pi / 180.0)
                ca = numpy.cos(kwargs["rotation"] * numpy.pi / 180.0)
                x2 = (xy[2] - xy[0]) * ca - (xy[3] - xy[1]) * sa + xy[0]
                y3 = (xy[4] - xy[0]) * sa + (xy[5] - xy[1]) * ca + xy[1]
            else:
                x2 = xy[2]
                y3 = xy[5]
            if kwargs["x_reflection"]:
                y3 = 2 * xy[1] - y3
            kwargs["spacing"] = (
                (x2 - xy[0]) / kwargs["columns"],
                (y3 - xy[1]) / kwargs["rows"],
            )
        else:
            kwargs["spacing"] = (
                (xy[2] - xy[0]) / kwargs["columns"],
                (xy[5] - xy[1]) / kwargs["rows"],
            )
        kwargs["ignore_missing"] = True
        ref = CellArray(**kwargs)
        ref.ref_cell = kwargs["ref_cell"]
        self._references.append(ref)
        return ref

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
        top = set(self)
        for cell in self:
            top.difference_update(cell.get_dependencies())
        return list(top)

    def rename_cell(self, cell, name, update_references=True):
        """
        Rename an existing cell in the library.

        Parameters
        ----------
        cell : `Cell` or string
            Cell to be renamed.  It must be present in the library.
        name : string
            New name for the cell.  It cannot be present in the library.
        update_references : bool
            If True, replace references using the old name with the new
            cell.

        Returns
        -------
        out : integer
            Number of updated references.
        """
        if isinstance(cell, Cell):
            old_name = cell.name
            if old_name not in self.cells:
                raise ValueError(
                    "[GDSPY] Cell named {0} not present in library.".format(old_name)
                )
            if self.cells[old_name] is not cell:
                raise ValueError(
                    "[GDSPY] Cell named {0} doesn't match library's.".format(old_name)
                )
        else:
            old_name = cell
            if old_name not in self.cells:
                raise ValueError(
                    "[GDSPY] Cell named {0} not present in library.".format(old_name)
                )
            cell = self.cells[old_name]
        if name in self.cells:
            raise ValueError(
                "[GDSPY] Cell named {0} already present in library.  "
                "Use `add` to overwrite cells.".format(name)
            )
        del self.cells[old_name]
        self.cells[name] = cell
        cell.name = name
        if update_references:
            return self.replace_references(old_name, cell)
        return 0

    def replace_references(self, old_cell, new_cell):
        """
        Replace cells in all references in the library.

        All `CellReference` and `CellArray` using the `old_cell` are
        updated to reference `new_cell`.  Matching with `old_cell` is
        by name only.

        Parameters
        ----------
        old_cell : `Cell` or string
            Cell to be replaced.
        new_cell : `Cell` or string
            Replacement cell.  If the cell name is passed and it is
            present in the library, the actual cell is used instead.

        Returns
        -------
        out : integer
            Number of replacements.
        """
        if isinstance(old_cell, Cell):
            old_name = old_cell.name
        else:
            old_name = old_cell
        if not isinstance(new_cell, Cell) and new_cell in self.cells:
            new_cell = self.cells[new_cell]
        replacements = 0
        for cell in self.cells.values():
            for ref in cell.references:
                if isinstance(ref.ref_cell, Cell):
                    if ref.ref_cell.name == old_name:
                        ref.ref_cell = new_cell
                        replacements += 1
                elif ref.ref_cell == old_name:
                    ref.ref_cell = new_cell
                    replacements += 1
        return replacements

    def extract(self, cell, overwrite_duplicate=False):
        """
        Extract a cell from the this GDSII file and include it in the
        current global library, including referenced dependencies.

        .. deprecated:: 1.5
           `extract` is deprecated and will be removed in a future
           version of Gdspy.  Gdspy no longer uses a global library.

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
        warnings.warn(
            "[GDSPY] extract and the use of the global library is deprecated.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        import gdspy

        cell = self.cells.get(cell, cell)
        gdspy.current_library.add(
            cell, include_dependencies=True, overwrite_duplicate=overwrite_duplicate
        )
        return cell


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
    outfile : file, string or Path
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

    __slots__ = "_outfile", "_close", "_res"

    def __init__(
        self, outfile, name="library", unit=1.0e-6, precision=1.0e-9, timestamp=None
    ):
        self._close = True
        if hasattr(outfile, "__fspath__"):
            self._outfile = open(outfile.__fspath__(), "wb")
        elif isinstance(outfile, (basestring, Path)):
            self._outfile = open(outfile, "wb")
        else:
            self._outfile = outfile
            self._close = False
        self._res = unit / precision
        now = datetime.datetime.today() if timestamp is None else timestamp
        if len(name) % 2 != 0:
            name = name + "\0"
        self._outfile.write(
            struct.pack(
                ">5H12h2H",
                6,
                0x0002,
                0x0258,
                28,
                0x0102,
                now.year,
                now.month,
                now.day,
                now.hour,
                now.minute,
                now.second,
                now.year,
                now.month,
                now.day,
                now.hour,
                now.minute,
                now.second,
                4 + len(name),
                0x0206,
            )
            + name.encode("ascii")
            + struct.pack(">2H", 20, 0x0305)
            + _eight_byte_real(precision / unit)
            + _eight_byte_real(precision)
        )

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
        cell.to_gds(self._outfile, self._res, timestamp)
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
        self._outfile.write(struct.pack(">2H", 4, 0x0400))
        if self._close:
            self._outfile.close()


def get_gds_units(infile):
    """
    Return the unit and precision used in the GDS stream file.

    Parameters
    ----------
    infile : file, string or Path
        GDSII stream file to be queried.

    Returns
    -------
    out : 2-tuple
        Return ``(unit, precision)`` from the file.
    """
    close = True
    if hasattr(infile, "__fspath__"):
        infile = open(infile.__fspath__(), "rb")
    elif isinstance(infile, (basestring, Path)):
        infile = open(infile, "rb")
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
    infile : file, string, or Path
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
    close = True
    if hasattr(infile, "__fspath__"):
        infile = open(infile.__fspath__(), "rb")
    elif isinstance(infile, (basestring, Path)):
        infile = open(infile, "rb")
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
                    name = data[4:-1].decode("ascii")
                else:
                    name = data[4:].decode("ascii")
            else:
                if data[-1] == "\0":
                    name = data[4:-1]
                else:
                    name = data[4:]
        # ENDSTR
        elif rec_type == 0x07:
            cell_data.append(data)
            cells[name] = b"".join(cell_data)
            cell_data = None
        elif cell_data is not None:
            cell_data.append(data)
    if close:
        infile.close()
    return cells
