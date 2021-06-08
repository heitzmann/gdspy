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
import warnings

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
import struct

from gdspy.gdsiiformat import _eight_byte_real


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
    properties : {integer: string} dictionary
        Properties for these elements.

    Examples
    --------
    >>> label = gdspy.Label('Sample label', (10, 0), 'sw')
    >>> myCell.add(label)
    """

    _anchor = {
        "nw": 0,
        "top left": 0,
        "upper left": 0,
        "n": 1,
        "top center": 1,
        "upper center": 1,
        "ne": 2,
        "top right": 2,
        "upper right": 2,
        "w": 4,
        "middle left": 4,
        "o": 5,
        "middle center": 5,
        "e": 6,
        "middle right": 6,
        "sw": 8,
        "bottom left": 8,
        "lower left": 8,
        "s": 9,
        "bottom center": 9,
        "lower center": 9,
        "se": 10,
        "bottom right": 10,
        "lower right": 10,
    }

    __slots__ = (
        "layer",
        "texttype",
        "text",
        "position",
        "anchor",
        "rotation",
        "magnification",
        "x_reflection",
        "properties",
    )

    def __init__(
        self,
        text,
        position,
        anchor="o",
        rotation=None,
        magnification=None,
        x_reflection=False,
        layer=0,
        texttype=0,
    ):
        self.layer = layer
        self.text = text
        self.position = numpy.array(position)
        self.anchor = Label._anchor.get(anchor.lower(), None)
        if self.anchor is None:
            raise ValueError(
                "[GDSPY] Label anchors must be one of: '"
                + "', '".join(Label._anchor.keys())
                + "'."
            )
        self.rotation = rotation
        self.magnification = magnification
        self.x_reflection = x_reflection
        self.texttype = texttype
        self.properties = {}

    def __repr__(self):
        return 'Label("{0}", ({1[0]}, {1[1]}), {2}, {3}, {4}, {5}, {6})'.format(
            self.text,
            self.position,
            self.rotation,
            self.magnification,
            self.x_reflection,
            self.layer,
            self.texttype,
        )

    def __str__(self):
        return 'Label ("{0}", at ({1[0]}, {1[1]}), rotation {2}, magnification {3}, reflection {4}, layer {5}, texttype {6})'.format(
            self.text,
            self.position,
            self.rotation,
            self.magnification,
            self.x_reflection,
            self.layer,
            self.texttype,
        )

    def to_gds(self, outfile, multiplier):
        """
        Convert this label to a GDSII structure.

        Parameters
        ----------
        outfile : open file
            Output to write the GDSII.
        multiplier : number
            A number that multiplies all dimensions written in the GDSII
            structure.
        """
        outfile.write(
            struct.pack(
                ">4Hh2Hh2Hh",
                4,
                0x0C00,
                6,
                0x0D02,
                self.layer,
                6,
                0x1602,
                self.texttype,
                6,
                0x1701,
                self.anchor,
            )
        )
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
        text = self.text
        if len(text) % 2 != 0:
            text = text + "\0"
        outfile.write(
            struct.pack(
                ">2H2l2H",
                12,
                0x1003,
                int(round(self.position[0] * multiplier)),
                int(round(self.position[1] * multiplier)),
                4 + len(text),
                0x1906,
            )
        )
        outfile.write(text.encode("ascii"))
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
        transform = "scale(1 -1) translate({} {})".format(
            numpy.format_float_positional(
                scaling * self.position[0], trim="0", precision=precision
            ),
            numpy.format_float_positional(
                -scaling * self.position[1], trim="0", precision=precision
            ),
        )
        if self.rotation is not None:
            transform += " rotate({})".format(
                numpy.format_float_positional(
                    -self.rotation, trim="0", precision=precision
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
        ta = ["start", "middle", "end"][self.anchor % 4]
        da = ["text-before-edge", "central", "text-after-edge"][self.anchor // 4]
        if sys.version_info.major < 3:
            text = (
                self.text.decode("utf8")
                .translate({38: "&amp;", 60: "&lt;", 62: "&gt;"})
                .encode("utf8")
            )
        else:
            text = self.text.translate({38: "&amp;", 60: "&lt;", 62: "&gt;"})
        outfile.write(
            '<text class="l{}t{}" text-anchor="{}" dominant-baseline="{}" '
            'transform="{}">{}</text>\n'.format(
                self.layer, self.texttype, ta, da, transform, text
            )
        )

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
