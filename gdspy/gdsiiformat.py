######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
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

import struct
import numpy
import hashlib


def _record_reader(stream):
    """
    Generator for complete records from a GDSII stream file.

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
        size, rec_type = struct.unpack(">HH", header)
        data_type = rec_type & 0x00FF
        rec_type = rec_type // 256
        data = None
        if size > 4:
            if data_type == 0x01:
                data = numpy.array(
                    struct.unpack(
                        ">{0}H".format((size - 4) // 2), stream.read(size - 4)
                    ),
                    dtype="uint",
                )
            elif data_type == 0x02:
                data = numpy.array(
                    struct.unpack(
                        ">{0}h".format((size - 4) // 2), stream.read(size - 4)
                    ),
                    dtype="int",
                )
            elif data_type == 0x03:
                data = numpy.array(
                    struct.unpack(
                        ">{0}l".format((size - 4) // 4), stream.read(size - 4)
                    ),
                    dtype="int",
                )
            elif data_type == 0x05:
                data = numpy.array(
                    [
                        _eight_byte_real_to_float(stream.read(8))
                        for _ in range((size - 4) // 8)
                    ]
                )
            else:
                data = stream.read(size - 4)
                if str is not bytes:
                    if data[-1] == 0:
                        data = data[:-1].decode("ascii")
                    else:
                        data = data.decode("ascii")
                elif data[-1] == "\0":
                    data = data[:-1]
        yield [rec_type, data]


def _raw_record_reader(stream):
    """
    Generator for complete records from a GDSII stream file.

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
        size, rec_type = struct.unpack(">HH", header)
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
        return b"\x00\x00\x00\x00\x00\x00\x00\x00"
    if value < 0:
        byte1 = 0x80
        value = -value
    else:
        byte1 = 0x00
    fexp = numpy.log2(value) / 4
    exponent = int(numpy.ceil(fexp))
    if fexp == exponent:
        exponent += 1
    mantissa = int(value * 16.0 ** (14 - exponent))
    byte1 += exponent + 64
    byte2 = mantissa // 281474976710656
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
    short1, short2, long3 = struct.unpack(">HHL", value)
    exponent = (short1 & 0x7F00) // 256 - 64
    mantissa = (
        ((short1 & 0x00FF) * 65536 + short2) * 4294967296 + long3
    ) / 72057594037927936.0
    if short1 & 0x8000:
        return -mantissa * 16.0 ** exponent
    return mantissa * 16.0 ** exponent


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
    with open(filename, "rb") as fin:
        data = fin.read()
    contents = []
    start = pos = 0
    while pos < len(data):
        size, rec = struct.unpack(">HH", data[pos : pos + 4])
        if rec == 0x0502:
            start = pos + 28
        elif rec == 0x0700:
            contents.append(data[start:pos])
        pos += size
    h = hashlib.sha1() if engine is None else engine
    for x in sorted(contents):
        h.update(x)
    return h.hexdigest()
