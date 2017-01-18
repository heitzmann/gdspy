import numpy
import gdspy


def test_8b_f():
    f = gdspy._eight_byte_real_to_float
    assert f(b'\x00\x00\x00\x00\x00\x00\x00\x00') == 0
    assert f(b'\x41\x10\x00\x00\x00\x00\x00\x00') == 1
    assert f(b'\x41\x20\x00\x00\x00\x00\x00\x00') == 2
    assert f(b'\xC1\x30\x00\x00\x00\x00\x00\x00') == -3


def test_f_8b():
    g = gdspy._eight_byte_real
    assert b'\x00\x00\x00\x00\x00\x00\x00\x00' == g(0)
    assert b'\x41\x10\x00\x00\x00\x00\x00\x00' == g(1)
    assert b'\x41\x20\x00\x00\x00\x00\x00\x00' == g(2)
    assert b'\xC1\x30\x00\x00\x00\x00\x00\x00' == g(-3)


def test_twoway():
    f = gdspy._eight_byte_real_to_float
    g = gdspy._eight_byte_real
    for x in [0, 1.5, -numpy.pi, 1/3.0e12, -1.0e12/7, 1.1e75, -0.9e-78]:
        assert x == f(g(x))
