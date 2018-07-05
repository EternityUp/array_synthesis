cimport cython

cdef extern from "math.h":
    float cosf(float theta)
    float sinf(float theta)


cdef float _my_evaluate(float a, float b):
    cdef float pi = 3.14159265
    cdef float x = pi / 180.0
    cdef float c = sinf(a * x) + cosf(b * x)
    cdef float r = sinf(c*a) + cosf(c*b)
    return r


def my_evaluate(float a,float b):
    cdef float x = _my_evaluate(a, b)
    return x
