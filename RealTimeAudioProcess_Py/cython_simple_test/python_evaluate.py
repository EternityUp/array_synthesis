import math


def my_evaluate(a, b):
    x = math.pi / 180.0
    c = math.sin(a * x) + math.cos(b * x)
    r = math.sin(c * a) + math.cos(c * b)
    return r

