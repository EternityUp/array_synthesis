# -*- encoding=utf-8 -*-

from ctypes import *

print cdll.LoadLibrary("libc.so.6")

libc = CDLL("libc.so.6")

print dir(libc)
print libc.printf

print libc.time(None)


print c_int()

print c_char_p("Hello, World!")

print c_ushort(-3)

i = c_int(42)
print i, i.value
i.value = -99
print i.value

s = "hello world"
c_s = c_char_p(s)
print c_s.value
c_s.value = "hi there"
print c_s.value

p = create_string_buffer(3)
print sizeof(p), repr(p.raw)

p = create_string_buffer("Hello")
print sizeof(p), repr(p.raw)

p = create_string_buffer("Hello", 10)
print sizeof(p), repr(p.raw)

p.value = "Hi"
print sizeof(p), repr(p.raw)


printf = libc.printf
printf("Hello, %s\n", "World!")
printf("Hello, %S\n", u"World!")
printf("%d bottles of beer\n", 42)
printf("%f bottles of beer\n", c_double(42.5))
printf("An int %d, a double %f\n", 1234, c_double(3.14))

fv = c_float(3.14)
print fv, fv.value


class Bottles(object):
    def __init__(self,number):
        self._as_parameter_ = number


bottles = Bottles(42)
printf("%d bottles of beer\n", bottles)

strchr = libc.strchr
strchr.restype = c_char_p
strchr.argtype = [c_char_p, c_char]
print strchr("abcdef", "d")
print strchr("abcdef", "x")

i = c_int()
f = c_float()
s = create_string_buffer('\000' * 32)
print i.value, f.value, repr(s.value)

libc.sscanf("1 3.14 Hello", "%d %f %s", byref(i), byref(f), s)
print i.value, f.value, repr(s.value)


class POINT(Structure):
    _fields_ = [("x", c_int),
                ("y", c_int)]


point = POINT(10, 20)
print point.x, point.y
point = POINT(y=5)
print point.x, point.y


class POINT3(Structure):
    _fields_ = [("x", c_int),
                ("y", c_int),
                ("z", c_int)]


point = POINT3(1,2,3)
print point.x, point.y, point.z

class RECT(Structure):
    _fields_ = [("upperleft", POINT),
                ("lowerright", POINT)]


rc = RECT(POINT(10, 20))
print rc.upperleft.x, rc.upperleft.y, rc.lowerright.x, rc.lowerright.y

r = RECT(POINT(1, 2), POINT(3, 4))
r = RECT((1, 2), (3, 4))

print POINT.x
print POINT.y


class Int(Structure):
    _fields_ = [("first_16", c_int, 16),
                ("second_16", c_int, 16)]


print Int.first_16
print Int.second_16


TenPointsArrayType = POINT * 10


class MyStruct(Structure):
    _fields_ = [("a", c_int),
                ("b", c_float),
                ("point_array", POINT * 4)]


print len(MyStruct().point_array)

arr = TenPointsArrayType()
for pt in arr:
    print pt.x, pt.y


TenIntegers = c_int * 10
ii = TenIntegers(1,2,3,4,5,6,7,8,9,10)
print ii
for i in ii:
    print i


i = c_int(42)
pi = pointer(i)
print pi.contents

print pi.contents is i

print pi.contents is pi.contents

i = c_int(99)
pi.contents = i
print pi.contents
print pi[0], i
pi[0] = 22
print pi[0], i

PI = POINTER(c_int)
print PI

print PI(c_int(42))

null_ptr = POINTER(c_int)()
print bool(null_ptr)
# null_ptr[0] = 1     # ValueError: NULL pointer access


class Bar(Structure):
    _fields_ = [("count", c_int),
                ("values", POINTER(c_int))]


bar = Bar()
bar.values = (c_int * 3)(1,2,3)
bar.count = 3
for i in range(bar.count):
    print bar.values[i]


bar = Bar()
bar.values = cast((c_byte * 4)(), POINTER(c_int))
print bar.values[0]


IntArray5 = c_int * 5
ia = IntArray5(5,1,7,33,99)
qsort = libc.qsort
qsort.restype = None

CMPFUNC = CFUNCTYPE(c_int, POINTER(c_int), POINTER(c_int))


def py_cmp_func(a, b):
    print "py_cmp_func", a, b
    return 0


cmp_func = CMPFUNC(py_cmp_func)

qsort(ia, len(ia), sizeof(c_int), cmp_func)


def py_cmp_func(a, b):
    print "py_cmp_func", a[0], b[0]
    return 0


cmp_func = CMPFUNC(py_cmp_func)
qsort(ia, len(ia), sizeof(c_int), cmp_func)


def py_cmp_func(a, b):
    print "py_cmp_func", a[0], b[0]
    return a[0] - b[0]


cmp_func = CMPFUNC(py_cmp_func)
qsort(ia, len(ia), sizeof(c_int), cmp_func)


for i in ia:
    print i


short_array = (c_short * 4)()
print sizeof(short_array)

resize(short_array, 32)
print sizeof(short_array)
print short_array

print short_array[:]

from ctypes.util import find_library
print find_library("m")
print find_library("c")


p = byref(c_int(1))

buffer = create_string_buffer("hello")
buffer2 = create_string_buffer("hello world")

print buffer2.value











