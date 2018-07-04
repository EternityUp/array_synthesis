#! /usr/bin/env python
# -*- encoding=utf-8 -*-

import wave
from ctypes import *
import os
import sys
import numpy as np

env_dist = os.environ


if 'LD_LIBRARY_PATH' not in os.environ:
    os.environ['LD_LIBRARY_PATH'] = os.getcwd() + u"/libs/"
    try:
        os.execv(sys.argv[0], sys.argv)
    except Exception, exc:
        print 'Failed re-exec:', exc
        sys.exit(1)

env_dist = os.environ
print env_dist.get('LD_LIBRARY_PATH')


lib = CDLL(os.getcwd() + u"/libs/libprint_2d_ptr.so")


int_array = c_int * 10
print int_array

a = int_array(1, 2 ,3)
print a[0], a[1], a[2], a[3]

x = c_char_p()  # 一级指针
print x, type(x)
x.value = 'x'

y = pointer(x)  # 二级指针
print y, type(y)

print x.value
print y.contents.value


a = np.asarray(range(16), dtype=np.int32).reshape([4,4])
print a.flags

if not a.flags['C_CONTIGUOUS']:
    a = np.ascontiguous(a, dtype=a.dtype)

a_ctypes_ptr = cast(a.ctypes.data, POINTER(c_int))

for i in range(16):
    print a_ctypes_ptr[i]


b = np.array(range(10))
print b, b.dtype
c = b.reshape(2,-1)
print b.reshape(2, -1), c[1][1]
print b.shape, c.shape


x = (c_int * 4)(1,2,3,4)
px = pointer(x)
print px, type(px)
print px.contents[0]


cp = cast(c.ctypes.data, POINTER(c_int))
print cp, type(cp)

# 创建二维指针数组
int4_array = (c_int * 4)  # 首先创建1维数组类型：元素数目为4的int数组
x1 = int4_array(1, 2, 3, 4)
x2 = int4_array(5, 6, 7, 8)
x3 = int4_array(9, 10, 11, 12)

# 创建3个指向先前创建的1维数组的指针
x1_p = pointer(x1)
x2_p = pointer(x2)
x3_p = pointer(x3)
print x1_p.contents[:]
print x2_p.contents[:]
print x3_p.contents[:]

# 创建指向一维数组的指针数组
px = (POINTER(int4_array) * 3)(x1_p, x2_p, x3_p)
print "px:{},type(px):{}".format(px, type(px))
print px[0], x1, px[0].contents[:], x1[:]

for i in range(3):
    print "px[i]:{}".format(px[i])
    print px[i].contents[:]

px1D = POINTER(int4_array)  # 指向1维int数组的指针
px2D = POINTER(px1D)  # 指向（指向1维int数组的指针）的指针

print px1D, px2D

px2D.contents = px

print px2D.contents, px2D.contents[0][0][:]


print_2d_ptr = lib.print_2d_ptr

print_2d_ptr(px, 3, 4)


int_array_frame_len = (c_int * 5)
int_2D_array_ptr = (POINTER(int_array_frame_len) * 4)
int4_5 = np.random.randint(low=0, high=10, size=(4, 5), dtype=np.int32)
print int4_5


print int4_5[0].tolist(), type(int4_5[0].tolist())


ex_2d_ptr0 = pointer(int_array_frame_len(* int4_5[0].tolist()))
ex_2d_ptr1 = pointer(int_array_frame_len(* int4_5[1].tolist()))
ex_2d_ptr2 = pointer(int_array_frame_len(* int4_5[2].tolist()))
ex_2d_ptr3 = pointer(int_array_frame_len(* int4_5[3].tolist()))

ex_2d_ptr = int_2D_array_ptr(ex_2d_ptr0, ex_2d_ptr1, ex_2d_ptr2, ex_2d_ptr3)

print ex_2d_ptr, int_2D_array_ptr

print_2d_ptr(ex_2d_ptr, 4, 5)


ex_2d_ptr = int_2D_array_ptr()
for i in range(4):
    ex_2d_ptr[i] = pointer(int_array_frame_len(* int4_5[i].tolist()))


print_2d_ptr(ex_2d_ptr, 4, 5)







