# -*- encoding=utf-8 -*-

import numpy as np
import sys
import eigen
from ctypes import *

# 查看python模块默认安装路径
print sys.path


A = np.random.random((2000, 50))
B = eigen.MatrixXd(A)

n = np.linalg.norm(B) # Implicit conversion to numpy object

vector_1D = eigen.VectorXd(3)
vector_1D.coeff(1, 2)
print type(vector_1D)
print vector_1D


matrix_1D = eigen.MatrixXd(3, 1)
matrix_1D.setZero()
print matrix_1D






