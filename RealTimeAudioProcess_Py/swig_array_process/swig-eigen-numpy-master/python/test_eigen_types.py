# -*- coding=utf-8 -*-

import eigen
import pyinverter
import numpy as np

print dir(pyinverter)


inv = pyinverter.Inverter()
print inv
print dir(inv)

m = np.array([[2.0, 0.0],
              [0.0, 2.0]])

m_inv = inv.getInverse(m)
print m_inv
print inv.getInversePlus1(m)


m_inv_list = inv.getInverseList([np.eye(2), np.eye(3)])
print m_inv_list
print m_inv_list[0]
print m_inv_list[1]

print(pyinverter.templatedInverse(np.eye(2)))

print(inv.getInverseRef(4.0 * np.eye(4)))


print pyinverter.inverter_wrapper

print dir(pyinverter.inverter_wrapper)

x = pyinverter.inverter_wrapper.vectorMatrixXd
y = pyinverter.inverter_wrapper.vectorVectorXd

print x
print y

print dir(x)
print dir(y)

x1 = x()
print x1
y1 = y()
print y1

print dir(x1)
print dir(y1)


x2 = x(10)
print x2[0]

print dir(x2[0])


print dir(eigen)

s = eigen.MatrixXd(3, 4)
print s

r = eigen.VectorXd(10)
print r

print dir(r)

r.__setitem__(0, 20)
print r


print type(s)
print type(r)

