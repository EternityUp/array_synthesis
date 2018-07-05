#! /usr/bin/python
# -*- encoding=utf-8 -*-
import example
print example.fact(4)

print example.cvar.My_variable

# Set the value of a C global variable
example.cvar.density = 0.8442
print example.cvar.density
# Use in a math operation
example.cvar.density = example.cvar.density * 1.10
print example.cvar.density

print example.cvar.path
print example.PI, example.VERSION, example.LAGER
print example.FOO

v = example.Vector()

v.x = 3.5
v.y = 7.2
print v.x, v.y, v.z

u = example.Foo()
print u.xx, u.yy, u.path
# u.xx = 2  # 不可更改
u.yy = 3
print u.xx, u.yy, u.path

# example.cvar.path = "/home/xzc/git-hub"  # 不可修改的常量

b = example.Bar()
print b, b.x

c = example.Bar()
c.x = b.x  
print b.x


