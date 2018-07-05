import timeit


a, b = 30, 60
num = 8000000000

t_python = timeit.Timer("python_evaluate.my_evaluate(%f,%f)" % (a,b),"import python_evaluate")
t_cython = timeit.Timer("cython_evaluate.my_evaluate(%f,%f)" % (a,b),"import cython_evaluate")
print "python function", t_python.timeit(10000), "sec"
print "cython function", t_cython.timeit(10000), "sec"