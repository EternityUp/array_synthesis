cdef extern from "cpp_matrixxfpy.h":
    cdef cppclass MatrixXf:
        MatrixXf()
    cdef cppclass MatrixXfPy:
        MatrixXfPy()
        MatrixXfPy(int d1, int d2)
        MatrixXfPy(MatrixXf &other)
        int rows()
        int cols()
        float coeff(int, int)