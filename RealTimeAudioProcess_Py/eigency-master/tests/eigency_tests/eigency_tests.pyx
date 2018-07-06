# distutils: language = c++
# distutils: sources = eigency_tests/eigency_tests_cpp.cpp

from eigency.core cimport *
# cimport eigency.conversions
# from eigency_tests.eigency cimport *


# import eigency
# include "../eigency.pyx"
 
cdef extern from "eigency_tests/eigency_tests_cpp.h":

     cdef long _function_w_vec_arg "function_w_vec_arg"(Map[VectorXd] &)

     cdef long _function_w_1darr_arg "function_w_1darr_arg"(Map[ArrayXi] &)

     cdef void _function_w_vec_arg_no_map1 "function_w_vec_arg_no_map1"(Map[VectorXd])

     cdef void _function_w_vec_arg_no_map2 "function_w_vec_arg_no_map2"(Map[VectorXd] &)

     cdef void _function_w_mat_arg "function_w_mat_arg"(Map[MatrixXd] &)

     cdef void _function_w_complex_mat_arg "function_w_complex_mat_arg"(Map[MatrixXcd] &)

     cdef void _function_w_fullspec_arg "function_w_fullspec_arg" (FlattenedMap[Array, double, Dynamic, _1] &)

     cdef VectorXd _function_w_vec_retval "function_w_vec_retval" ()

     cdef Matrix3d _function_w_mat_retval "function_w_mat_retval" ()

     cdef PlainObjectBase _function_w_mat_retval_full_spec "function_w_mat_retval_full_spec" ()

     cdef Map[ArrayXXd] &_function_filter1 "function_filter1" (Map[ArrayXXd] &)

     cdef PlainObjectBase _function_filter2 "function_filter2" (FlattenedMapWithOrder[Array, double, Dynamic, Dynamic, RowMajor])

     cdef PlainObjectBase _function_filter3 "function_filter3" (FlattenedMapWithStride[Array, double, Dynamic, Dynamic, ColMajor, Unaligned, _1, Dynamic])

     cdef PlainObjectBase _function_type_double "function_type_double" (Map[ArrayXXd] &)
     cdef PlainObjectBase _function_type_float "function_type_float" (Map[ArrayXXf] &)
     cdef PlainObjectBase _function_type_long "function_type_long" (FlattenedMap[Array, long, Dynamic, Dynamic] &)
     cdef PlainObjectBase _function_type_ulong "function_type_ulong" (FlattenedMap[Array, unsigned long, Dynamic, Dynamic] &)
     cdef PlainObjectBase _function_type_int "function_type_int" (Map[ArrayXXi] &)
     cdef PlainObjectBase _function_type_uint "function_type_uint" (FlattenedMap[Array, unsigned int, Dynamic, Dynamic] &)
     cdef PlainObjectBase _function_type_short "function_type_short" (FlattenedMap[Array, short, Dynamic, Dynamic] &)
     cdef PlainObjectBase _function_type_ushort "function_type_ushort" (FlattenedMap[Array, ushort, Dynamic, Dynamic] &)
     cdef PlainObjectBase _function_type_signed_char "function_type_signed_char" (FlattenedMap[Array, signed char, Dynamic, Dynamic] &)
     cdef PlainObjectBase _function_type_unsigned_char "function_type_unsigned_char" (FlattenedMap[Array, unsigned char, Dynamic, Dynamic] &)
     cdef PlainObjectBase _function_type_complex_double "function_type_complex_double" (Map[ArrayXXcd] &)
     cdef PlainObjectBase _function_type_complex_float "function_type_complex_float" (Map[ArrayXXcf] &)

     cdef PlainObjectBase _function_single_col_matrix "function_single_col_matrix" (Map[ArrayXXd] &)

     cdef cppclass _FixedMatrixClass "FixedMatrixClass":
         _FixedMatrixClass () except +
         Matrix3d &get_matrix()
         const Matrix3d &get_const_matrix()

     cdef cppclass _DynamicArrayClass "DynamicArrayClass":
         _DynamicArrayClass (Map[ArrayXXd] &) except +
         ArrayXXd &get_array()
         ArrayXXd get_array_copy()

     cdef cppclass _DynamicRowMajorArrayClass "DynamicRowMajorArrayClass":
         _DynamicRowMajorArrayClass (FlattenedMapWithOrder[Array, double, Dynamic, Dynamic, RowMajor] &) except +
         PlainObjectBase &get_array()
         PlainObjectBase get_array_copy()

# Function with vector argument. 
def function_w_vec_arg(np.ndarray[np.float64_t] array):
    return _function_w_vec_arg(Map[VectorXd](array))

# Function with vector argument. 
def function_w_1darr_arg(np.ndarray[np.int32_t] array):
    return _function_w_1darr_arg(Map[ArrayXi](array))

# Function with vector argument - no map - no ref. 
def function_w_vec_arg_no_map1(np.ndarray[np.float64_t] array):
    return _function_w_vec_arg_no_map1(Map[VectorXd](array))

# Function with vector argument - no map. 
def function_w_vec_arg_no_map2(np.ndarray[np.float64_t] array):
    return _function_w_vec_arg_no_map2(Map[VectorXd](array))

# Function with matrix argument. 
def function_w_mat_arg(np.ndarray[np.float64_t, ndim=2] array):
    return _function_w_mat_arg(Map[MatrixXd](array))

# Function with complex matrix argument. 
def function_w_complex_mat_arg(np.ndarray[np.complex128_t, ndim=2] array):
    return _function_w_complex_mat_arg(Map[MatrixXcd](array))

# Function using a full Map specification, rather than the convenience typedefs
# Note that since cython does not support nested fused types, the Map has been
# flattened to include all arguments at once
def function_w_fullspec_arg(np.ndarray[np.float64_t] array):
    return _function_w_fullspec_arg(FlattenedMap[Array, double, Dynamic, _1](array))

# Function returning vector (copy is made)
def function_w_vec_retval():
    return ndarray(_function_w_vec_retval())

# Function returning matrix (copy is made)
def function_w_mat_retval():
    return ndarray(_function_w_mat_retval())

# Function returning matrix (copy is made)
def function_w_mat_retval_full_spec():
    return ndarray(_function_w_mat_retval_full_spec())

# Function both taking array as argument and returning it
def function_filter1(np.ndarray array):
    return ndarray(_function_filter1(Map[ArrayXXd](array)))

# Function both taking array as argument and returning it - RowMajor order
def function_filter2(np.ndarray[np.float64_t, ndim=2] array):
    return ndarray(_function_filter2(FlattenedMapWithOrder[Array, double, Dynamic, Dynamic, RowMajor](array)))

# Function both taking array as argument and returning it - RowMajor stride
def function_filter3(np.ndarray[np.float64_t, ndim=2] array):
    return ndarray(_function_filter3(FlattenedMapWithStride[Array, double, Dynamic, Dynamic, ColMajor, Unaligned, _1, Dynamic](array)))

# Functions with different matrix types: float64
def function_type_float64(np.ndarray[np.float64_t, ndim=2] array):
    return ndarray(_function_type_double(Map[ArrayXXd](array)))

# Functions with different matrix types: float32
def function_type_float32(np.ndarray[np.float32_t, ndim=2] array):
    return ndarray(_function_type_float(Map[ArrayXXf](array)))

# Functions with different matrix types: long
def function_type_long(np.ndarray[np.long_t, ndim=2] array):
    return ndarray(_function_type_long(FlattenedMap[Array, long, Dynamic, Dynamic](array)))

# Functions with different matrix types: ulong
def function_type_ulong(np.ndarray[np.ulong_t, ndim=2] array):
    return ndarray(_function_type_ulong(FlattenedMap[Array, ulong, Dynamic, Dynamic](array)))

# Functions with different matrix types: int
def function_type_intc(np.ndarray[np.int32_t, ndim=2] array):
    return ndarray(_function_type_int(Map[ArrayXXi](array)))

# Functions with different matrix types: uint
def function_type_uintc(np.ndarray[np.uint32_t, ndim=2] array):
    return ndarray(_function_type_uint(FlattenedMap[Array, uint, Dynamic, Dynamic](array)))

# Functions with different matrix types: short
def function_type_short(np.ndarray[np.int16_t, ndim=2] array):
    return ndarray(_function_type_short(FlattenedMap[Array, short, Dynamic, Dynamic](array)))

# Functions with different matrix types: ushort
def function_type_ushort(np.ndarray[np.uint16_t, ndim=2] array):
    return ndarray(_function_type_ushort(FlattenedMap[Array, ushort, Dynamic, Dynamic](array)))

# Functions with different matrix types: signed char
def function_type_int8(np.ndarray[np.int8_t, ndim=2] array):
    return ndarray(_function_type_signed_char(FlattenedMap[Array, schar, Dynamic, Dynamic](array)))

# Functions with different matrix types: unsigned char
def function_type_uint8(np.ndarray[np.uint8_t, ndim=2] array):
    return ndarray(_function_type_unsigned_char(FlattenedMap[Array, uchar, Dynamic, Dynamic](array)))

# Functions with different matrix types: complex128
def function_type_complex128(np.ndarray[np.complex128_t, ndim=2] array):
    return ndarray(_function_type_complex_double(Map[ArrayXXcd](array)))

# Functions with different matrix types: complex64
def function_type_complex64(np.ndarray[np.complex64_t, ndim=2] array):
    return ndarray(_function_type_complex_float(Map[ArrayXXcf](array)))

# Functions testing a matrix with only one column
def function_single_col_matrix(np.ndarray[np.float64_t, ndim=2] array):
    return ndarray(_function_single_col_matrix(Map[ArrayXXd](array)))

# Functions testing that map properly holds a reference to python objects.
def function_map_holds_reference(np.ndarray[np.float64_t, ndim=2] array):
    # Hold a reference to a copy of an array.
    cdef Map[ArrayXXd] eigency_map = Map[ArrayXXd](array.copy(order="K"))

    # Do some nontrivial operation so that array_copy might be clobbered.
    array_doubled = 2 * array

    # Use the reference to the copy held by eigency_map.
    return ndarray(_function_type_double(eigency_map))


cdef class FixedMatrixClass:
    cdef _FixedMatrixClass *thisptr;
    def __cinit__(self):
        self.thisptr = new _FixedMatrixClass()
    def __dealloc__(self):
        del self.thisptr
    def get_matrix(self):
        return ndarray(self.thisptr.get_matrix())
    def get_const_matrix(self):
        return ndarray(self.thisptr.get_const_matrix())
    def get_const_matrix_force_view(self):
        return ndarray_view(self.thisptr.get_const_matrix())

cdef class DynamicArrayClass:
    cdef _DynamicArrayClass *thisptr;
    def __cinit__(self, np.ndarray[np.float64_t, ndim=2] array):
        self.thisptr = new _DynamicArrayClass(Map[ArrayXXd](array))
    def __dealloc__(self):
        del self.thisptr
    def get_array(self):
        return ndarray(self.thisptr.get_array())

cdef class DynamicRowMajorArrayClass:
    cdef _DynamicRowMajorArrayClass *thisptr;
    def __cinit__(self, np.ndarray[np.float64_t, ndim=2] array):
        self.thisptr = new _DynamicRowMajorArrayClass(FlattenedMapWithOrder[Array, double, Dynamic, Dynamic, RowMajor](array))
    def __dealloc__(self):
        del self.thisptr
    def get_array(self):
        return ndarray(self.thisptr.get_array())
    def get_array_copy(self):
        return ndarray(self.thisptr.get_array_copy())
