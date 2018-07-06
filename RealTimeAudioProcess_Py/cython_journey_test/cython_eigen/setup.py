from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


setup(
   name = 'array_synthesis_test',
   ext_modules=[
      Extension("_array_synthesis",
	         sources=["array_process.pyx",],
	         extra_compile_args=['-std=c++11',],
	         include_dirs = ['./inc/', '/home/xzc/eigen3.3.4/'],
	         language="c++"),
   ],
   cmdclass = {'build_ext':build_ext},
)
