import os
from os.path import basename, join
from setuptools import setup, find_packages
from setuptools.extension import Extension

import eigency
import numpy as np

__package_name__ = "eigency"
__eigen_dir__ = eigency.__eigen_dir__
__eigen_lib_dir__ = join(basename(__eigen_dir__), 'Eigen')

# Not all users may have cython installed.  If they only want this as a means
# to access the Eigen header files to compile their own C++ code, then they
# may not have cython already installed.  Therefore, only require cython
# for cases where the user will need to build the .cpp files from the .pyx
# files (typically from a git clone) and not for other pip installations.
# cf. discussion in PR #26.

# Follow the pattern recommended here:
# http://cython.readthedocs.io/en/latest/src/reference/compilation.html#distributing-cython-modules
try:
    from Cython.Build import cythonize
    # Maybe make this a command line option?
    USE_CYTHON = True
    ext = '.pyx'
except ImportError:
    USE_CYTHON = False
    ext = '.cpp'

extensions = [
    Extension("eigency.conversions", ["eigency/conversions"+ext],
              include_dirs = [np.get_include(), __eigen_dir__],
              language="c++"
    ),
    Extension("eigency.core", ["eigency/core"+ext],
              include_dirs = [np.get_include(), __eigen_dir__],
              language="c++"
    )
]

if USE_CYTHON:
    extensions = cythonize(extensions)

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except(IOError, ImportError):
    long_description = open('README.md').read()

eigen_data_files = []
for root, dirs, files in os.walk(join(__eigen_dir__, 'Eigen')):
    for f in files:
        if f.endswith('.h'):
            eigen_data_files.append(join(root, f))

dist = setup(
    name = __package_name__,
    description = "Cython interface between the numpy arrays and the Matrix/Array classes of the Eigen C++ library",
    long_description=long_description,
    classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: C++',
        'Programming Language :: Cython',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    license = "MIT",
    author = "Wouter Boomsma",
    author_email = "wb@di.ku.dk",
    url = "https://github.com/wouterboomsma/eigency",
    use_scm_version = False,
    ext_modules = extensions,
    packages = find_packages(),
    include_package_data=True,
    package_data = {__package_name__: [
        '*.h', '*.pxd', '*.pyx',
        join(__eigen_lib_dir__, '*'),
    ] + eigen_data_files},
    exclude_package_data = {__package_name__: [join(__eigen_lib_dir__, 'CMakeLists.txt')]},
    install_requires = ['numpy']
)

