# Run as:
#    python setup.py build_ext --inplace

import sys
sys.path.insert(0, "..")

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

#ext_modules = cythonize("primes.pyx")#, exclude="numpy_*.pyx")

# Only compile the following if numpy is installed.
try:
    from numpy.distutils.misc_util import get_numpy_include_dirs
    #print help(Extension)
    numpy_demo = [Extension("*",
                            ["listings_loops_cython.pyx"],
                            include_dirs=get_numpy_include_dirs())]
    ext_modules = cythonize(numpy_demo)
except ImportError:
    pass

setup(
  name = 'Demos',
  ext_modules = ext_modules,
)
