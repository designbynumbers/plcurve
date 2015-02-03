from distutils.core import setup
from Cython.Distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import os

ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
import Cython

extensions = [
    Extension("pdcode.*", [os.path.join(ROOT_DIR, "libpl/pdcode/*.pyx")],
              include_dirs = ["../src",".",os.path.join(numpy.__path__[0], "core", "include")],
              library_dirs = ["../src/.libs"],
              libraries = ["gsl", "plCurve", "gsl", "gslcblas"]),

    Extension("*", [os.path.join(ROOT_DIR, "libpl/*.pyx")],
              include_dirs = ["../src",".",os.path.join(numpy.__path__[0], "core", "include")],
              library_dirs = ["../src/.libs"],
              libraries = ["gsl", "plCurve", "gsl", "gslcblas"]),
]

print "Using Cython version "+Cython.__version__
setup(
    name = "libpl",
    packages = ["libpl"],
    package_dir={'libpl': os.path.join(ROOT_DIR, 'libpl')},
    package_data={"libpl": ["*.pxd", "data/*.txt",]},
    ext_modules = cythonize(extensions),
)
