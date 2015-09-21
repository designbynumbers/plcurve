from setuptools import setup
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
              libraries = ["gsl", "plCurve", "gsl", "gslcblas", "planarmap"]),

    Extension("*", [os.path.join(ROOT_DIR, "libpl/*.pyx")],
              include_dirs = ["../src",".",os.path.join(numpy.__path__[0], "core", "include")],
              library_dirs = ["../src/.libs"],
              libraries = ["gsl", "plCurve", "gsl", "gslcblas"]),
]

print "Using Cython version "+Cython.__version__
setup(
    name = "libpl",
    packages = ["libpl", "libpl.pdcode"],
    package_dir={'libpl': os.path.join(ROOT_DIR, 'libpl'), 'libpl.pdcode': os.path.join(ROOT_DIR, 'libpl', 'pdcode')},
    package_data={"libpl": ["*.pxd",
                            os.path.join("data", "*.txt"),
                            os.path.join("pdcode", "*.pxd")]},
    ext_modules = cythonize(extensions),
)
