from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import os

ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
extensions = [
    Extension("*", [os.path.join(ROOT_DIR, "libpl/*.pyx")],
              include_dirs = ["../src",".",os.path.join(numpy.__path__[0], "core", "include")],
              library_dirs = ["../src/.libs"],
              libraries = ["gsl", "plCurve", "gsl", "gslcblas"]),
    Extension("pdcode.*", [os.path.join(ROOT_DIR, "libpl/pdcode/*.pyx")],
              include_dirs = ["../src",".",os.path.join(numpy.__path__[0], "core", "include")],
              library_dirs = ["../src/.libs"],
              libraries = ["gsl", "plCurve", "gsl", "gslcblas"]),

#    Extension("pdcode", ["libpl/pdcode.pyx"],
#              include_dirs = ["..","."],
#              libraries = ["gsl","plCurve"]),
#    Extension("tsmcmc", ["libpl/plcurve.pyx", "libpl/tsmcmc.pyx",],
#              include_dirs = ["..","."],
#              libraries = ["gsl","plCurve"])
]

setup(
    name = "libpl",
    packages = ["libpl"],
    package_dir={'libpl': os.path.join(ROOT_DIR, 'libpl')},
    package_data={"libpl": ["*.pxd", "data/*.txt",]},
    ext_modules = cythonize(extensions),
)
