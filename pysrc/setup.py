from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import os

extensions = [
    Extension("*", ["libpl/*.pyx"],
              include_dirs = ["../src",".",os.path.join(numpy.__path__[0], "core", "include")],
              library_dirs = ["../src/.libs"],
              libraries = ["gsl","plCurve"]),
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
    package_data={"libpl": ["*.pxd", "data/*.txt", "util/*.py"]},
    ext_modules = cythonize(extensions),
)
