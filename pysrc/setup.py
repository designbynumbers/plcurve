from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("*", ["libpl/*.pyx"],
              include_dirs = ["../src","."],
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
    package_data={"libpl": ["*.pxd"]},
    ext_modules = cythonize(extensions),
)
