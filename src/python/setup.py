from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("plcurve", ["plcurve/plcurve.pyx"],
              include_dirs = [".."],
              libraries = ["gsl","plCurve"]),
    Extension("pdcode", ["plcurve/pdcode.pyx"],
              include_dirs = [".."],
              libraries = ["gsl","plCurve"])
]

setup(
    name = "plcurve",
    packages = ["plcurve"],
    package_data={"plcurve": ["*.pxd"]},
    ext_modules = cythonize(extensions),
)
