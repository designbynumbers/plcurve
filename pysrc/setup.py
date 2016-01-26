import os, sys
import os.path

try:
    from setuptools import setup, Command, Extension
except ImportError:
    print "libpl requires setuptools in order to build"
    sys.exit(1)

def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in ('.pyx', '.py'):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    print extensions
    return extensions

try:
    import Cython
    print "Using Cython version "+Cython.__version__
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

if USE_CYTHON:
    #from Cython.Distutils.extension import Extension
    from Cython.Build import cythonize
    yes_or_no_cythonize = cythonize
else:
    #from setuptools.extension import Extension
    yes_or_no_cythonize = no_cythonize

import numpy


ROOT_DIR = os.path.abspath(os.path.dirname(__file__))

LIBPL_PDCODE_PYX = (
    "components.pyx",
    "diagram.pyx",
    "homfly.pyx",
    "isomorphism.pyx",
    "pdstor.pyx",
    "plctopology.pyx")

LIBPL_PYX = (
    "pdstor.pyx",
    "plcurve.pyx",
    "tsmcmc.pyx")

EXTENSIONS = [
    Extension("pdcode.*", [os.path.join(ROOT_DIR, "libpl", "pdcode", pyx_filename) for
                           pyx_filename in LIBPL_PDCODE_PYX],
              include_dirs=["../src", ".", os.path.join(numpy.__path__[0], "core", "include")],
              library_dirs=["../src/.libs"],
              libraries=["gsl", "plCurve", "gsl", "gslcblas", "planarmap"]),

    Extension("*", [os.path.join(ROOT_DIR, "libpl", pyx_filename) for
                    pyx_filename in LIBPL_PYX],
              include_dirs=["../src", ".", os.path.join(numpy.__path__[0], "core", "include")],
              library_dirs=["../src/.libs"],
              libraries=["gsl", "plCurve", "gsl", "gslcblas"]),
]

class CythonizeSources(Command):
    """Custom command which only Cythonizes .pyx into .c (and doesn't
    build). Useful for, e.g., ``make dist``

    """

    description = "Cythonize sources (i.e. compile .pyx to .c)"

    user_options = [
    ]

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        cythonize(EXTENSIONS)

if __name__ == "__main__":    
    setup(
        name="libpl",
        version="1.0.0",
        author="Harrison Chapman",
        author_email="hchapman@uga.edu",
        description=("Python bindings for Jason Cantarella's libplCurve, "
                     "a piecewise-linear curve and planar diagram library."),
        keywords="piecewise linear curves planar diagrams knot theory",
        packages=["libpl", "libpl.pdcode"],
        url="http://jasoncantarella.com/wordpress/software/plcurve",
        package_dir={'libpl': os.path.join(ROOT_DIR, 'libpl'),
                     'libpl.pdcode': os.path.join(ROOT_DIR, 'libpl', 'pdcode')},
        package_data={"libpl": ["*.pxd",
                                os.path.join("data", "*.txt"),
                                os.path.join("pdcode", "*.pxd")]},
        ext_modules=yes_or_no_cythonize(EXTENSIONS),
        cmdclass={
            "cythonize": CythonizeSources,
        },
        test_suite="tests",
    )
