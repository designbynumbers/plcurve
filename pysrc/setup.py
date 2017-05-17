import os, sys
import os.path
import shlex

CFLAG = []
CFLAG_INCLUDE = []
if 'CFLAGS' in os.environ or 'CPPFLAGS' in os.environ:
    flags = shlex.split(os.environ.get('CFLAGS',"") +" "+ os.environ.get('CPPFLAGS',""))
    if 'CPPFLAGS' in os.environ:
        del os.environ['CPPFLAGS']
    for flag in flags:
        if flag[:2] == "-I":
            CFLAG_INCLUDE.append(flag[2:])
        else:
            CFLAG.append(flag)
    os.environ['CFLAGS'] = " ".join(CFLAG)

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

USE_CYTHON = None
if "--no-cython" in sys.argv:
    USE_CYTHON = False
    sys.argv.remove("--no-cython")

if USE_CYTHON is None:
    try:
        import Cython
        print "Using Cython version "+Cython.__version__
        USE_CYTHON = True
    except ImportError:
        USE_CYTHON = False

if USE_CYTHON:
    #from Cython.Distutils.extension import Extension
    import Cython.Distutils.build_ext as _build_ext
    from Cython.Build import cythonize
    yes_or_no_cythonize = cythonize
else:
    #from setuptools.extension import Extension
    from setuptools.command. build_ext import build_ext as _build_ext
    yes_or_no_cythonize = no_cythonize

import numpy


ROOT_DIR = os.path.abspath(os.path.dirname(__file__))

LIBPL_PDCODE_PYX = (
    "components.pyx",
    "diagram.pyx",
    "homfly.pyx",
    "isomorphism.pyx",
    "pdstor.pyx",
    "plctopology.pyx",
    "pd_invariants.pyx",
    "planarmap.pyx",)

LIBPL_PYX = (
    "pdstor.pyx",
    "plcurve.pyx",
    "tsmcmc.pyx")

INCLUDE_DIRS = [
    "../src",
    ".",
    os.path.join(numpy.__path__[0], "core", "include")]
LIBRARY_DIRS = ["../src/.libs"]
PLCURVE_LIBRARIES = ["gsl", "plCurve", "gsl", "gslcblas"]
PDCODE_LIBRARIES = ["gsl", "plCurve", "gsl", "gslcblas", "planarmap"]

LIBPL_PDCODE_EXTENSIONS = [
    Extension("libpl.pdcode.%s"%os.path.splitext(pyx_filename)[0],
              [os.path.join(ROOT_DIR, "libpl", "pdcode", pyx_filename)],
              include_dirs=INCLUDE_DIRS, library_dirs=LIBRARY_DIRS,
              libraries=PDCODE_LIBRARIES) for
    pyx_filename in LIBPL_PDCODE_PYX]

LIBPL_EXTENSIONS = [
    Extension("libpl.%s"%os.path.splitext(pyx_filename)[0],
              [os.path.join(ROOT_DIR, "libpl", pyx_filename)],
              include_dirs=INCLUDE_DIRS, library_dirs=LIBRARY_DIRS,
              libraries=PLCURVE_LIBRARIES) for
    pyx_filename in LIBPL_PYX]

EXTENSIONS = LIBPL_EXTENSIONS + LIBPL_PDCODE_EXTENSIONS

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

class build_ext(_build_ext, object):
    def finalize_options(self):
        super(build_ext, self).finalize_options()
        self.include_dirs.extend(CFLAG_INCLUDE)

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
            "build_ext": build_ext,
            "cythonize": CythonizeSources,
        },
        test_suite="tests",
    )
