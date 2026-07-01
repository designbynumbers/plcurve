#!/bin/sh
#
# autogen.sh : regenerate the autotools build system (configure, Makefile.in, config.h.in).
#
# Run this once after a fresh git checkout of the source, then ./configure && make. Because
# releases are shipped as pure source tarballs (no generated files), packagers run this too;
# it needs autoconf/automake/libtool/pkg-config but NONE of libplCurve's own dependencies.
#
set -e
autoreconf --install --verbose
echo
echo "Bootstrap complete. Now run: ./configure && make && make check"
