if BUILD_PYTHON

PY_FILES = \
	libpl/__init__.py \
	libpl/data.py \
	libpl/graphs.py \
	libpl/pdcode/__init__.py

PYX_FILES = \
	libpl/plcurve.pyx \
	libpl/tsmcmc.pyx \
	libpl/pdstor.pyx \
	libpl/pdcode/diagram.pyx \
	libpl/pdcode/components.pyx \
	libpl/pdcode/pdstor.pyx \
	libpl/pdcode/homfly.pyx \
	libpl/pdcode/isomorphism.pyx \
	libpl/pdcode/plctopology.pyx

CY_FILES = $(PYX_FILES:.pyx=.c)
CY_H_FILES = \
	libpl/pdcode/diagram.h \
	libpl/pdcode/diagram_api.h

PXD_FILES = \
	libpl/plcurve.pxd \
	libpl/plCurve_h.pxd \
	libpl/tsmcmc.pxd \
	libpl/gsl_rng.pxd \
	libpl/pdcode/__init__.pxd \
	libpl/pdcode/diagram.pxd \
	libpl/pdcode/components.pxd \
	libpl/pdcode/homfly.pxd \
	libpl/pdcode/plctopology.pxd \
	libpl/pdcode/pdisomorphism.pxd

CPY_FILES = $(PY_FILES) $(PYX_FILES) $(PXD_FILES)
SETUP_PY = $(srcdir)/setup.py

EXTRA_DIST = $(SETUP_PY) tests $(CPY_FILES) $(CY_FILES) $(CY_H_FILES)

if BUILD_CYTHON
BUILT_SOURCES = $(CY_FILES) $(CY_H_FILES)
endif

if INSTALL_PYTHON_USER
PY_INSTALL_FLAGS = --user
else
if INSTALL_PYTHON_PREFIX
PY_INSTALL_FLAGS = --prefix=$(DESTDIR)/$(prefix)
endif
endif

all-local:
	CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)" \
		$(PYTHON) $(SETUP_PY) build_ext \
			-I$(top_srcdir)/src:$(top_srcdir)/planarmap \
			-L"$(abs_top_builddir)"/src/.libs:"$(abs_top_builddir)"/planarmap/.libs

PYLIBS = tsmcmc.so pdcode.so plcurve.so

clean-local:
	$(PYTHON) $(SETUP_PY) clean -a
	rm -r built_doc 2> /dev/null || true
	rm -r $(builddir)/libpl.egg-info 2> /dev/null || true

install-exec-local:
	$(PYTHON) $(SETUP_PY) install $(PY_INSTALL_FLAGS) \
		--single-version-externally-managed \
		--verbose \
		--record $(builddir)/install_files.txt

uninstall-local:
	cat $(builddir)/install_files.txt | xargs rm -rf
	rm $(builddir)/install_files.txt 2> /dev/null || true

#check-local:
#	$(PYTHON) $(SETUP_PY) test

html-local: install
	pushd $(srcdir)/doc ; make html BUILDDIR="$(abs_builddir)"/built_doc ; popd

if BUILD_CYTHON
libpl/pdcode/diagram_api.h: $(srcdir)/libpl/pdcode/diagram.c
libpl/pdcode/diagram.h: $(srcdir)/libpl/pdcode/diagram.c
.pyx.c:
	cython $<
endif

endif
