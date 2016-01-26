if BUILD_PYTHON

PY_FILES = \
	__init__.py \
	data.py \
	graphs.py \
	pdcode/__init__.py

PYX_FILES = \
	plcurve.pyx \
	tsmcmc.pyx \
	pdstor.pyx \
	pdcode/diagram.pyx \
	pdcode/components.pyx \
	pdcode/homfly.pyx \
	pdcode/pdisomorphism.pyx

PXD_FILES = \
	plcurve.pxd \
	plCurve_h.pxd \
	tsmcmc.pxd \
	gsl_rng.pxd \
	pdcode/__init__.pxd \
	pdcode/diagram.pxd \
	pdcode/components.pxd \
	pdcode/homfly.pxd \
	pdcode/plctopology.pxd \
	pdcode/pdisomorphism.pxd

CPY_FILES = $(PY_FILES) $(PYX_FILES) $(PXD_FILES)
SETUP_PY = $(srcdir)/setup.py

EXTRA_DIST = libpl $(SETUP_PY) test

if INSTALL_PYTHON_USER
PY_INSTALL_FLAGS = --user
else
if INSTALL_PYTHON_PREFIX
PY_INSTALL_FLAGS = --prefix=$(prefix)
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

install-exec-local:
	$(PYTHON) $(SETUP_PY) install $(PY_INSTALL_FLAGS)

check-local:
	$(PYTHON) $(SETUP_PY) test

html-local: install
	pushd $(srcdir)/doc ; make html BUILDDIR="$(abs_builddir)"/built_doc ; popd

libpl:
	mkdir -p libpl
pydir:
	rsync --del -tapE $(srcdir)/libpl/ libpl/
pdcode:
	mkdir -p libpl/pdcode
$(CPY_FILES): %: $(srcdir)/libpl/%
	rsync -up $< libpl/$@
data:
	rsync -aup $(srcdir)/libpl/data libpl
util:
	rsync -aup $(srcdir)/libpl/util libpl

endif
