if BUILD_PYTHON

PY_FILES = __init__.py data.py
PYX_FILES = plcurve.pyx pdcode.pyx tsmcmc.pyx pdstor.pyx
PXD_FILES = plcurve.pxd pdcode.pxd plCurve_h.pxd tsmcmc.pxd gsl_rng.pxd

CPY_FILES = $(PY_FILES) $(PYX_FILES) $(PXD_FILES)
SETUP_PY = $(srcdir)/setup.py

all-local: $(CPY_FILES) util
	$(PYTHON) $(SETUP_PY) build_ext -I$(top_srcdir)/src -L"$(abs_top_builddir)"/src/.libs

PYLIBS = tsmcmc.so pdcode.so plcurve.so

clean-local:
	$(PYTHON) $(SETUP_PY) clean -a
	pushd libpl ; \
	rm $(PYLIBS) 2> /dev/null || true ; \
	popd
	rm -r built_doc 2> /dev/null || true

install-exec-local: data
	$(PYTHON) $(SETUP_PY) install --prefix=$(prefix)

check-local:
	pushd $(srcdir)/test ; $(PYTHON) suite.py ; popd

html-local: install
	pushd $(srcdir)/doc ; make html BUILDDIR="$(abs_builddir)"/built_doc ; popd

libpl:
	mkdir -p libpl
$(CPY_FILES): %: $(srcdir)/libpl/%
	rsync -up $< libpl/$@
data:
	rsync -aup $(srcdir)/libpl/data libpl
util:
	rsync -aup $(srcdir)/libpl/util libpl

endif
