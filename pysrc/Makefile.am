PY_FILES = __init__.py
PYX_FILES = plcurve.pyx pdcode.pyx tsmcmc.pyx
PXD_FILES = plcurve.pxd pdcode.pxd plCurve_h.pxd tsmcmc.pxd gsl_rng.pxd

CPY_FILES = $(PY_FILES) $(PYX_FILES) $(PXD_FILES)

all-local: setup.py $(CPY_FILES)
	$(PYTHON) setup.py build_ext -I$(top_srcdir)/src -L$(abs_top_builddir)/src/.libs

PYLIBS = tsmcmc.so pdcode.so plcurve.so

clean-local: setup.py
	$(PYTHON) setup.py clean -a
	pushd libpl ; \
	rm $(PYLIBS) 2> /dev/null || true ; \
	popd
	rm -r built_doc 2> /dev/null || true

install-exec-local:
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(prefix)

check-local:
	pushd test ; $(PYTHON) suite.py ; popd

html-local: install
	pushd $(srcdir)/doc ; make html BUILDDIR=$(abs_builddir)/built_doc ; popd

libpl:
	mkdir -p libpl
setup.py:
	rsync -up $(srcdir)/setup.py .
$(CPY_FILES): %: $(srcdir)/libpl/%
	rsync -up $< libpl/
