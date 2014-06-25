all-local:
	mkdir build 2> /dev/null || true
	pushd $(srcdir) ; \
	$(PYTHON) setup.py build_ext -L$(abs_top_builddir)/src/.libs -b$(abs_builddir) -t$(abs_builddir)/build ; \
	popd
	cp $(srcdir)/libpl/__init__.py $(builddir)/libpl

PYLIBS = tsmcmc.so pdcode.so plcurve.so

clean-local:
	pushd $(srcdir) ; \
	$(PYTHON) setup.py clean -b$(abs_builddir) -t$(abs_builddir)/build --build-lib=$(abs_builddir) ; \
	popd
	pushd $(builddir)/libpl ; \
	rm $(PYLIBS) 2> /dev/null || true ; \
	popd

install-exec-local:
	pushd $(srcdir) ; \
	$(PYTHON) setup.py install -b$(abs_builddir) --prefix=$(prefix) ; \
	popd
