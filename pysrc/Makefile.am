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
	libpl/pdcode/homfly.pyx \
	libpl/pdcode/pdisomorphism.pyx

CY_FILES = $(PYX_FILES:.pyx=.c)

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

EXTRA_DIST = $(SETUP_PY) tests $(PY_FILES) $(PYX_FILES) $(PXD_FILES)
BUILT_SOURCES = $(CY_FILES)

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

libpl:
	mkdir -p libpl
pydir:
	rsync --del -tapE $(srcdir)/libpl/ libpl/
pdcode:
	mkdir -p libpl/pdcode
$(CPY_FILES):
	rsync -up $< $@
data:
	rsync -aup $(srcdir)/libpl/data libpl
util:
	rsync -aup $(srcdir)/libpl/util libpl

cy-stamp: $(CY_FILES)
	@rm -f cy-temp
	@touch cy-temp
	$(PYTHON) $(SETUP_PY) cythonize
	@mv -f cy-temp $@

$(CY_FILES): cy-stamp
## Recover from the removal of $@
	@dry=; for f in x $$MAKEFLAGS; do \
	  case $$f in \
	    *=*|--*);; \
	    *n*) dry=:;; \
	  esac; \
	done; \
	if test -f $@; then :; else \
	  $$dry trap 'rm -rf cy-lock cy-stamp' 1 2 13 15; \
	  if $$dry mkdir cy-lock 2>/dev/null; then \
## This code is being executed by the first process.
	    $$dry rm -f cy-stamp; \
	    $(MAKE) $(AM_MAKEFLAGS) cy-stamp; \
	    $$dry rmdir cy-lock; \
	  else \
## This code is being executed by the follower processes.
## Wait until the first process is done.
	    while test -d cy-lock && test -z "$$dry"; do \
	      sleep 1; \
	    done; \
## Succeed if and only if the first process succeeded.
	    $$dry test -f cy-stamp; exit $$?; \
	  fi; \
	fi

endif
