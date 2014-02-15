if WITH_GUI
   MAYBE_GUI = ffpgui 
endif

bin_SCRIPTS = ffpvprof ffpreprof ffpdf ffpsubnam $(MAYBE_GUI)

do_subst = sed -e 's,\[@\]VERSION\[@\],$(PACKAGE_VERSION),g' \
	       -e 's,\[@\]DOC\[@\],$(datadir)/doc/@PACKAGE@,g' \
	       -e 's,\[@\]EMAIL\[@\],$(PACKAGE_BUGREPORT),g' \
	       -e 's,\[@\]ISOSX\[@\],$(ISOSX),g'


if WITH_GUI
ffpgui : ffpgui.in Makefile

	$(do_subst) < $(srcdir)/ffpgui.in > ffpgui.tmp
	chmod +x ffpgui.tmp
	mv ffpgui.tmp ffpgui
endif

ffpdf : ffpdf.in Makefile

	$(do_subst) < $(srcdir)/ffpdf.in > ffpdf.tmp
	chmod +x ffpdf.tmp
	mv ffpdf.tmp ffpdf

ffpreprof : ffpreprof.in Makefile

	$(do_subst) < $(srcdir)/ffpreprof.in > ffpreprof.tmp
	chmod +x ffpreprof.tmp
	mv ffpreprof.tmp ffpreprof

ffpsubnam : ffpsubnam.in Makefile

	$(do_subst) < $(srcdir)/ffpsubnam.in > ffpsubnam.tmp
	chmod +x ffpsubnam.tmp
	mv ffpsubnam.tmp ffpsubnam

ffpvprof : ffpvprof.in Makefile

	$(do_subst) < $(srcdir)/ffpvprof.in > ffpvprof.tmp
	chmod +x ffpvprof.tmp
	mv ffpvprof.tmp ffpvprof



CLEANFILES=\
	  $(MAYBE_GUI) \
	  ffpdf \
	  ffpreprof \
	  ffpsubnam \
	  ffpvprof

EXTRA_DIST= \
	  ffpgui.in \
	  ffpdf.in \
	  ffpreprof.in \
	  ffpsubnam.in \
	  ffpvprof.in 

