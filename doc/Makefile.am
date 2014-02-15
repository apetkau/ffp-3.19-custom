docdir = $(datadir)/doc/@PACKAGE@
doc_DATA = manual.pdf tutorial.pdf

do_subst = sed -e 's,\[@\]VERSION\[@\],$(PACKAGE_VERSION),g' \
	       -e "s|\[@\]DATE\[@\]|$(DATE)|g" \
	       -e "s|\[@\]COPY\[@\]|$(COPY)|g" 

TUTORIAL_FILES = \
	tutorial.in

MANUAL_FILES = \
	../man/ffpry.1      \
	../man/ffpaa.1      \
	../man/ffptxt.1     \
	../man/ffpcol.1     \
	../man/ffprwn.1     \
	../man/ffpboot.1    \
	../man/ffpjsd.1     \
	../man/ffpdf.1      \
	../man/ffpre.1      \
	../man/ffpsubnam.1  \
	../man/ffpvocab.1   \
	../man/ffpcomplex.1 \
	../man/ffpfilt.1    \
	../man/ffpmerge.1   \
	../man/ffpreprof.1  \
	../man/ffptree.1    \
	../man/ffpvprof.1      


manual.pdf: $(MANUAL_FILES) Makefile

	groff -t -e -mandoc -Tps $(MANUAL_FILES) > manual.tmp
	ps2pdf14 manual.tmp manual.pdf
	rm -fr manual.tmp


tutorial.pdf: $(TUTORIAL_FILES) Makefile

	$(do_subst) < $(TUTORIAL_FILES) > tutorial.tmp1
	eqn -d'&&' tutorial.tmp1 | groff -t -ms -Tps  > tutorial.tmp2
	ps2pdf14 tutorial.tmp2 tutorial.pdf
	rm -fr tutorial.tmp*

CLEANFILES = \
	     manual.pdf \
	     tutorial.pdf

EXTRA_DIST = tutorial.in
