man_MANS = \
	ffpaa.1 \
	ffptree.1 \
	ffpboot.1 \
	ffpjsd.1  \
	ffpmerge.1  \
	ffpre.1  \
	ffpreprof.1  \
	ffprwn.1  \
	ffpry.1  \
	ffpvocab.1  \
	ffpvprof.1 \
	ffpcol.1 \
	ffptxt.1 \
	ffpfilt.1 \
	ffpcomplex.1 \
	ffpsubnam.1 \
	ffpdf.1

COPY=$$(date +'%Y')
DATE=$$(date +'%B %-d, %Y')
do_subst = sed -e 's,\[@\]VERSION\[@\],$(PACKAGE_VERSION),g' \
	       -e "s|\[@\]DATE\[@\]|$(DATE)|g" \
	       -e "s|\[@\]COPY\[@\]|$(COPY)|g" 

ffpaa.1 : ffpaa.1.in Makefile
	$(do_subst) < $(srcdir)/ffpaa.1.in > ffpaa.1


ffptree.1 : ffptree.1.in Makefile

	$(do_subst) < $(srcdir)/ffptree.1.in > ffptree.1

ffpboot.1 : ffpboot.1.in Makefile

	$(do_subst) < $(srcdir)/ffpboot.1.in > ffpboot.1


ffpjsd.1 : ffpjsd.1.in Makefile

	$(do_subst) < $(srcdir)/ffpjsd.1.in > ffpjsd.1


ffpmerge.1 : ffpmerge.1.in Makefile

	$(do_subst) < $(srcdir)/ffpmerge.1.in > ffpmerge.1


ffpre.1 : ffpre.1.in Makefile

	$(do_subst) < $(srcdir)/ffpre.1.in > ffpre.1

ffpreprof.1 : ffpreprof.1.in Makefile

	$(do_subst) < $(srcdir)/ffpreprof.1.in > ffpreprof.1


ffprwn.1 : ffprwn.1.in Makefile

	$(do_subst) < $(srcdir)/ffprwn.1.in > ffprwn.1

ffpry.1 : ffpry.1.in Makefile

	$(do_subst) < $(srcdir)/ffpry.1.in > ffpry.1

ffpvocab.1 : ffpvocab.1.in Makefile

	$(do_subst) < $(srcdir)/ffpvocab.1.in > ffpvocab.1
	
ffpvprof.1 : ffpvprof.1.in Makefile

	$(do_subst) < $(srcdir)/ffpvprof.1.in > ffpvprof.1

ffpcol.1 : ffpcol.1.in Makefile

	$(do_subst) < $(srcdir)/ffpcol.1.in > ffpcol.1

ffptxt.1 : ffptxt.1.in Makefile

	$(do_subst) < $(srcdir)/ffptxt.1.in > ffptxt.1


ffpfilt.1 : ffpfilt.1.in Makefile

	$(do_subst) < $(srcdir)/ffpfilt.1.in > ffpfilt.1

ffpcomplex.1 : ffpcomplex.1.in Makefile

	$(do_subst) < $(srcdir)/ffpcomplex.1.in > ffpcomplex.1

ffpdf.1 : ffpdf.1.in Makefile

	$(do_subst) < $(srcdir)/ffpdf.1.in > ffpdf.1

ffpsubnam.1 : ffpsubnam.1.in Makefile

	$(do_subst) < $(srcdir)/ffpsubnam.1.in > ffpsubnam.1

EXTRA_DIST = \
	ffpaa.1.in \
	ffptree.1.in \
	ffpboot.1.in \
	ffpjsd.1.in  \
	ffpmerge.1.in  \
	ffpre.1.in  \
	ffpreprof.1.in  \
	ffprwn.1.in  \
	ffpry.1.in  \
	ffpvocab.1.in  \
	ffpvprof.1.in \
	ffpcol.1.in \
	ffptxt.1.in \
	ffpfilt.1.in \
	ffpcomplex.1.in \
	ffpsubnam.1.in \
	ffpdf.1.in

CLEANFILES = \
	ffpaa.1 \
	ffptree.1 \
	ffpboot.1 \
	ffpjsd.1  \
	ffpmerge.1  \
	ffpre.1  \
	ffpreprof.1  \
	ffprwn.1  \
	ffpry.1  \
	ffpvocab.1  \
	ffpvprof.1 \
	ffpcol.1 \
	ffptxt.1 \
	ffpfilt.1 \
	ffpcomplex.1 \
	ffpsubnam.1 \
	ffpdf.1


