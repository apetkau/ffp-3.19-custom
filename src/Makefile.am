# what flags you want to pass to the C compiler & linker
#AM_CFLAGS = --pedantic -Wall -std=c99 -O3  -pg
AM_CPPFLAGS = --pedantic -Wall -std=c99 -O3  #-pg
AM_LDFLAGS = #-pg
#
# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = ffpry ffpaa ffprwn ffpjsd ffpboot ffpvocab ffpre ffpmerge ffpcol ffptxt ffpfilt ffpcomplex ffptree #ffpgui2
ffpry_SOURCES  = ffpry.c ffpry.h hashroll.c hashroll.h mask.c mask.h utils.c utils.h vstring.h sighandle.c sighandle.h parse_features.c parse_features.h 
ffpaa_SOURCES  = ffpaa.c hashroll.c hashroll.h mask.c mask.h utils.h utils.c vstring.h sighandle.c sighandle.h parse_features.h parse_features.c
ffprwn_SOURCES = ffprwn.c utils.c utils.h vstring.h sighandle.c sighandle.h
ffpjsd_SOURCES = ffpjsd.c utils.c utils.h vstring.h vstring.h sighandle.c sighandle.h
ffpboot_SOURCES = ffpboot.c utils.c utils.h vstring.h  sighandle.c sighandle.h
ffpvocab_SOURCES = ffpvocab.c vstring.h utils.c utils.h sighandle.c sighandle.h
ffpre_SOURCES = ffpre.c hashroll.c hashroll.h utils.c utils.h vstring.h sighandle.c sighandle.h
ffpmerge_SOURCES = ffpmerge.c hash.c hash.h utils.c utils.h vstring.h sighandle.c sighandle.h
ffpcol_SOURCES = ffpcol.c hash.c hash.h utils.c utils.h vstring.h sighandle.c sighandle.h
ffptxt_SOURCES = ffptxt.c hashroll.c hashroll.h utils.c utils.h vstring.h sighandle.c sighandle.h parse_features.c parse_features.h
ffpfilt_SOURCES = ffpfilt.c hash.c hash.h utils.c utils.h vstring.h cdfmacros.h sighandle.c sighandle.h
ffpcomplex_SOURCES = ffpcomplex.c hash.c hash.h utils.c utils.h vstring.h cdfmacros.h  sighandle.c sighandle.h
ffptree_SOURCES = ffptree.c  utils.c utils.h sighandle.c sighandle.h
#ffpgui2_SOURCES = tcl.c


# Binary specific libraries
# ffpgui2_LDADD = -ltk8.5 -ltcl8.5


# added this line otherwise received errors using 'make dist'
noinst_HEADERS = ffpry.h  hash.h mask.h parse_features.h utils.h codon.h vstring.h sighandle.h 

