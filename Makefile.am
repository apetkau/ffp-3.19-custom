AUTOMAKE_OPTIONS = foreign
SUBDIRS = src scripts tests man doc examples

dist-hook: 
	cd admin && ./add_license.sh
	cd admin && ./fixDate.sh
	cd admin && ./fixversion.sh
