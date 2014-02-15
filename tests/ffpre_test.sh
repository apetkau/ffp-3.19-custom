#!/usr/bin/env sh

src="../src"
# Output should be:
#-0.431630

# have also gotten on linux AMD 32 bit
#-0.430877   chksum = 19136

if [ $($src/ffpre -l 3 test*.fna | sum | cut -f1 -d" ") = 19136 ]
	then
	exit 0
else
	exit 1
	fi



