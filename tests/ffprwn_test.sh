#!/usr/bin/env bash

SRC="../src"

#Should look like:
#0.00e+00        0.00e+00        1.00e+00        0.00e+00        0.00e+00
#6.43e-01        7.14e-02        2.86e-01        0.00e+00        0.00e+00
#0.00e+00        0.00e+00        0.00e+00        5.00e-01        5.00e-01

echo "Testing ffprwn normalization" 2>&1

if [ $($SRC/ffpry -l 5 test{1,2,3}.fna | $SRC/ffpcol | $SRC/ffprwn   | sum | cut -f1 -d" ") = 25409 ]
	then
	exit 0
else
	exit 1
	fi

