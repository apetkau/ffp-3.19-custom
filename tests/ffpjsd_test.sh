#!/usr/bin/env bash

if [ $(../src/ffpry -l 5 test{1..3}.fna | ../src/ffpcol |  ../src/ffprwn| ../src/ffpjsd | sum | cut -f1 -d" ") = 08044  ]
	then
	exit 0
else
	exit 1
	fi

