#!/usr/bin/env bash

#20	2	1	16	1	14

if [ $( ../src/ffpry -l 4 test{1,2,3}.fna | ../src/ffpcol |  ../src/ffpmerge | sum | cut -f1 -d" ") != 07986 ]	
	then
	exit 1
	fi

exit 0	




