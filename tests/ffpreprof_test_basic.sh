#!/usr/bin/env sh

echo "ffpreprof: Testing basic functionality" 2>&1
if [ $( ../scripts/ffpreprof -s 3 -e 8 -p ../src test5.fna | sum | cut -f1 -d" ") =  61901 ]
	then
	exit 0
else
	exit 1
	fi
