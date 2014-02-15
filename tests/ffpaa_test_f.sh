#!/usr/bin/env sh

# Use checksum comparison to perform unit tests

echo "ffpaa: Testing option -f, --feature-list" 2>&1 
# Test -f option

#output should be:
#1	34	1	1	0	

[ $(../src/ffpaa -l 9 -f faalist.txt test1.faa | sum | cut -f1 -d" ") != 18413 ] && exit 1

#In combination with -w
#1	1	34	1	
# Notice that b/c of the redundancy in features (after masking) that 
# one of the features is elimininated.

[ $(../src/ffpaa -l 9 -w 111000111 -f faalist.txt -q test1.faa | sum | cut -f1 -d" ") != 48943 ] && exit 1
exit 0

