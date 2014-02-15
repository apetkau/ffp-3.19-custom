#!/usr/bin/env sh

# Use checksum comparison to perform unit tests

echo "ffpaa: Testing -w,--mask  and -q,--quiet" 2>&1 
# Test -w option and option -q
#Output should look like:
#SNPI	1	AGKD	1	SGHH	1	KAAA	40	KASG	2	AAGK	1	SASS	1	NPIK	1	GHHH	1	HHHI	1	HHIK	1	ASSF	1	ASGH	1	SSFK	1	FSNP	1	ASAS	1	HIKF	1	PIKK	1	KDSC	1	KKKA	2	IKKK	1	FKAS	3	SFKA	1	KFSN	1	

[ $(../src/ffpaa -l 4 -w "0101" -q test1.faa | sum | cut -f1 -d" ") != 53937 ] && exit 1
exit 0

