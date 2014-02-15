#!/usr/bin/env bash
SUM=$(../src/ffpry -l 20 -r  ecoli.2 | sed 's/\([[:digit:]]\)\t/\1\n/g' | sort -k 1 | sum | cut -f1 -d" ")

[ "$SUM" != "58435" ] && exit 1

exit 0

# validated run should have checksum of 04150

