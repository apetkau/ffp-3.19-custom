#!/usr/bin/env bash
SUM=$(../src/ffpry -l 10 -r -d ecoli | sed 's/\([[:digit:]]\)\t/\1\n/g' | sort -k 1 | sum | cut -f1 -d" ")

[ "$SUM" != "15826" ] && exit 1

exit 0

# validated run should have checksum of 04150

