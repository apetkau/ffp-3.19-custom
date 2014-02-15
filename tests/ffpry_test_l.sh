#!/usr/bin/env sh

src="../src"

echo "ffpry: Testing option -l, --length" 2>&1
#Output should be:
#RY      14      YR      11
[ $($src/ffpry -l 2 test1.fna | sum | cut -f1 -d" ") = 43576 ] || exit 1

exit 0




