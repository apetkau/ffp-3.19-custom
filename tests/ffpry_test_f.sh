#!/usr/bin/env sh

src="../src"

echo "ffpry: Testing option -f, --feature-list" 2>&1
#Should look like:
#10      0       0       0
[ $($src/ffpry -l 3 -d -f fnalist.txt test1.fna | sum | cut -f1 -d" ") = 01325 ] || exit 1

exit 0

