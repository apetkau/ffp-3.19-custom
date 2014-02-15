#!/usr/bin/env sh

src="../src"


# Output should be:
#YRYR    9       RYRY    11
echo "ffpry: Testing option -r, --disable-rev" 2>&1
[ $($src/ffpry -l 4 -r test1.fna | sum | cut -f1 -d" ") = 17342 ] || exit 1

exit 0


