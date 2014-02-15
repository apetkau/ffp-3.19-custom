#!/usr/bin/env sh

src="../src"

echo "ffpry: Testing option -d, --disable" 2>&1
# Output
#ATGC    10      CATG    4       TGCA    4

[ $($src/ffpry -l 4 -d test1.fna | sum | cut -f1 -d" ") = 60920 ] || exit 1

exit 0


