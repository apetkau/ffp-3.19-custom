#!/usr/bin/env sh

src="../src"

echo "ffpry: Testing option -w, --mask" 2>&1
# Should look like:
#ATC     3       TCA     1
[ $($src/ffpry -l 3 -d -w 110 -q test3.fna  2> /dev/null | sum | cut -f1 -d" ") = 33236 ] || exit 1

exit 0

