#!/usr/bin/env sh

src="../src"

echo "ffpry: Testing option -m, --multiple" 2>&1
# Output should be:
#YRYR    4       RYRY    5
#RRRR    19      RRRY    1
[ $( $src/ffpry -l 4 -m  test2.fna| sum | cut -f1 -d" ") = 21944 ] || exit 1

exit 0



