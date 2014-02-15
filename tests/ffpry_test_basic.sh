#!/usr/bin/env sh



src="../src"

echo "ffpry: Testing basic functionality." 2>&1

# Output should be:
#YRYRYRYRYR      3       RYRYRYRYRY      5
[ $($src/ffpry test1.fna | sum | cut -f1 -d" ") = 59338 ] || exit 1


exit 0
