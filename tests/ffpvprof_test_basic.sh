#!/usr/bin/env sh
# output should look like this:

#3 4.000000e+00
#4 1.000000e+01
#5 1.500000e+01

echo "ffpvprof: Testing basic functionality" 2>&1
[ $( ../scripts/ffpvprof -s 3 -e 5 -p ../src test5.fna | sum | cut -f1 -d" ") = 18942 ] || exit 1

# Output should look like this:
#3 8.000000e+00
#4 1.400000e+01
#5 2.300000e+01

[ $( ../scripts/ffpvprof -s 3 -e 5 -r -p ../src test5.fna | sum | cut -f1 -d" ") =  19598  ] || exit 1

# Output should look like this
#3 1.000000e+00
#4 1.000000e+00
#5 1.000000e+00
[ $( ../scripts/ffpvprof -s 3 -e 5 -a -p ../src test1.faa | sum | cut -f1 -d" ") =  51196 ] || exit 1


exit 0

