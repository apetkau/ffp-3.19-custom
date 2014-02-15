#!/usr/bin/env bash

function cleanup() {
rm -fr $TMP_FILE
exit $1
}


echo "ffpvocab: Comparing to expected output." 2>&1 

# should look like: 
#1.400000e+00	
if [ $( ../src/ffpry -l 3 test*.fna | ../src/ffpvocab -f 3 | sum | cut -f1 -d" ") = 55055 ] ; then
	exit 0
else
	exit 1
fi

echo "ffpvocab: Comparing stdin vs file arg input" 2>&1 
TMP_FILE=$(mktemp)
BIN=../src
$BIN/ffpry *.fna > $TMP_FILE
diff <( cat $TMP_FILE | $BIN/ffpvocab ) <( $BIN/ffpvocab $TMP_FILE ) &> /dev/null || cleanup 1

cleanup 0

