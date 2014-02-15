#!/usr/bin/env bash

function cleanup() {
rm -fr $TMP_FILE
exit $1
}


echo "ffpcol: Comparing stdin vs file arg input" 2>&1 
TMP_FILE=$(mktemp)
BIN=../src
$BIN/ffpry -l 3 *.fna > $TMP_FILE
diff <( cat $TMP_FILE | $BIN/ffpcol ) <( $BIN/ffpcol $TMP_FILE ) &> /dev/null || cleanup 1

cleanup 0
