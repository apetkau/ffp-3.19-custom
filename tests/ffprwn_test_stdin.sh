#!/usr/bin/env bash

function cleanup() {
rm -fr $TMP_FILE
exit $1
}


echo "ffprwn: Comparing stdin vs file arg input" 2>&1 
TMP_FILE=$(mktemp)
BIN=../src
$BIN/ffpry -l 3 *.fna | $BIN/ffpcol > $TMP_FILE
diff <( cat $TMP_FILE | $BIN/ffprwn ) <( $BIN/ffprwn $TMP_FILE ) &> /dev/null || cleanup 1

cleanup 0
