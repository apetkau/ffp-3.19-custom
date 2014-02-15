#!/usr/bin/env bash

function cleanup() {
rm -fr $TMP_FILE
exit $1
}

TMP_FILE=$(mktemp)
BIN=../src
cat test*.fna > $TMP_FILE

echo "ffpry: Comparing stdin vs file arg input" 2>&1 
#process substitution
diff <( $BIN/ffpry <(cat test*.fna) ) <( cat test*.fna | $BIN/ffpry ) &> /dev/null || cleanup 1
#directly vi a file
diff <( $BIN/ffpry $TMP_FILE ) <( cat test*.fna | $BIN/ffpry ) &> /dev/null || cleanup 1


cleanup 0

