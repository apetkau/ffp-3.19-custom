#!/usr/bin/env bash

function cleanup() {
rm -fr $TMP_FILE
exit $1
}

echo "ffpaa: Comparing stdin vs file arg input" 2>&1 
TMP_FILE=$(mktemp)
BIN=../src
cat test*.faa > $TMP_FILE
diff <( $BIN/ffpaa <(cat test*.faa) ) <( cat test*.faa | $BIN/ffpaa ) &> /dev/null || cleanup 1
diff <( $BIN/ffpaa $TMP_FILE ) <( cat test*.faa | $BIN/ffpaa ) &> /dev/null || cleanup 1
cleanup 0

