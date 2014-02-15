#!/usr/bin/env bash
# Use checksum comparison to perform unit tests

echo "ffpaa: Testing option -d, --disable-classes" 2>&1 
# Test -d option
#Output should look like
#TNPI    1       GRDS    1       TGHH    1       ATGH    1       TWRA    1       AAAA    39      AAAG    1       HHIK    1       RATG    1       GHHH    1   KAAA     1       KQKA    1       QKAA    1       IKYT    1       AAGR    1       RDSC    1       PIKQ    1       TATT    1       HIKY    1       TTWR    1   KYTN     1       NPIK    1       IKQK    1       ATAT    1       YTNP    1       AGRD    1       HHHI    1       WRAT    1       ATTW    1

[ $(../src/ffpaa -l 4 -d test1.faa | sum | cut -f1 -d" ") != 25432 ] && exit 1

#Output should look like:
#TNPI    1       GRDS    1       YPAA    1       TGHH    1       ATGH    1       TWRA    1       AAAA    32      AAAG    1       AAAT    1       HHIK    1   RATG     1       AATY    1       GHHH    1       KAAA    1       KQKA    1       QKAA    1       IKYT    1       AAGR    1       RDSC    1       PIKQ    1   TATT     1       HIKY    1       ATYP    1       TTWR    1       KYTN    1       NPIK    1       IKQK    1       ATAT    1       YTNP    1       PAAA    1   TYPA     1       AGRD    1       HHHI    1       WRAT    1       ATTW    1
[ $(../src/ffpaa -l 4 -d test2.faa | sum | cut -f1 -d" ") !=  65022 ] && exit 1
#Output should look like:
#TNPI	1	GRDS	1	TGHH	1	ATGH	1	TWRA	1	AAAA	39	AAAG	1	HHIK	1	RATG	1	GHHH	1	KAAA	1	KQKA	1	QKAA	1	IKYT	1	AAGR	1	RDSC	1	PIKQ	1	TATT	1	HIKY	1	TTWR	1	KYTN	1	NPIK	1	IKQK	1	ATAT	1	YTNP	1	AGRD	1	HHHI	1	WRAT	1	ATTW	1	
#TNPI	1	GRDS	1	YPAA	1	TGHH	1	ATGH	1	TWRA	1	AAAA	32	AAAG	1	AAAT	1	HHIK	1	RATG	1	AATY	1	GHHH	1	KAAA	1	KQKA	1	QKAA	1	IKYT	1	AAGR	1	RDSC	1	PIKQ	1	TATT	1	HIKY	1	ATYP	1	TTWR	1	KYTN	1	NPIK	1	IKQK	1	ATAT	1	YTNP	1	PAAA	1	TYPA	1	AGRD	1	HHHI	1	WRAT	1	ATTW	1	
[ $(../src/ffpaa -l 4 -d test{1,2}.faa | sum | cut -f1 -d" ") != 24460  ] && exit 1

exit 0

