#!/usr/bin/env bash
# Use checksum comparison to perform unit tests
echo "ffpaa: Testing basic functionality" 2>&1
# Output should look like:
#KKKA	1	HHIK	1	SFKA	1	PIKK	1	SASS	1	GHHH	1	AAAA	39	AAAG	1	AAGK	1	KASG	1	KKAA	1	IKFS	1	IKKK	1	FSNP	1	ASAS	1	ASGH	1	HIKF	1	SSFK	1	KAAA	1	KDSC	1	FKAS	1	SNPI	1	SGHH	1	HHHI	1	GKDS	1	AGKD	1	NPIK	1	KFSN	1	ASSF	1	
[ $(../src/ffpaa test1.faa | sum | cut -f1 -d" ") != 19101 ] && exit 1


# Output should look like:
#KKKA	1	HHIK	1	SFKA	1	PIKK	1	SASS	1	FPAA	1	GHHH	1	AAAA	32	AAAG	1	AAAS	1	AAGK	1	KASG	1	ASFP	1	KKAA	1	IKFS	1	IKKK	1	FSNP	1	ASAS	1	ASGH	1	HIKF	1	SFPA	1	SSFK	1	KAAA	1	KDSC	1	FKAS	1	AASF	1	SNPI	1	SGHH	1	HHHI	1	GKDS	1	AGKD	1	NPIK	1	KFSN	1	PAAA	1	ASSF	1	
[ $(../src/ffpaa test2.faa | sum | cut -f1 -d" ") != 35837 ] && exit 1

# Output should look like:
#KKKA	1	HHIK	1	SFKA	1	PIKK	1	SASS	1	GHHH	1	AAAA	39	AAAG	1	AAGK	1	KASG	1	KKAA	1	IKFS	1	IKKK	1	FSNP	1	ASAS	1	ASGH	1	HIKF	1	SSFK	1	KAAA	1	KDSC	1	FKAS	1	SNPI	1	SGHH	1	HHHI	1	GKDS	1	AGKD	1	NPIK	1	KFSN	1	ASSF	1	
#KKKA	1	HHIK	1	SFKA	1	PIKK	1	SASS	1	FPAA	1	GHHH	1	AAAA	32	AAAG	1	AAAS	1	AAGK	1	KASG	1	ASFP	1	KKAA	1	IKFS	1	IKKK	1	FSNP	1	ASAS	1	ASGH	1	HIKF	1	SFPA	1	SSFK	1	KAAA	1	KDSC	1	FKAS	1	AASF	1	SNPI	1	SGHH	1	HHHI	1	GKDS	1	AGKD	1	NPIK	1	KFSN	1	PAAA	1	ASSF	1	

[ $(../src/ffpaa test{1,2}.faa | sum | cut -f1 -d" ") != 00296 ] && exit 1



exit 0

