#!/usr/bin/env bash

# Use checksum comparison to perform unit tests

echo "ffpaa: Testing -m,--multiple option " 2>&1 
# Test -l option 
# Output should look like:
#ASA	1	ASG	1	ASS	1	PIK	1	KAA	1	KAS	1	HIK	1	GKD	1	SSF	1	AAA	40	AAG	1	DSC	1	AGK	1	SFK	1	KDS	1	NPI	1	IKF	1	IKK	1	GHH	1	FSN	1	SAS	1	SGH	1	SNP	1	KFS	1	KKA	1	KKK	1	HHH	1	HHI	1	FKA	1	
#ASA	1	ASF	1	ASG	1	ASS	1	PIK	1	KAA	1	KAS	1	HIK	1	GKD	1	SSF	1	FPA	1	AAA	34	AAG	1	AAS	1	DSC	1	AGK	1	SFK	1	SFP	1	PAA	1	KDS	1	NPI	1	IKF	1	IKK	1	GHH	1	FSN	1	SAS	1	SGH	1	SNP	1	KFS	1	KKA	1	KKK	1	HHH	1	HHI	1	FKA	1	

[ $(../src/ffpaa -l 3 -m test{1,2}.faa | sum | cut -f1 -d" ") != 31357 ] && exit 1

exit 0

    1
