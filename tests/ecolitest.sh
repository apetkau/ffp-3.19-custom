#!/usr/bin/env bash
set -x

# This script exhaustively tests every aspect of $SRC/ffpjsd
# including all distance measures, and different methods 
# of row normalization.

TEST_DIR=$(dirname $(readlink -f $0) )
SRC="${TEST_DIR%/*}/src"

function clean_up() {
	echo "Encountered unexpected output"
	rm -fr $TMPDIR;
	exit $1
}

# SETUP actions for early termination
trap "clean_up 1" ABRT HUP INT QUIT PIPE



function getFtpFiles() {
	local FTP=ftp://ftp.ncbi.nlm.nih.gov
	local DIR=genomes/Bacteria
	echo "Grabbing files..." 1>&2
	while read FILE ; do
		echo "$FILE" 1>&2
		wget "$FTP/$DIR/$FILE" --quiet
	done <<-EOF
	Escherichia_coli_536_uid58531/NC_008253.fna
	Escherichia_coli_55989_uid59383/NC_011748.fna
	Escherichia_coli_APEC_O1_uid58623/NC_008563.fna
	Escherichia_coli_ATCC_8739_uid58783/NC_010468.fna
	Escherichia_coli_BW2952_uid59391/NC_012759.fna
	Escherichia_coli_B_REL606_uid58803/NC_012967.fna
	Escherichia_coli_CFT073_uid57915/NC_004431.fna
	Escherichia_coli_E24377A_uid58395/NC_009801.fna
	Escherichia_coli_ED1a_uid59379/NC_011745.fna
	Escherichia_coli_HS_uid58393/NC_009800.fna
	Escherichia_coli_IAI1_uid59377/NC_011741.fna
	Escherichia_coli_IAI39_uid59381/NC_011750.fna
	Escherichia_coli_K_12_substr__DH10B_uid58979/NC_010473.fna
	Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna
	Escherichia_coli_O103_H2_12009_uid41013/NC_013353.fna
	Escherichia_coli_O111_H__11128_uid41023/NC_013364.fna
	Escherichia_coli_O127_H6_E2348_69_uid59343/NC_011601.fna
	Escherichia_coli_O157_H7_EC4115_uid59091/NC_011353.fna
	Escherichia_coli_O157_H7_EDL933_uid57831/NC_002655.fna
	Escherichia_coli_O157_H7_Sakai_uid57781/NC_002695.fna
	Escherichia_coli_O157_H7_TW14359_uid59235/NC_013008.fna
	Escherichia_coli_O26_H11_11368_uid41021/NC_013361.fna
	Escherichia_coli_O55_H7_CB9615_uid46655/NC_013941.fna
	Escherichia_coli_S88_uid62979/NC_011742.fna
	Escherichia_coli_SE11_uid59425/NC_011415.fna
	Escherichia_coli_SMS_3_5_uid58919/NC_010498.fna
	Escherichia_coli_UMN026_uid62981/NC_011751.fna
	Escherichia_coli_UTI89_uid58541/NC_007946.fna
	Escherichia_coli__BL21_Gold_DE3_pLysS_AG__uid59245/NC_012947.fna
	Escherichia_fergusonii_ATCC_35469_uid59375/NC_011740.fna
	Shigella_boydii_CDC_3083_94_uid58415/NC_010658.fna
	Shigella_boydii_Sb227_uid58215/NC_007613.fna
	Shigella_dysenteriae_Sd197_uid58213/NC_007606.fna
	Shigella_flexneri_2a_2457T_uid57991/NC_004741.fna
	Shigella_flexneri_2a_301_uid62907/NC_004337.fna
	Shigella_flexneri_5_8401_uid58583/NC_008258.fna
	Shigella_sonnei_Ss046_uid58217/NC_007384.fna
	EOF
}

TMPDIR=$(mktemp -d)
cd $TMPDIR
getFtpFiles

echo "Building FFP distance matrix and trees" 1>&2

\ls *.fna | cut -f1 -d"." > species.txt
$SRC/ffpry -l 17 *.fna | $SRC/ffpcol | tee ffp.col | $SRC/ffprwn | tee ffp.row | $SRC/ffpjsd -p species.txt | tee ffp.jsd | $SRC/ffptree > trees 2> progress

#Todo tee and save all trees -- then compare the trees produced with treedist

echo "Comparing JSD matrix to expected output"  1>&2 
[ $(sum ffp.jsd | cut -f1 -d" ") != 56601 ] && clean_up 1

echo "Comparing tree to expected output"  1>&2 
[ $(sum trees | cut -f1 -d" ") != 32573 ] && clean_up 1

echo "Comparing tree build progress to expected output"  1>&2 
[ $(sum progress | cut -f1 -d" ") != 19202 ] && clean_up 1

echo "Building Euclidean tree"
[ $($SRC/ffpjsd --euclid -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 53921 ] && clean_up 1

echo "Building Euclidean-squared tree"
[ $($SRC/ffpjsd --euclid2 -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 51531 ] && clean_up 1

echo "Building Euclidean-Norm=4 tree"
[ $($SRC/ffpjsd --euclid -n 4 -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 19301 ] && clean_up 1

echo "Building Cosine distance tree"
[ $($SRC/ffpjsd --cosine -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 27746 ] && clean_up 1

echo "Building Manhattan distance tree"
[ $($SRC/ffpjsd --manhattan -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 03569 ] && clean_up 1

echo "Building Pearson distance tree"
[ $($SRC/ffpjsd --pearson -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 20523 ] && clean_up 1

echo "*** Repeating with row normalized ffp"

echo "Building Euclidean tree"
[ $($SRC/ffpjsd --euclid -p species.txt ffp.row | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 52475 ] && clean_up 1

echo "Building Euclidean-squared tree"
[ $($SRC/ffpjsd --euclid2 -p species.txt ffp.row | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 49538 ] && clean_up 1

echo "Building Euclidean-Norm=4 tree"
[ $($SRC/ffpjsd --euclid -n 4 -p species.txt ffp.row | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 15520 ] && clean_up 1

echo "Building Cosine distance tree"
[ $($SRC/ffpjsd --cosine -p species.txt ffp.row | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 38883 ] && clean_up 1

echo "Building Manhattan distance tree"
[ $($SRC/ffpjsd --manhattan -p species.txt ffp.row | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 00574 ] && clean_up 1

echo "Building Pearson distance tree"
[ $($SRC/ffpjsd --pearson -p species.txt ffp.row | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 44745 ] && clean_up 1

echo "*** Repeating with normalization by largest row"

$SRC/ffprwn --largest-row ffp.col > ffp.rwn 

echo "Building JSD tree"
[ $($SRC/ffpjsd -p species.txt ffp.rwn | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 24892 ] && clean_up 1

echo "Building Euclidean tree"
[ $($SRC/ffpjsd --euclid -p species.txt ffp.rwn | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 32620 ] && clean_up 1

echo "Building Euclidean-squared tree"
[ $($SRC/ffpjsd --euclid2 -p species.txt ffp.rwn | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 01496 ] && clean_up 1

echo "Building Euclidean-Norm=4 tree"
[ $($SRC/ffpjsd --euclid -n 4 -p species.txt ffp.rwn | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 06552 ] && clean_up 1

echo "Building Cosine distance tree"
[ $($SRC/ffpjsd --cosine -p species.txt ffp.rwn | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 37382 ] && clean_up 1

echo "Building Manhattan distance tree"
[ $($SRC/ffpjsd --manhattan -p species.txt ffp.rwn | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 36128 ] && clean_up 1

echo "Building Pearson distance tree"
[ $($SRC/ffpjsd --pearson -p species.txt ffp.rwn | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 48476 ] && clean_up 1

echo "Building Chebyshev distance tree"
[ $($SRC/ffpjsd --chebyshev -p species.txt ffp.rwn | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 26851 ] && clean_up 1

echo "Building Canberra distance tree"
[ $($SRC/ffpjsd --canberra -p species.txt ffp.rwn | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 11965 ] && clean_up 1

echo "Building Evolutionary distance tree"
[ $($SRC/ffpjsd --evol -p species.txt ffp.rwn | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 00231 ] && clean_up 1

$SRC/ffpry -d -l 10 *.fna | $SRC/ffpcol -d > ffp.col

echo "Building Hamming distance tree"
[ $($SRC/ffpjsd --hamming -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 57877 ] && clean_up 1

echo "Building Matching distance tree"
[ $($SRC/ffpjsd --matching -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 06947 ] && clean_up 1


echo "Building Jaccard distance tree"
[ $($SRC/ffpjsd --jaccard -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 16975 ] && clean_up 1


echo "Building Tanimoto distance tree"
[ $($SRC/ffpjsd --tanimoto -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 62002 ] && clean_up 1


echo "Building Dice distance tree"
[ $($SRC/ffpjsd --dice -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 54897 ] && clean_up 1


echo "Building AntiDice distance tree"
[ $($SRC/ffpjsd --antidice -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 13413 ] && clean_up 1


echo "Building Sneath-Sokal distance tree"
[ $($SRC/ffpjsd --sneath -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 40524 ] && clean_up 1


echo "Building Hamman distance tree"
[ $($SRC/ffpjsd --hamman -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 14768 ] && clean_up 1


echo "Building Pearson Phi distance tree"
[ $($SRC/ffpjsd --phi -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 10865 ] && clean_up 1


echo "Building Anderberg distance tree"
[ $($SRC/ffpjsd --anderberg -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 03720 ] && clean_up 1


echo "Building Gower distance tree"
[ $($SRC/ffpjsd --gower -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 65378 ] && clean_up 1


echo "Building Russel-Rao distance tree"
[ $($SRC/ffpjsd --russel -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 16975 ] && clean_up 1


echo "Building Yule distance tree"
[ $($SRC/ffpjsd --yule -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 04408 ] && clean_up 1


echo "Building Ochiai distance tree"
[ $($SRC/ffpjsd --ochiai -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 07115 ] && clean_up 1


echo "Building Kulczynski distance tree"
[ $($SRC/ffpjsd --kulczynski -p species.txt ffp.col | $SRC/ffptree 2> /dev/null | sum | cut -f1 -d" ") != 33072 ] && clean_up 1

cd ..


clean_up 0


