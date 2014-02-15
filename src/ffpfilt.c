/*****************************************************
* This code is distributed under a Non-commercial use 
* license.  For details see LICENSE.  Use of this
* code must be properly attributed to its author
* Gregory E. Sims provided that its use or derivative 
* use is non-commercial in nature.  Proper attribution        
* can be made by citing:
*
* Sims GE, et al (2009) Alignment-free genome 
* comparison with feature frequency profiles (FFP) and 
* optimal resolutions. Proc. Natl. Acad. Sci. USA.
* 106, 2677-82.
*
* Gregory E. Sims (C) 2010-2012
*
*****************************************************/
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <sys/stat.h>
#include "hash.h"
#include "utils.h"
#include "vstring.h"
#include "string.h"
#include "cdfmacros.h"
#include "sighandle.h"
#include "../config.h"

char PROG_NAME[FILENAME_MAX];
#define MAX_WORD_LENGTH 40 /**< Maximum feature length */
#define DEFAULT_PRECISION 2

void hashCol(FILE * fp);
int isKeyBased(FILE * fp);
int getKeyLength(FILE * fp);


char usage_str[] = "Usage: %s [OPTIONS]... [FILE]...\n\
This program filters FFP profiles based on numeric properties\n\n\
\t-s, --stats\tPrint stats to stderr\n\
\t-z, --zeros\tInclude zeros in statistics calculations\n\
\t-u FLOAT, --upper=FLOAT\tupper threshold limit\n\
\t-l FLOAT, --lower=FLOAT\tlower threshold limit\n\
\t-f, --freq\tExclude raw frequencies outside of -u and -l limits\n\
\t-n, --norm\tUse Normal CDF to exclude features out -u and -l probabilities\n\
\t-e, --evd\tUse Extreme value CDF to exclude features out -u and -l probabilities\n\
\t-d, --disable\tDisable classing of AAs and NAs\n\
\t-a, --amino\tInput is an amino acid sequence\n\
\t-t, --text\tInput is a text file\n\
\t-k, --key\tForce key output.\n\
\t-v, --version\tPrint version information\n\
\t-h, --help\tThis text.\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";


char *weightVector;

int Length = MAX_WORD_LENGTH;
			    /**< Feature length */
/**< It might be clearer if a distribution mode enum were used */

bool sflag = false;
bool zflag = false;
bool fflag = false;
bool nflag = false;
bool dflag = false;
bool aflag = false;
bool cflag = false;
bool tflag = false;
bool eflag = false;
bool lflag = false;
bool uflag = false;
bool kflag = false;

int pvalue = DEFAULT_PRECISION;
float upper = INT_MAX;
float lower = INT_MIN;
int main(int argc, char **argv)
{
    FILE *fp = NULL;
    int opt;
    int mode = nucleotide;
    int option_index = 0;

    static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"print-stats", no_argument, 0, 's'},
	{"version", no_argument, 0, 'v'},
	{"precision", required_argument, 0, 'p'},
	{"upper", required_argument, 0, 'u'},
	{"lower", required_argument, 0, 'l'},
	{"freq", no_argument, 0, 'f'},
	{"normal", no_argument, 0, 'n'},
	{"evd", no_argument, 0, 'e'},
	{"zeros", no_argument, 0, 'z'},
	{"disable", no_argument, 0, 'd'},
	{"keys", no_argument, 0, 'k'},
	{"amino", no_argument, 0, 'a'},
	{"text", no_argument, 0, 't'},
	{0, 0, 0, 0}
    };

    initSignalHandlers();

    strcpy(PROG_NAME,basename( argv[0] ));

    while ((opt = getopt_long(argc, argv, "svfnu:dhl:etak",
			      long_options, &option_index)) != -1)
	switch (opt) {
	case 's':		// print stats to stderr
	    sflag = !sflag;
	    break;
	case 'z':		// include zeros in stat calc
	    zflag = !zflag;
	    break;
	case 'l':
	    lflag = !lflag;
	    lower = atof(optarg);
	    break;
	case 'u':
	    uflag = !uflag;
	    upper = atof(optarg);
	    break;
	case 'e':		// evd filtering two args lower upper
	    eflag = !eflag;
	    break;
	case 'n':		// norm filtering two args lower upper
	    nflag = !nflag;
	    break;
	case 'f':		// freq filtering two args lower upper 
	    fflag = !fflag;
	    break;
	case 't':
	    tflag = !tflag;
	case 'a':
	    aflag = !aflag;
	    break;
	case 'd':
	    dflag = !dflag;
	    break;
	case 'k':
	    kflag = !kflag;
	    break;
	case 'v':
	    printVersion();
	    exit(EXIT_SUCCESS);
	    break;
	case 'h':
	    printUsageStr();
	    exit(EXIT_SUCCESS);
	    break;
	default:
	    printErrorUsageStr();
	    break;
	}

    if ((argc - optind) > 1)
	fatal_msg("Specify only one File argument\n");

    if (eflag + nflag + fflag > 1)
	fatal_msg("Specify -f, -n or -e, but not more than 1\n");


    if (aflag)
	mode = amino;

    if (tflag)
	mode = text;


    initHash(0, mode, !dflag);
    struct stat input;


    argv += optind;

    do {
	fp = stdin;
	if (*argv) {
	    if (!strcmp(*argv, "-"))
		fp = stdin;
	    else if ((fp = fopen(*argv, "r")) == NULL)
		fatal_msg("%s: %s.\n", *argv,strerror(errno));

	    if ( isDirectory(*argv) ) 
	    	fatal_msg("%s: %s\n",*argv, strerror(EISDIR));

	    argv++;
	} else if (isatty(STDIN_FILENO))
	    printErrorUsageStr();

	if (!isRegularFile(fp)) 
	    fp = convertPipeToFile(fp);

	if (!isKeyBased(fp))
	    fatal_msg("%s: Not a key valued FFP.", *argv);

	Length = getKeyLength(fp);
	hashCol(fp);

	if (fp != stdin)
	    fclose(fp);

    } while (*argv);


    return EXIT_SUCCESS;
}


/**
 * Convert a (key,value) FFP to a columnar FFP
 *
 * FFP in fp must be columnar FFP.
 *
 * @param fp A file pointer to a (key,value) FFP
 * @return none
 */



void hashCol(FILE * fp)
{
    unsigned val;
    char s[Length];
    char **keys;
    unsigned i;
    char c[2];
    int d;
    double x1 = 0;
    double x2 = 0;
    double sigma;
    long unsigned int ctr = 0;
    float alpha, beta;
    float zscore;
    int num;
    int row = 0;
    bool deleteKey;

    // Find all keys in the file.

    while (!feof(fp)) {
	d = fscanf(fp, "%s %u%*[\t]%[\n]", s, &val, c);
	s[Length] = '\0';
	//calc stats here.
	x1 += val;
	x2 += val * val;
	hashMax(s, val);	// Add to hash if greater than existing val
	ctr++;
	if (d == 3) {
	    row++;
	}
    }

    if (zflag)
	ctr = pow(2, Length) * row;	// if RY coded

    x1 /= ctr;
    x2 /= ctr;

    sigma = sqrt(x2 - x1 * x1);

    beta = sigma / PI_SQRT_6;
    alpha = x1 - beta * EULER_MASCHERONI;


    if (sflag)
	fprintf(stderr, "%.*e\t%.*e%.*e%.*e\n", pvalue, x1, pvalue, sigma,
		pvalue, alpha, pvalue, beta);


    hashKeys(&keys);

    num = numKeys();
    // delete any keys outside the range 
    for (i = 0; i < num; i++) {
	val = hashval(keys[i]);
	deleteKey = false;
	zscore = ((float) val - x1) / sigma;

	if (uflag) {
	    if (fflag && val > upper)
		deleteKey = true;
	    if (nflag && norm_cdf(zscore) > upper)
		deleteKey = true;
	    if (eflag && evd_cdf(alpha, beta, val) > upper)
		deleteKey = true;
	}
	if (lflag) {
	    if (fflag && val < lower)
		deleteKey = true;
	    if (nflag && norm_cdf(zscore) < lower)
		deleteKey = true;
	    if (eflag && evd_cdf(alpha, beta, val) < lower)
		deleteKey = true;
	}

	if (!deleteKey)
	    hashAssign(keys[i], 1);
	else
	    hashDel(keys[i]);

    }

    // grab updated set of keys
    hashKeys(&keys);
    rewind(fp);
    d = 0;
    while (d != EOF) {
	while ((d = fscanf(fp, "%s %u%[\r\n]", s, &val, c)) != EOF) {
	    s[Length] = '\0';
	    if (hashval(s)) {
		hashAssign(s, val);
	    }

	    if (d == 3)
		break;
	}
	if (d == 3) {
	    if (kflag) { 
	    	for (i = 0; i < numKeys()-1; i++) {
			printf("%s\t%u\t",keys[i], hashval(keys[i]));
			hashAssign(keys[i], 1);
	    	}
		printf("%s\t%u\n",keys[i], hashval(keys[i]));
		hashAssign(keys[i], 1);

	    } else {	    
	    	for (i = 0; i < numKeys()-1; i++) {
			printf("%u\t", hashval(keys[i]));
			hashAssign(keys[i], 1);
	    	}
	    printf("%u\n", hashval(keys[i]));
	    hashAssign(keys[i], 1);
	    }
	}
    }
}


// Add maximum word length given base function
