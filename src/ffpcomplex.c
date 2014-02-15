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
#include <float.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include "hash.h"
#include "utils.h"
#include "vstring.h"
#include "cdfmacros.h"
#include "sighandle.h"
#include "../config.h"

char PROG_NAME[FILENAME_MAX];
#define MAX_WORD_LENGTH 40 /**< Maximum feature length */

void hashCol(FILE * fp);
int isKeyBased(FILE * fp);
int getKeyLength(FILE * fp);
float complexity(char *key);

char usage_str[] = "Usage: %s [OPTIONS]... [FILE]...\n\
This program filters words by complexity\n\n\
Given no options, the default behavior of the program is to\n\
assume Nucleic acids FFP input\n\
\t-a, --amino\tInput is amino acids\n\
\t-t, --text\tInput is text\n\
\t-d, --disable\tDisable classing of AAs and Nuc.\n\
\t-n, --normal\tAssume normal CDF distribution of complexity.\n\
\t-l, --lower\tLower limit for filtration, with -n this is a\n\
\t\t\tprobability.\n\
\t-u, --upper\tUpper limit for filtration, with -n this is a\n\
\t\t\tprobability.\n\
\t-s, --stats\tPrint complexity stats to stderr\n\
\t-h, --help\tThis text.\n\
\t-v, --version\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";

char *weightVector;
int Length = MAX_WORD_LENGTH;
			    /**< Feature length */

bool flagK = false;
bool flagA = false;
bool flagT = false;
bool flagD = false;
bool flagF = true;
bool flagU = false;
bool flagL = false;
bool flagN = false;
bool flagS = false;

float lower = -FLT_MAX;
float upper = FLT_MAX;

int main(int argc, char **argv)
{
    FILE *fp = NULL;
    struct stat input;
    int opt;

    int mode = nucleotide;
    int option_index = 0;

    static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"amino", no_argument, 0, 'a'},
	{"text", no_argument, 0, 't'},
	{"disable", no_argument, 0, 'd'},
	{"lower", required_argument, 0, 'l'},
	{"upper", required_argument, 0, 'u'},
	{"version", no_argument, 0, 'v'},
	{"normal", no_argument, 0, 'n'},
	{"stats", no_argument, 0, 's'},
	{0, 0, 0, 0}
    };


    initSignalHandlers();

    strcpy(PROG_NAME,basename( argv[0] ));

    while ((opt = getopt_long(argc, argv, "hatdvu:l:ns",
			      long_options, &option_index)) != -1)
	switch (opt) {
	case 'a':
	    flagA = !flagA;
	    break;
	case 't':
	    flagT = !flagT;
	    break;
	case 'd':
	    flagD = !flagD;
	    break;
	case 'n':
	    flagN = !flagN;
	    flagF = !flagF;
	    break;
	case 'l':
	    flagL = !flagL;
	    lower = atof(optarg);
	    break;
	case 'u':
	    flagU = !flagU;
	    upper = atof(optarg);
	    break;
	case 's':
	    flagS = !flagS;
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


    if (flagA && flagT)
	fatal_msg("Option -a or -t not both\n");


    if (flagA)
	mode = amino;
    else if (flagT)
	mode = text;

    initHash(0, mode, !flagD);

    if ((argc - optind) > 1)
	fatal_msg("Specify only one File argument\n");

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

    	// process if its a pipe

   	 if (!isRegularFile(fp)) {
		fp = convertPipeToFile(fp);
    	}
    	// gather complexity stats.

	if (!isKeyBased(fp))
	   fatal_msg("%s: Not a key valued FFP.\n", *argv);

	Length = getKeyLength(fp);
	hashCol(fp);

	fclose(fp);

    } while (*argv);

    return EXIT_SUCCESS;
}

/**
 * Calculate the complexity of a feature, with length n.
 * The frequencies are calcualted for all k-mers from
 * length 1 to n-1 which are sub-words of the feature. 
 * From this the log2 based Shannon Entropy is calculated.  
 * Lower values indicate lower complexity.
 *
 * @param key an amino acid, nucleotide or text feature 
 * @return float Value indicates the Shannon entropy 
 * @todo It would simplify the use of complexity for filtration if this function had its own
 *       Internal hash variable
 */


float complexity(char *key)
{
    int L, i, N;
    char s[Length];
    float entropy = 0;
    char **keys;
    N = 0;
    L = 1;

    while (L < Length) {
	for (i = 0; i < Length - L + 1; i++) {
	    strncpy(s, &key[i], L);
	    s[L] = '\0';
	    hashInc(s);
	}
	L++;
	N += L;
    }

    hashKeys(&keys);

    for (i = 0; i < numKeys(); i++)
	entropy -=
	    (float) hashval(keys[i]) / N * log2((float) hashval(keys[i]) / N);

    freeHash();

    return entropy;
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
    char c[2];
    bool deleteKey;
    float cval;
    float cvalMin=FLT_MAX;
    float cvalMax=0;
    double x1 = 0;
    double x2 = 0;
    double sigma;
    double zscore;
    unsigned ctr = 0;
    unsigned items;
    unsigned lineno=0;
    unsigned line_offset=0;

    errno=0;
    while ( (items=fscanf(fp, "%s %u%[\r\n]", s, &val,c)) != EOF ) { 
	//fprintf(stderr,"%s %d\n",strerror(errno),items);
	if (errno) 
		fatal_msg("Parse error at line %u char %ld: %s.\n",
			lineno,ftell(fp)-line_offset, strerror(errno) );
	
	switch (items) {
	        case 3:
			line_offset=ftell(fp);
			lineno++;	  
		case 2:
			cval = complexity(s);
			x1 += cval;
			x2 += cval * cval;
			ctr++;
		break;	
		default:
			fatal_msg("Parse error at line %u char %ld.\n",
			lineno,ftell(fp)-line_offset);
		break;
	}


	if (flagS) {
		if (cval < cvalMin) 
			cvalMin=cval;
		if (cval > cvalMax)
			cvalMax=cval;
	}

    }


    x1 /= ctr;
    x2 /= ctr;

    sigma = sqrt(x2 - x1 * x1);

    if (flagS){
	fprintf(stderr, "Avg\tstddev\tmin\tmax\n");
	fprintf(stderr, "%lf\t%lf\t%f\t%f\n", x1, sigma,cvalMin,cvalMax);
    }

    rewind(fp);

    // perform tests here to decide whether to keep ffp.

    errno=0;
    while ( (items=fscanf(fp, "%s %u%[\r\n]", s, &val,c)) != EOF ) { 
	//fprintf(stderr,"%s %d\n",strerror(errno),items);
	if (errno) 
		fatal_msg("Parse error at line %u char %ld: %s.\n",
			lineno,ftell(fp)-line_offset, strerror(errno) );
	
	    cval = complexity(s);
	    zscore = ((float) cval - x1) / sigma;

	    deleteKey = false;

	    if (flagU) {
		if (flagF && cval > upper)
		    deleteKey = true;
		if (flagN && norm_cdf(zscore) > upper)
		    deleteKey = true;
	    }

	    if (flagL) {
		if (flagF && cval < lower)
		    deleteKey = true;
		if (flagN && norm_cdf(zscore) < lower)
		    deleteKey = true;
	    }

	    if (!deleteKey)
		printf("%s %u", s, val);


	switch (items) {
	        case 3:
			printf("\n");
			line_offset=ftell(fp);
			lineno++;
	     	break;		
		case 2:
			printf("\t");
		break;	
		default:
			fatal_msg("Parse error at line %u char %ld.\n",
			lineno,ftell(fp)-line_offset);
		break;
	}
    }

}
