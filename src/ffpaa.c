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
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include "hashroll.h"
#include "mask.h"
#include "utils.h"
#include "vstring.h"
#include "sighandle.h"
#include "parse_features.h"
#include "../config.h"
#define DEFAULT_WORD_LENGTH 4 /**< Default Feature length if not given by opt -l */
#define CHAR_BUFFER_SIZE 10000 /**< Default buffer size for reading sequence files */


/** @Todo implement a variety of character classes 
 SE-B(14)	A, C, D, EQ, FY, G, H, IV, KR, LM, N, P, ST, W
 SE-B(10)	AST, C, DN, EQ, FY, G, HW, ILMV, KR, P
 SE-V(10)	AST, C, DEN, FY, G, H, ILMV, KQR, P, W
 Li-A(10)	AC, DE, FWY, G, HN, IV, KQR, LM, P, ST
 Li-B(10)	AST, C, DEQ, FWY, G, HN, IV, KR, LM, P
 Solis-D(10)	AM, C, DNS, EKQR, F, GP, HT, IV, LY, W
 Solis-G(10)	AEFIKLMQRVW, C, D, G, H, N, P, S, T, Y
 Murphy(10)	A, C, DENQ, FWY, G, H, ILMV, KR, P, ST
 SE-B(8)	AST, C, DHN, EKQR, FWY, G, ILMV, P
 SE-B(6)	AST, CP, DEHKNQR, FWY, G, ILMV
 Dayhoff(6)	AGPST, C, DENQ, FWY, HKR, ILMV

 These are from Edgar, R (2004, NAR, 32:380-385

 */


char PROG_NAME[FILENAME_MAX];

static void parseFile(HASH *, FILE *);
static void loopFeatureList(HASH * h, FILE * fp);
static void loopRaw(HASH * h, FILE * fp);

char *weightVector;  /**< A mask to allow mismatches in features */
int Length = DEFAULT_WORD_LENGTH;
				/**< Feature length to use if not specified by opt -l */
int Buffsize = CHAR_BUFFER_SIZE;/**< File input buffer, increase value for larger genomes. opt -b can be used to change this value */
int maxWordSize = MAX_WORD_SIZE;


bool wflag = false;	     /**<-w Supply a mismatch character mask */
char *wvalue = NULL; /**<-w ptr to mismatch character mask str */
bool zflag = false;	     /**<-z Create random mismatch character mask */
int zflagN = 0;	     /**<-z Maximum number of mismatches to allow */
bool sflag = false;	     /**<-s Fix random seed */
int sflagN = 0;	     /**<-s Random seed value */
bool qflag = false;	     /**<-q Force quiet operation */
bool fflag = false;	     /**<-f Restrict counting to a list of features*/
char *fvalue = NULL; /**<-f ptr to name of featmer mask file */
bool dflag = false;	     /**<-d disable classing of amino acids */
bool mflag = false;	     /**<-m option, FNA file contains multiple sequences */

char usage_str[] = "Usage: %s [OPTION] ... [FILE] ... \n\
This program generates an FFP vector of amino acid features\n\n\
Given no options the default behavior of the program is to\n\
generate a Feature Frequency Profile (FFP) using features of\n\
length 4.\n\
\t-l LEN, --length\n\
\t-f FILE, --feature-list\n\
\t-w STR, --mask\n\
\t-z INT, --rand-mask=INT\n\
\t-s INT, --rand-seed\n\
\t-q, --quiet\n\
\t-d, --disable-classes\n\
\t-m, --multiple\n\
\t-h, --help\n\
\t-v, --version\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";

int main(int argc, char **argv)
{
    int opt;
    FILE *fp;
    HASH h;

    int option_index = 0;

    static struct option long_options[] = {
	{"length", required_argument, 0, 'l'},
	{"disable-classes", no_argument, 0, 'd'},
	{"feature-list", required_argument, 0, 'f'},
	{"mask", required_argument, 0, 'w'},
	{"rand-mask", required_argument, 0, 'z'},
	{"rand-seed", required_argument, 0, 's'},
	{"quiet", no_argument, 0, 'q'},
	{"help", no_argument, 0, 'h'},
	{"multiple", no_argument, 0, 'm'},
	{"version", no_argument, 0, 'v'},
	{0, 0, 0, 0}
    };

  initSignalHandlers();

    while ((opt = getopt_long(argc, argv, "l:dw:z:s:qf:h?mv",
			      long_options, &option_index)) != -1)

	switch (opt) {
	case 'l':
	    Length = atoi(optarg);
	    break;
	case 'w':
	    wflag = !wflag;
	    wvalue = optarg;
	    break;
	case 'z':
	    zflag = !zflag;
	    zflagN = atoi(optarg);
	    break;
	case 's':
	    sflag = !sflag;
	    sflagN = atoi(optarg);
	    break;
	case 'q':
	    qflag = !qflag;
	    break;
	case 'd':
	    dflag = !dflag;
	    break;
	case 'f':
	    fflag = !fflag;
	    fvalue = optarg;
	    break;
	case 'm':
	    mflag = !mflag;
	    break;
	case 'h':
	    printUsageStr();
	    exit(EXIT_SUCCESS);
	    break;
	case 'v':
	    printVersion();
	    exit(EXIT_SUCCESS);
	    break;
	default:
	    printErrorUsageStr();
	    exit(EXIT_FAILURE);
	    break;
	}


    if ( getenv("MAX_WORD_SIZE") ) {
	maxWordSize=atoi(getenv("MAX_WORD_SIZE"));	    
    }	    

    if (Length > maxWordSize) 
	    fatal_msg("%d: Max Word size is : %d",Length,maxWordSize);

    if (sflag)
	srand((unsigned) sflagN);
    else
	srand((unsigned) time(NULL)*getpid());


    if (wflag) {
	weightVector = (char *) malloc(sizeof(char) * (Length + 1));
	strcpy(weightVector, wvalue);
	if (!qflag)
	    warn_msg("USING FEATURE MASK: %s\n", weightVector);

    }


    if (zflag) {
	while (zflagN > Length - 1)
	    zflagN--;
	weightVector = (char *) malloc(sizeof(char) * (Length + 1));
	randweight(weightVector, Length, zflagN);
	if (!qflag)
	    warn_msg("USING FEATURE MASK: %s\n", weightVector);
    }
    // Initialize the rolling hash
    // reverse is not applicable, therefore 0
    init(&h, (zflag || wflag), amino, !dflag, 0, Length);


    // If provided a feature list read it and store in hash

    if (fflag) 
	parseFeatureList(&h,fvalue,Length,amino);

// Must now process file arguments
    argv += optind;

    do {
	fp = stdin;
	if (*argv) {
	    if (!strcmp(*argv, "-"))
		fp = stdin;
	    else if ((fp = fopen(*argv, "r")) == NULL)
	    	fatal_msg("%s: %s\n", *argv, strerror(errno));

	    if ( isDirectory(*argv) ) 
	    	fatal_msg("%s: %s\n",*argv, strerror(EISDIR));

	    argv++;
	} else if (isatty(STDIN_FILENO))
	    printErrorUsageStr();

	parseFile(&h, fp);

	if (fp != stdin)
	    fclose(fp);

    } while (*argv);

    //freeHash(&h);

    return EXIT_SUCCESS;
}




/**
 *
 *  Parses FASTA format FNA file
 *
 *  This function can read from stdin or from files.
 *  The speed of the function is affected by the size
 *  of the input buffer.  Multiple FFPs are printed
 *  out with the -m option for each >header
 *  in the file. Another flag which affects this function
 *  is -f
 *
 *  @param h The hash table to use.
 *  @param fp The FNA file to parse
 *  @return none
 *
 */


static void parseFile(HASH * h, FILE * fp)
{

// to eliminate a test w/n a tight loop
// First test here.
    if (fflag)
	loopFeatureList(h, fp);
    else
	loopRaw(h, fp);

}


// For finding only features w/n a feature list
static void loopFeatureList(HASH * h, FILE * fp)
{
    struct stat fattr;
    static size_t optimal_size;
    ssize_t nr;
    char *buf;
    bool firstRecord=true;

    //Determine optimal buffer size for file
    if (fstat(fileno(fp), &fattr))
	fatal_msg("fd %s: File stat error.\n",fileno(fp));

    //Find optimal size for Disk IO 
    optimal_size = (fattr.st_blksize >= BUFSIZ) ? fattr.st_blksize : BUFSIZ;

    	if ((buf = (char *) malloc(optimal_size))==NULL) 
		fatal_at_line("%s\n",strerror(errno));

    // replace w/ function pointers
    while ((nr = fread(buf, sizeof(char), optimal_size, fp)) != -1 && nr != 0) {

	if (zflag || wflag)
	    chkpushaaw(h, buf, nr,firstRecord);
	else
	    chkpushaa(h, buf, nr,firstRecord);

	firstRecord=false;
    }

    printFeatures(h);
    free(buf);
    freeHash(h);
}



// For finding all features

static void loopRaw(HASH * h, FILE * fp)
{
    struct stat fattr;
    static size_t optimal_size;
    ssize_t nr;
    char *buf;
    bool firstRecord=true;


    //Determine optimal buffer size for file
    if (fstat(fileno(fp), &fattr))
	fatal_msg("fd %s: File stat error.\n",fileno(fp));

    //Find optimal size for Disk IO 
    optimal_size = (fattr.st_blksize >= BUFSIZ) ? fattr.st_blksize : BUFSIZ;

    	if ((buf = (char *) malloc(optimal_size))==NULL) 
		fatal_at_line("%s\n",strerror(errno));

    // replace w/ function pointers
    while ((nr = fread(buf, sizeof(char), optimal_size, fp)) != -1 && nr != 0) {

	if (zflag || wflag)
	    pushaaw(h, buf, nr,firstRecord);
	else
	    pushaa(h, buf, nr,firstRecord);

	firstRecord=false;
    }

    printFeatures(h);
    freeHash(h);
    free(buf);
}

