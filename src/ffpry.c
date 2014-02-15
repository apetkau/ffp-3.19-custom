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

/* Current revision: $Revision: 1.37 $ 
 * On Tag name:  $Name:  $ 
 * Date commited: $Date: 2012-02-26 01:57:45 $
 * Author: $Author: gsims $ 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include "hashroll.h"
#include "mask.h"
#include "ffpry.h"
#include "utils.h"
#include "vstring.h"
#include "sighandle.h"
#include "parse_features.h"
#include "../config.h"


#define DEFAULT_WORD_LENGTH 10  /**< Default Feature length if not given by opt -l */
#define CHAR_BUFFER_SIZE 10000	/**< Default buffer size for reading sequence files */


/* Function prototypes */

static void parseFile(HASH *, FILE *);
static void loopFeatureList(HASH * h, FILE * fp);
static void loopRaw(HASH * h, FILE * fp);


/* Global variables */

char PROG_NAME[FILENAME_MAX];
char *weightVector;               /**< A mask to allow mismatches in features */
int Length = DEFAULT_WORD_LENGTH; /**< Feature length to use if not specified by opt -l */
int maxWordSize = MAX_WORD_SIZE;  /**< Maximum allowed length for opt -l */
long Buffsize = CHAR_BUFFER_SIZE; /**< File input buffer, increase value for larger genomes. opt -b can be used to change this value */


/* Option flags and option arguments */
char *cvalue = NULL;    /**< Option argument string specifying character mask, should change to w */
char *fvalue = NULL;    /**< Option argument file name specifying a file of containing features */
bool zflag = false;	/**< Create a random character mask, -z */
int zflagN = 0;	        /**< Option argument to -z */
bool sflag = false;	/**< Specify a random seed, -s */
int sflagN = 0;	        /**< Option argument to -s */
bool qflag = false;	/**< Use quiet option */
bool wflag = false;	/**< Use a character mask, -w */
bool fflag = false;	/**< Use a list of features specified in a file, -f */
bool dflag = false;	/**< Disable RY coding, -d */
bool mflag = false;	/**< Multiple sequences in one file, -m */
bool rflag = true;	/**< Do reverse complement */


char usage_str[] = "Usage: %s [OPTIONS]... [FILE]... \n\
This program generates an FFP vector of RY-coded features\n\n\
Given no options the default behavior of the program is to\n\
generate a Feature Frequency Profile (FFP) using features of length\n\
10.\n\
\t-l LEN, --length=LEN\n\
\t-f FILE, --feature-list=FILE\n\
\t-w STR, --mask=STR\n\
\t-z K, --random-mask=K\n\
\t-s INT, --rand-seed=INT\n\
\t-q, --quiet\n\
\t-d, --disable\n\
\t-r, --disable-rev\n\
\t-m, --multiple\n\n\
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
	{"mask", required_argument, 0, 'w'},
	{"disable", no_argument, 0, 'd'},
	{"random-mask", required_argument, 0, 'z'},
	{"quiet", no_argument, 0, 'q'},
	{"random-seed", required_argument, 0, 's'},
	{"help", no_argument, 0, 'h'},
	{"feature-list", required_argument, 0, 'f'},
	{"multiple", no_argument, 0, 'm'},
	{"disable-rev", no_argument, 0, 'r'},
	{"version", no_argument, 0, 'v'},
	{0, 0, 0, 0}
    };


  initSignalHandlers();

  strcpy(PROG_NAME,basename( argv[0] ));

    while ((opt = getopt_long(argc, argv, "l:dw:z:s:qf:hmvr",
			      long_options, &option_index)) != -1)
	switch (opt) {
	case 'l':
	    Length = atoi(optarg);
	    break;
	case 'w':
	    if ((weightVector = (char *) malloc(sizeof(char) * strlen(optarg)))==NULL) 
		fatal_at_line("%s\n",strerror(errno));
	    strcpy(weightVector, optarg);
	    wflag = !wflag;
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
	case 'm':
	    mflag = !mflag;
	    break;
	case 'f':
	    fflag = !fflag;
	    fvalue = optarg;
	    break;
	case 'r':
	    rflag = !rflag;
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
	    exit(EXIT_FAILURE);
	    break;
	}

    if ( getenv("MAX_WORD_SIZE") ) 
	maxWordSize=atoi(getenv("MAX_WORD_SIZE"));	    

    if (Length > maxWordSize) 
	    fatal_msg("%d: Max Word size is: %d",Length,maxWordSize);

    if (sflag)
	srand((unsigned) sflagN);
    else
	srand((unsigned) time(NULL) * getpid());


    if (wflag) {
	if (strlen(weightVector) != Length)
	    fatal_msg("length of mask not equal to feature length.\n");
	if (!qflag)
	    warn_msg("Using feature mask: %s\n", weightVector);
    }



    if (zflag) {
	while (zflagN > Length - 1)
	    zflagN--;
	if ((weightVector = (char *) malloc(sizeof(char) * (Length + 1))) == NULL) 
		fatal_at_line("%s\n",strerror(errno));
	randweight(weightVector, Length, zflagN);
	if (!qflag)
	    warn_msg("Using feature mask: %s\n", weightVector);

    }

    //Initialize a hash
    init(&h, (zflag || wflag), nucleotide, !dflag, rflag, Length);

    //Fill with keys if restricting to a feature list
    if (fflag)
        parseFeatureList(&h,fvalue,Length,nucleotide);

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

    }
    while (*argv);

    //freeHash(&h);

    return EXIT_SUCCESS;
}


/**
 *
 *  Parses FASTA format FNA file
 *
 *  This function can read from stdin or from fasta files.
 *  The hash h, is populated with features differently
 *  if there is a feature list passed via the -f option.
 *
 *
 *  @param h The hash table to use.
 *  @param fp A file pointer to the nucleotide fasta file to parse.
 *  @return none
 *
 */


static void parseFile(HASH * h, FILE * fp)
{
    if (fflag)
	loopFeatureList(h, fp);
    else
	loopRaw(h, fp);

}


/**
 *
 *  Populates a hash with feature counts from a Nucleotide fasta FILE
 *  ptr restricting the search to features defined in the feature list
 *
 *  This function can read from fp which points to stdin or a fasta file.
 *  Function behaves differently if a spaced seed (or mask) is defined with
 *  -z or -w.
 *
 *
 *  @param h The hash table to use.
 *  @param fp A file pointer to the nucleotide fasta file to parse.
 *  @return none
 *
 */

static void loopFeatureList(HASH * h, FILE * fp)
{
    struct stat fattr;
    size_t optimal_size;
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
	    chkpushatgcw(h, buf, nr,firstRecord);
	else
	    chkpushatgc(h, buf, nr,firstRecord);

	firstRecord=false;
    }
    printFeatures(h);
    free(buf);
    freeHash(h);
}


/**
 *
 *  Populates a hash with feature counts from a Nucleotide fasta FILE
 *  ptr.
 *
 *  This function can read from fp which points to stdin or a fasta file.
 *  Function behaves differently if a spaced seed (or mask) is defined with
 *  -z or -w.
 *
 *  @param h The hash table to use.
 *  @param fp A file pointer to the nucleotide fasta file to parse.
 *  @return none
 *
 */

static void loopRaw(HASH * h, FILE * fp)
{
    struct stat fattr;
    size_t optimal_size;
    ssize_t nr;
    char *buf;
    bool firstRecord=true;
    

    //Determine optimal buffer size for file
    if (fstat(fileno(fp), &fattr))
	fatal_msg("fd %s : File stat error.\n",fileno(fp));

    //Find optimal size for Disk IO 
    optimal_size = (fattr.st_blksize >= BUFSIZ) ? fattr.st_blksize : BUFSIZ;

    if ((buf = (char *) malloc(optimal_size))==NULL) 
		fatal_at_line("%s\n",strerror(errno));

    // replace w/ function pointers
    while ((nr = fread(buf, sizeof(char), optimal_size, fp)) != -1 && nr != 0) {

	if (zflag || wflag)
	    pushatgcw(h, buf, nr,firstRecord);
	else
	    pushatgc(h, buf, nr,firstRecord);

	firstRecord=false;
    	
    }


    printFeatures(h);
    free(buf);
    freeHash(h);
}


