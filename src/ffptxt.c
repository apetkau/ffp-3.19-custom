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

/** @todo Need to implement stop word removal i.e. 'the' 'a' 'and' */

char PROG_NAME[FILENAME_MAX];
#define DEFAULT_WORD_LENGTH 6 /**< Default Feature length if not given by opt -l */

void parseFile(HASH *, FILE *);
void loopFeatureList(HASH * h, FILE * fp);
void loopRaw(HASH * h, FILE * fp);

char usage_str[] = "Usage: %s [OPTION] ... [FILE] ...\n\
This program generates an FFP vector from text data\n\n\
Given no options the default behavior of the program is to\n\
generate a Feature Frequency Profile (FFP) using features of\n\
length\n\
\t-l LEN, --length\n\
\t-f FILE, --feature-list\n\
\t-v, --version\n\
\t-h, --help\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";


int Length = DEFAULT_WORD_LENGTH;
				/**< Feature length to use if not specified by opt -l */
int maxWordSize = MAX_WORD_SIZE;
char *weightVector;



char *fvalue = NULL;  /**< -f File name to load features from */
char fflag = 0;	      /**< -f Option to only count features listed in a file */
char mflag = 0;       /**< unused but needed to link with hashroll.c */
bool zflag = false;
bool wflag = false;

int main(int argc, char **argv)
{
    int opt;
    HASH h;
    FILE *fp;

    int option_index = 0;

    static struct option long_options[] = {
	{"length", required_argument, 0, 'l'},
	{"feature-list", required_argument, 0, 'f'},
	{"help", no_argument, 0, 'h'},
	{"version", no_argument, 0, 'v'},
	{0, 0, 0, 0}
    };

    initSignalHandlers();

    strcpy(PROG_NAME,basename( argv[0] ));

    while ((opt = getopt_long(argc, argv, "l:f:hv",
			      long_options, &option_index)) != -1)

	switch (opt) {
	case 'l':
	    Length = atoi(optarg);
	    break;
	case 'f':
	    fflag = 1;
	    fvalue = optarg;
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

    if ( getenv("MAX_WORD_SIZE") )
	maxWordSize=atoi(getenv("MAX_WORD_SIZE"));	    

    if (Length > maxWordSize) 
	    fatal_msg("%d: Max Word size is : %d",Length,maxWordSize);


    init(&h, 0, text, 0, 0, Length);


    // If provided a feature list read it and store in hash

    if (fflag) 
	    parseFeatureList(&h,fvalue,Length,text);
    

// Must now process file arguments
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

	parseFile(&h, fp);

	if (fp != stdin)
	    fclose(fp);

    } while (*argv);

    freeHash(&h);

    return EXIT_SUCCESS;
}



void parseFile(HASH * h, FILE * fp)
{

// to eliminate a test w/n a tight loop
// First test here.
    if (fflag)
	loopFeatureList(h, fp);
    else
	loopRaw(h, fp);

}


// For finding only features w/n a feature list
void loopFeatureList(HASH * h, FILE * fp)
{
    struct stat fattr;
    static size_t optimal_size;
    ssize_t nr;
    static char *buf = NULL;

    //Determine optimal buffer size for file
    if (fstat(fileno(fp), &fattr))
	fatal_msg("fstat error.\n");

    //Find optimal size for Disk IO 
    optimal_size = (fattr.st_blksize >= BUFSIZ) ? fattr.st_blksize : BUFSIZ;

    buf = (char *) malloc(optimal_size);


    // replace w/ function pointers
    while ((nr = fread(buf, sizeof(char), optimal_size, fp)) != -1 && nr != 0) {

	    chkpushtxt(h, buf, nr);
    }

    printFeatures(h);
    freeHash(h);
}



// For finding all features

void loopRaw(HASH * h, FILE * fp)
{
    struct stat fattr;
    static size_t optimal_size;
    ssize_t nr;
    static char *buf = NULL;

    //Determine optimal buffer size for file
    if (fstat(fileno(fp), &fattr))
	fatal_msg("fstat error.\n");

    //Find optimal size for Disk IO 
    optimal_size = (fattr.st_blksize >= BUFSIZ) ? fattr.st_blksize : BUFSIZ;

    buf = (char *) malloc(optimal_size);


    // replace w/ function pointers
    while ((nr = fread(buf, sizeof(char), optimal_size, fp)) != -1 && nr != 0) {

	    pushtxt(h, buf, nr);
    }

    printFeatures(h);
    freeHash(h);
}
