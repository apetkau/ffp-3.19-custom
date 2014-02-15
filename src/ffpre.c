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
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include "hashroll.h"
#include "ffpry.h"
#include "utils.h"
#include "vstring.h"
#include "sighandle.h"
#include "../config.h"

char PROG_NAME[FILENAME_MAX];
#define DEFAULT_WORD_LENGTH 10
#define MAX_LENGTH 30
#define CHAR_BUFFER_SIZE 2000

/**< @todo this needs to be check thoroughly to make sure it produces the correct output
*/
void reNuc(FILE * fp);
void reAA(FILE * fp);
void reText(FILE * fp);
void loopRawRe(HASH * h, HASH * h1, HASH * h2, FILE * fp);

char usage_str[] = "Usage: %s [OPTION] ... [FILE]...\n\
This program prints out the relative entropy between an\n\
l-2 Markov model estimation of frequency and the observed\n\
frequency of features\n\n\
Given no options the default behavior of the program is to\n\
generate the relative entropy value (Kullbach-Leibler Divergence)\n\
of feature frequencies using l=10 and estimates of the frequencies\n\
using l=9 and l=8 to provide the estimate\n\
\t-l LEN, --length=LEN\n\
\t-d, --disable\n\
\t-a, --amino\n\
\t-t, --text\n\
\t-v, --version\n\
\t-v, --help\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";


char *weightVector;
int Length = DEFAULT_WORD_LENGTH;
int maxWordSize = MAX_WORD_SIZE;
int Buffsize = CHAR_BUFFER_SIZE;
unsigned int lcount; /**< Number of features counted of length l */
unsigned int l1count; /**< Number of features counted of length l-1 */
unsigned int l2count; /**< Number of featurs counted of length l-2 */

char aflag = 0;
char dflag = 0;
char tflag = 0;
char rflag = 1;
char fflag = 0;			// not used but needed to link with hashroll.c fix
char mflag = 0;			// not used currently, same as fflag

int main(int argc, char **argv)
{
    char **keys;
    int mode = nucleotide;
    int i;
    int opt;
    FILE *fp;
    double Ef;
    double N;
    double kld;
    HASH h, h1, h2;
    int option_index = 0;
    char *r;
    char *s;
    char *t;

    static struct option long_options[] = {
	{"length", required_argument, 0, 'l'},
	{"disable", no_argument, 0, 'd'},
	{"help", no_argument, 0, 'h'},
	{"amino", no_argument, 0, 'a'},
	{"text", no_argument, 0, 't'},
	{"no-reverse", no_argument, 0, 'r'},
	{"version", no_argument, 0, 'v'},
	{0, 0, 0, 0}
    };

    initSignalHandlers();

    strcpy(PROG_NAME,basename( argv[0] ));

    while ((opt = getopt_long(argc, argv, "l:dhatvr",
			      long_options, &option_index)) != -1)

	switch (opt) {
	case 'l':
	    Length = atoi(optarg);
	    break;
	case 'd':
	    dflag = 1;
	    break;
	case 'r':
	    rflag = !rflag;
	    break;
	case 'a':
	    aflag = 1;
	    mode = amino;
	    break;
	case 't':
	    tflag = 1;
	    mode = text;
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

    if ( getenv("MAX_WORD_SIZE") ) {
	maxWordSize=atoi(getenv("MAX_WORD_SIZE"));	    
    }

  if (Length > maxWordSize) 
        fatal_msg("%d: Max Word size is : %d",Length,maxWordSize);


// check -l length to make sure it is valid

    if (Length < 3)
	fatal_msg("Feature Length must be 3 or greater\n");

    init(&h, 0, mode, !dflag, rflag, Length);
    init(&h1, 0, mode, !dflag, rflag, Length - 1);
    init(&h2, 0, mode, !dflag, rflag, Length - 2);
// Must now process file arguments

    argv += optind;

    do {
	fp = stdin;
	if (*argv) {
	    if (!strcmp(*argv, "-"))
		fp = stdin;
	    else if ((fp = fopen(*argv, "r")) == NULL)
		fatal_msg("%s: No Such file", *argv);
	    argv++;
	} else if (isatty(STDIN_FILENO))
	    printErrorUsageStr();

	loopRawRe(&h, &h1, &h2, fp);

	if (fp != stdin)
	    fclose(fp);

    } while (*argv);


    r = (char *) malloc(sizeof(char) * Length);
    s = (char *) malloc(sizeof(char) * Length);
    t = (char *) malloc(sizeof(char) * Length);

//@todo fix to work with amino acids and text

    lcount = sumValues(&h);
    l1count = sumValues(&h1);
    l2count = sumValues(&h2);

    hashKeys(&h, &keys);

    N = (double) l2count / l1count / l1count;


    kld = 0.0;
    for (i = 0; i < h.keyN; i++) {
	r = strncpy(r, keys[i] + 1, Length - 1);
	r[Length - 1] = '\0';
	s = strncpy(s, keys[i], Length - 1);
	s[Length - 1] = '\0';
	t = strncpy(t, keys[i] + 1, Length - 1);
	t[Length - 2] = '\0';
	Ef = (double) hashValNuc(&h1, r) * hashValNuc(&h1, s) / hashValNuc(&h2,
									   t) *
	    N;
	if (isnormal(Ef)) {
	    kld -= (double) Ef *log2(hashValNuc(&h, keys[i]) / Ef / lcount);
	}
    }
    printf("%lf\n", kld);
    // printf("\n");
    // should free keys too;

    // (a+b)/An   



    freeHash(&h);
    freeHash(&h1);
    freeHash(&h2);
    return EXIT_SUCCESS;
}

// Could change this to an array of hashes
void loopRawRe(HASH * h, HASH * h1, HASH * h2, FILE * fp)
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


    // Make hash aware of its own mode amino, nucleotide or text
    // replace w/ function pointers
    while ((nr = fread(buf, sizeof(char), optimal_size, fp)) != -1 && nr != 0) {
	if (aflag) {
	    pushaa(h,  buf, nr,false);
	    pushaa(h1, buf, nr,false);
	    pushaa(h2, buf, nr,false);
	} else if (tflag) {
	    pushtxt(h,  buf, nr);
	    pushtxt(h1, buf, nr);
	    pushtxt(h2, buf, nr);
	} else {
	    pushatgc(h,  buf, nr,false);
	    pushatgc(h1, buf, nr,false);
	    pushatgc(h2, buf, nr,false);
	}


    }
}
