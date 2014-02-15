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
#include <ctype.h>
#include <time.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include "vstring.h"
#include "utils.h"
#include "sighandle.h"
#include "../config.h"

char PROG_NAME[] = "ffpboot";
#define COL_SIZE 1000 /**< Initial guess for columns in FFP */
#define DEFAULT_JACK 0.36787944117144232159 /**< Default Jackknife probability of deletion: 1/exp(1) */
#define INIT_ROWSIZE 10

void bootstrap(FILE * fp);
void jacknife(FILE * fp, float);

char usage_str[] = "Usage: %s [OPTIONS] ... [FILE] ...\n\
This program performs bootstrap permutations of the FFP vector\n\n\
Given no options, the defeault behavior of the program is to\n\
calculate a bootstrap permutation.\n\
\t-j, --jackknife\n\
\t-p PROB, --delete-prob=PROB\n\
\t-s INT, --rand-seed=INT\n\
\t-v, --version\n\
\t-h, --help\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";


int main(int argc, char **argv)
{
    FILE *fp;
    int opt;
    char jflag = 0;
    char sflag = 0;
    unsigned svalue = 0;
    float pvalue = DEFAULT_JACK;
    int option_index = 0;

    static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"delete-prob", required_argument, 0, 'p'},
	{"jackknife", no_argument, 0, 'j'},
	{"rand-seed", required_argument, 0, 's'},
	{"version", no_argument, 0, 's'},
	{0, 0, 0, 0}
    };

    initSignalHandlers();

    while ((opt = getopt_long(argc, argv, "jp:s:vh",
			      long_options, &option_index)) != -1)

	switch (opt) {
	case 'p':
	    pvalue = atof(optarg);
	    break;
	case 'j':
	    jflag = 1;
	    break;
	case 's':
	    sflag = 1;
	    svalue = atoi(optarg);
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

    if (sflag)
	srand((unsigned) svalue);
    else
	srand((unsigned) time(NULL) * getpid());


    if (pvalue > 1 || pvalue < 0) {
	fatal_msg("Jacknife deletion probability must be betwen 0 and 1\n");
    }



    
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

	if (!isRegularFile(fp)) 
	    fp = convertPipeToFile(fp);

	if (isKeyBased(fp))
	    fatal_msg("Input is not columnar format\n");

	if (jflag)
	    jacknife(fp, pvalue);
	else
	    bootstrap(fp);

	if (fp != stdin)
	    fclose(fp);

    } while (*argv);

    return EXIT_SUCCESS;
}



/**
 *
 * Performs a bootsrap permutation of a columnar FFP.
 *
 * Column numbers are sampled with replacement up to the
 * number of columns in the original matrix and
 * the chosen columns are printed to stdout.
 *
 * @param fp A file pointer to an FFP.
 * @return void
 *
 */



void bootstrap(FILE * fp)
{
    long unsigned i, j;
    unsigned *randCols;
    int rows;
    unsigned *vals;
    unsigned cols;
    unsigned size = COL_SIZE;	/* Inital guess */
    int rowsize = INIT_ROWSIZE;	/* initial guess */
    char s[2];
    int d;


    vals = (unsigned *) malloc(sizeof(unsigned) * size);

    i = 0;
    rows = 0;
    while (!feof(fp)) {
	if (i < size) {
	    d = fscanf(fp, "%u%[\n\r]", &vals[i++], s);
	    if (d == 2) {
		rows++;
		if (rows > rowsize) {
		    rowsize += INIT_ROWSIZE;
		}
	    }
	} else {
	    size += COL_SIZE;
	    vals = (unsigned *) realloc(vals, sizeof(unsigned) * size);
	}
    }

    cols = (unsigned) i / rows;

    //Now we know the number of columns

    randCols = (unsigned *) malloc(sizeof(unsigned) * cols);

    for (i = 0; i < cols; i++)
	randCols[i] = rand() % cols;

    for (i = 0; i < rows; i++) {
	for (j = 0; j < cols-1; j++)
	    printf("%u\t", vals[i * cols + randCols[j]]);
	printf("%u\n", vals[i * cols + randCols[j]]);
    }

    free(vals);
}



/**
 *
 * Performs a jackknife resampling of a columnar FFP.
 *
 * Column numbers are sampled without replacement using
 * a deletion jackknife. Each column has a fixed probability
 * of being deleted equal to the second argument del. 
 * The remainin columns are printed to stdout.
 *
 * @param fp A file pointer to an FFP.
 * @param del The deletion probability.
 * @return void
 *
 */

void jacknife(FILE * fp, float del)
{
    long unsigned i, j;
    char *randCols;
    int rows;
    unsigned *vals;
    unsigned cols;
    unsigned size = COL_SIZE;	/* Inital guess of column number */
    int rowsize = INIT_ROWSIZE;	
    char s[2];
    int d;


    vals = (unsigned *) malloc(sizeof(unsigned) * size);

    i = 0;
    rows = 0;
    while (!feof(fp)) {
	if (i < size) {
	    d = fscanf(fp, "%u%[\n\r]", &vals[i++], s);
	    if (d == 2) {

		rows++;
		if (rows > rowsize) {
		    rowsize += INIT_ROWSIZE;
		}
	    }
	} else {
	    size += COL_SIZE;
	    vals = (unsigned *) realloc(vals, sizeof(unsigned) * size);
	}
    }

    cols = (unsigned) i / rows;

    //Now we know the number of columns

    randCols = (char *) malloc(sizeof(char) * cols);

    for (i = 0; i < cols; i++) {
	randCols[i] = ((float) rand() / INT_MAX > del);
    }

    for (i = 0; i < rows; i++) {
	for (j = 0; j < cols; j++)
	    if (randCols[j])
		printf("%u\t", vals[i * cols + j]);
	printf("%u\n", vals[i * cols + j]);
    }

    free(vals);
}
