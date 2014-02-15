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
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include "utils.h"
#include "vstring.h"
#include "sighandle.h"
#include "../config.h"

char PROG_NAME[FILENAME_MAX];
#define COL_SIZE 1000  /**< Initial guess for number of columns in FFP */
#define ROW_SIZE 10    /**< Initial guess for number of rows in FFP */
#define DEFAULT_PRECISION 2 /**< Default precision for formated printing of normalized FFP */

void rownorm(FILE * fp);
void rownorml(FILE * fp);

char usage_str[] = "usage: %s [OPTION] ... [FILE] ...\n\
This program performs row normalization of an FFP vector file\n\n\
Given no options, the defeault behavior of the program is to\n\
generate row normalized relative frequency vectors\n\
\t-n, --largest-row\tNormalized by largest row sum\n\
\t-d=INT, --precision=INT\tSpecify n digits of decimal precision\n\
\t-v, --version\n\
\t-h, --help\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";


int precision = DEFAULT_PRECISION;
				 /**< Precision in Decimal places for normalized FFP */



int main(int argc, char **argv)
{
    FILE *fp;
    int opt;
    char nflag = 0;
    char dflag = 0;
    int dvalue = 0;
    int option_index = 0;

    static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"largest-row", no_argument, 0, 'n'},
	{"precision", no_argument, 0, 'd'},
	{"version", no_argument, 0, 'v'},
	{0, 0, 0, 0}
    };

    initSignalHandlers();

    strcpy(PROG_NAME,basename( argv[0] ));

    while ((opt = getopt_long(argc, argv, "nd:hv",
			      long_options, &option_index)) != -1)
	switch (opt) {
	case 'n':
	    nflag = 1;
	    break;
	case 'd':
	    dflag = 1;
	    dvalue = atoi(optarg);
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


    if (dflag)
	precision = dvalue;


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

    // check if not a seekable pipe
        if (!isRegularFile(fp)) 
	   fp = convertPipeToFile(fp);

	if (isKeyBased(fp))
	   fatal_msg("%s: Not a columnar FFP - see ffpcol.\n", *argv);


	if (nflag)
	    rownorml(fp);
	else
	    rownorm(fp);

	if (fp != stdin)
	    fclose(fp);

    } while (*argv);

    return EXIT_SUCCESS;
}



/**
 *
 * Performs standard row normalization
 *
 * Each row is converted into a relative
 * frequency vector.
 *
 * @param fp A file pointer 
 * @return void
 *
 */



/* Reconsider the design of this function -- It may be faster to do a double scan
 * of the file as in rownorml */

void rownorm(FILE * fp)
{
    long unsigned *sum;
    long unsigned val;
    int rowsize = ROW_SIZE;	/* initial guess */
    char s[3];
    int d=0;
    int row=0;

    sum = (long unsigned *) calloc(rowsize, sizeof(long unsigned));
 
    sum[row]=0;
    while ( (d = fscanf(fp, "%lu%[\n\r]", &val, s)) != EOF  ) {
	sum[row] += val;
        if (d == 2) { 
		row++;
		if (row == rowsize) {
		    rowsize += ROW_SIZE;
		    sum = (long unsigned *) realloc(sum, 
						  sizeof(long unsigned) *
						  rowsize);
		}
		sum[row]=0;
	}
    }


    row=0;
    rewind(fp);
    while ( (d = fscanf(fp, "%lu%[\n\r]", &val, s)) != EOF  ) {
	if (sum[row] == 0) {
		if (d==2) {
			printf("%.*e\n", precision, 0.0);
			warn_msg("Row %ld has row sum of 0.0\n", row + 1);
			row++;
		} else {
			printf("%.*e\t", precision, 0.0);
		}
	} else {
        	if (d == 2) { 
			printf("%.*e\n", precision, (float) val / sum[row]);
			row++;
		} else {
			printf("%.*e\t", precision, (float) val / sum[row]);
		}	
	
    	}
    }

    free(sum);
}




/**
 *
 * Performs row normalization by largest row
 *
 * Each row is converted into a relative
 * frequency vector.
 *
 * @param fp A file pointer
 * @return void
 *
 */



void rownorml(FILE * fp)
{
    long unsigned sum;
    long unsigned val;
    char s[3];
    int d;
    long unsigned maxsum = 0;

   sum = 0;
   while ( (d=fscanf(fp, "%lu%[\r\n]", &val,s)) != EOF ) {
	    sum += val;
	if ( d == 2 ) {
		if (sum > maxsum) {
		    maxsum = sum;
		}
		sum=0;
	}
    }

    rewind(fp);

   while ( (d=fscanf(fp, "%lu%[\r\n]", &val,s)) != EOF ) {
	    sum += val;
	    if ( d == 1) 
	    	printf("%.*e\t", precision, (float) val / maxsum);
	    else 
	    	printf("%.*e\n", precision, (float) val / maxsum);
	
	    
    }
}
