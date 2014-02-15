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
#include <string.h>
#include <errno.h>
#include "utils.h"
#include "vstring.h"
#include "sighandle.h"
#include "../config.h"

char PROG_NAME[FILENAME_MAX];
#define VECTOR_SIZE 1000 /**< Initial guess for the number of columns in the FFP */
#define DEFAULT_THRESH 2 /**< The default frequency threshold for counting vocab features */

float vocab(FILE * fp, int threshold);

char usage_str[] = "Usage: %s [OPTION] ... [FILE] ...\n\
This program determines the number of features used in a vector\n\n\
Given no options, the defeault behavior of the program is to\n\
print out the number of features with frequencies greater than\n\
two in all rows.\n\
\t-f INT, --freq-thresh=INT\n\
\t-v, --version\n\
\t-h, --help\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";


int main(int argc, char **argv)
{
    FILE *fp;
    int opt;
    char fflag = 0;
    int threshold = DEFAULT_THRESH;
    int option_index = 0;

    static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"version", no_argument, 0, 'v'},
	{"freq-thresh", required_argument, 0, 'f'},
	{"amino", no_argument, 0, 'a'},
	{"text", no_argument, 0, 't'},
	{0, 0, 0, 0}
    };

    initSignalHandlers();

    strcpy(PROG_NAME,basename( argv[0] ));

    while ((opt = getopt_long(argc, argv, "hf:atv",
			      long_options, &option_index)) != -1)

	switch (opt) {
	case 'f':
	    fflag = 1;
	    threshold = atoi(optarg);
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

//add code to detect what kind of ffp is being given as input.


    if (threshold < 1)
	fatal_msg("Feature threshold must be greater than or equal to 1\n");


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

	if (!isKeyBased(fp))
	   fatal_msg("%s: Not a key valued FFP.\n", *argv);

	printf("%e\t\n", vocab(fp, threshold));

	if (fp != stdin)
	    fclose(fp);

    } while (*argv);
    return EXIT_SUCCESS;
}



/**
 *
 * Determine the number of features in the FFP which
 * occur more than threshold number of times
 *
 * Determine the feature usage (vocabulary usage) of
 * the FFP stored in the File pointed to by fp.
 *
 * @param fp A file pointer
 * @param threshold A frequency threshold
 * @return The number of features above the frequency threshold
 *
 */



float vocab(FILE * fp, int threshold)
{
    int rows = 0;
    rows = 0;
    char ch[2];
    unsigned numFeature = 0;
    unsigned val;
    unsigned lineno=0;
    unsigned line_offset=0;
    int items;

    errno=0;
    while ( (items=fscanf(fp, "%*s %u%[\r\n]", &val,ch)) != EOF ) { 
	//fprintf(stderr,"%s %d\n",strerror(errno),items);
	if (errno) 
		fatal_msg("Parse error at line %u char %ld: %s.\n",
			lineno,ftell(fp)-line_offset, strerror(errno) );

	if (val >= threshold)
		numFeature++;
	
	switch (items) {
	        case 2:
			line_offset=ftell(fp);
			lineno++;
	      		rows++;		
		case 1:
		break;	
		default:
			fatal_msg("Parse error at line %u char %ld.\n",
			lineno,ftell(fp)-line_offset);
		break;
	}
    }


    if (ferror(fp)) 
    	fatal_msg("Read Error: %s",strerror(errno));

    return ((float) numFeature / rows);
}
