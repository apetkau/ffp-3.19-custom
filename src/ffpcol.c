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
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "hash.h"
#include "utils.h"
#include "vstring.h"
#include "sighandle.h"
#include "../config.h"

#define MAX_WORD_LENGTH 40 /**< Maximum feature length */

char PROG_NAME[FILENAME_MAX];

void hashCol(FILE * fp);
int isKeyBased(FILE * fp);
int getKeyLength(FILE * fp);

char *weightVector;
int Length = MAX_WORD_LENGTH;
			    /**< Feature length */



char usage_str[] = "Usage: %s [OPTIONS]... [FILE]...\n\
This program creates columnar FFPs from key val output of ffpry,\n\
ffpaa, or ffptxt.\n\n\
Given no options, the defeault behavior of the program is to\n\
generate columnar FFPs for Nucleic acids.\n\
\t-a, --amino\tInput is Amino acid\n\
\t-t, --text\tInput is text\n\
\t-d, --disable\tDisable classing of AAs and Nuc.\n\
\t-V, --verbose\tBe more verbose.\n\
\t-h, --help\tThis text.\n\
\t-v, --version\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";


bool flagV = false;


int main(int argc, char **argv)
{
    FILE *fp = NULL;
    int opt;
    bool flagA = false;
    bool flagT = false;
    bool flagD = false;

    int mode = nucleotide;
    int option_index = 0;

    static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"amino", no_argument, 0, 'a'},
	{"text", no_argument, 0, 't'},
	{"disable", no_argument, 0, 'd'},
	{"version", no_argument, 0, 'v'},
	{"verbose", no_argument, 0, 'V'},
	{0, 0, 0, 0}
    };

    initSignalHandlers();
    
    strcpy( PROG_NAME, basename(argv[0]) );

    while ((opt = getopt_long(argc, argv, "hatdvV",
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
	case 'V':
	    flagV = !flagV;
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
    else {
    }

    initHash(0, mode, !flagD);


    if ((argc - optind) > 1)
	fatal_msg("Specify only one file argument.\n");

    argv += optind;

   do {
	fp = stdin;
	if (*argv) {
	    if (!strcmp(*argv, "-"))
		fp = stdin;
	    else if ((fp = fopen(*argv, "r")) == NULL){
		    fatal_msg("%s: %s.\n", *argv,strerror(errno));
		
	    }

	    if ( isDirectory(*argv) ) 
	    	fatal_msg("%s: %s\n",*argv, strerror(EISDIR));

	    argv++;
	} else if (isatty(STDIN_FILENO)) {
		printErrorUsageStr();
	}

	if (!isRegularFile(fp)) 
	    fp = convertPipeToFile(fp);

        if (!isKeyBased(fp))
	    fatal_msg("%s: Not a key valued FFP.\n", *argv);




	Length = getKeyLength(fp);
	hashCol(fp);

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
    //for accepting input via pipe      
    // Used in all cases
    unsigned val;
    char * s;
    char **keys;
    unsigned i;
    char c[3];
    int items;
    unsigned lineno=1;
    off_t line_offset=0;



    s=(char *) chkmalloc(sizeof(char),(Length+1));

//@TODO move this comment to description of convertFilePtr

// cut here
// If file is stdin we need to create a tmp file
// This is equivalent to to:
//  stdin | cat > /tmp/tmp.1a321asf 
// However we don't want a file passed to the
// program via  ffpcol < infile, to undergo this
// process, since a file opened in this manner
// is seekable and doesn't need to be copied.
// File opened this way ffpcol <(ffpry *.fna) 
// i.e. via process substitution is not seekable

    // Find all keys in the file.

    errno=0;

   

    //@TODO see if this can be sped up by using %*u to skip reading and conversion of val.
    while ( (items=fscanf(fp, "%s %u%[\r\n]", s, &val,c)) != EOF ) { 
	//if (errno) 
	//	fatal_msg("Parse error at line %u char %ld: %s. %d\n",
	//		lineno,ftell(fp)-line_offset, strerror(errno),errno );
	
	switch (items) {
	        case 3:
			if (flagV)
				fprintf(stderr,"Processed Line: %d\n",lineno);
			line_offset=ftell(fp);
			lineno++;	  
		case 2:
			hashInc(s);		// increments hash by value;
		break;	
		default:
			fatal_msg("Parse error at line %u char %ld.\n",
			lineno,ftell(fp)-line_offset);
		break;
	}
    }

    if (ferror(fp)) 
    	fatal_msg("Read Error: %s",strerror(errno));

    hashKeys(&keys);

    rewind(fp);

    fflush(fp); // Flush last character read


    //zero out the hash

    for (i = 0; i < numKeys(); i++) {
	hashAssign(keys[i], 0);
    }


	while ((items = fscanf(fp, "%s %u%[\r\n]", s, &val, c)) != EOF) {

	    switch (items) {
		    case 2:
	    		hashAssign(s, val);
		    break;
		    case 3:
	    		hashAssign(s, val);
	   		for (i = 0; i < numKeys()-1; i++) {
				printf("%u\t", hashval(keys[i]));
				hashAssign(keys[i], 0);
	   		}
	    		printf("%u\n",hashval(keys[i]));
	    		hashAssign(keys[i], 0);
		    break;
	    }
	}

  free(s);
  // No need to free keys, since program terminates
}


