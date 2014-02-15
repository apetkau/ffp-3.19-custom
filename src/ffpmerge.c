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
#include "hash.h"
#include "utils.h"
#include "vstring.h"
#include "sighandle.h"
#include "../config.h"

char PROG_NAME[FILENAME_MAX];
#define COL_SIZE 1000 /**< Intial guess of column number in the FP */
#define MAX_WORD_LENGTH 40 /**< Maximum feature length */

unsigned merge(FILE * fp, unsigned **vals, unsigned size);
void mergeHash(FILE * fp);
void mergeHashComplement(FILE * fp);
int isKeyBased(FILE * fp);
int getKeyLength(FILE * fp);


char usage_str[] = "Usage: %s [OPTIONS] ... [FILE]\n\
This program merges multiple rows of an FFP into a single row\n\n\
The behavior of the program is to merge all rows of the FFP.\n\
The FFP can be a columnar or key-valued FFP.\n\
\t-t, --text\tInput is text\n\
\t-a, --amino\tInput is amino acid\n\
\t-d, --disable\tDisable classing of input\n\
\t-k, --keys\tPrint key-value pairs.\n\
\t-c, --complement\tMerges keys with rev complement (must have -d).\n\
\t-v, --vesion\tDisplay version\n\
\t-h, --help\tThis messsage\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";


char *weightVector;
int Length = MAX_WORD_LENGTH;
			    /**< Feature Length */

int main(int argc, char **argv)
{
    FILE *fp;
    int opt;
    unsigned cols;
    unsigned i;
    unsigned *vals;
    char **keys;
    bool isKey = false;
    bool dflag = false;
    bool aflag = false;
    bool tflag = false;
    bool kflag = false;
    bool cflag = false;
    int option_index = 0;


    static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"disable", no_argument, 0, 'd'},
	{"keys", no_argument, 0, 'k'},
	{"amino", no_argument, 0, 'a'},
	{"text", no_argument, 0, 't'},
	{"complement", no_argument, 0, 'c'},
	{"version", no_argument, 0, 'v'},
	{0, 0, 0, 0}
    };

    initSignalHandlers();

    strcpy(PROG_NAME,basename( argv[0] ));

    while ((opt = getopt_long(argc, argv, "khdactv",
			      long_options, &option_index)) != -1)
	switch (opt) {
	case 'd':
	    dflag = !dflag;
	    break;
	case 'a':
	    aflag = !aflag;
	    break;
	case 't':
	    tflag = !tflag;
	    break;
	case 'v':
	    printVersion();
	    exit(EXIT_SUCCESS);
	    break;
	case 'k':
	    kflag = !kflag;
	    break;
	case 'c':
	    cflag = !cflag;
	    break;
	case 'h':
	    printUsageStr();
	    exit(EXIT_SUCCESS);
	    break;
	default:
	    printErrorUsageStr();
	    break;
	}

    if (aflag && dflag)
	fatal_msg("Error -a or -t, not both\n");

    if (cflag && !dflag)
	fatal_msg("Error -c must be used with -d");


    if (aflag)
	initHash(0, amino, !dflag);
    else if (tflag)
	initHash(0, text, !dflag);
    else
	initHash(0, nucleotide, !dflag);

// init to zero
    vals = (unsigned *) calloc(sizeof(unsigned), COL_SIZE);


// IT may be necessary to determine what kind of data is 
// being merged, add additional switches -disable-classes?

// masking should also be used in this function

    cols = COL_SIZE;


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

	isKey = isKeyBased(fp);

	if (!isKey)
	    cols = merge(fp, &vals, cols);
	else {
	    Length = getKeyLength(fp);
	    if (cflag)
	        mergeHashComplement(fp);
	    else
	        mergeHash(fp);
	}

	if (fp != stdin)
	    fclose(fp);

    } while (*argv);

    if (!isKey) {
	for (i = 0; i < cols-1; i++)
	    printf("%u\t", vals[i]);
	printf("%u\n", vals[i]);
    } else {
	hashKeys(&keys);
	if (kflag) {
	    for (i = 0; i < numKeys()-1; i++)
		printf("%s %u\t", keys[i], hashval(keys[i]));
	    printf("%s %u\n", keys[i], hashval(keys[i]));
	} else {
	    for (i = 0; i < numKeys()-1; i++)
		printf("%u\t", hashval(keys[i]));
	    printf("%u\n", hashval(keys[i]));
	}
	printf("\n");
	freeHash();
	//@todo Free keys?
    }

    free(vals);
    return EXIT_SUCCESS;
}

/**
 * Merges all the rows in a single column together
 *
 * All the frequencies in a column are added together
 * to the vector in vals.  The FFP in the file fp
 * must be stored in columnar format.  Space in vals
 * is dynamically reallocated as needed to accomodate
 * the size of the FFP.
 *
 * @param fp A file pointer to an FFP
 * @param vals A ptr to an array of unsigned integers
 * @param size An initial guess as to the size of vals.
 * @return The number of columns read.
 */


unsigned merge(FILE * fp, unsigned **vals, unsigned size)
{
    long unsigned i;
    int rows;
    unsigned val;
    unsigned cols = 0;
    char ch[2];
    int d;


    i = 0;
    rows = 0;
    while (!feof(fp)) {
	if (i < size) {
	    d = fscanf(fp, "%u%[\n\r]", &val, ch);
	    (*vals)[i++] += val;
	    if (d == 2) {
		cols = i;
		rows++;
		i = 0;
	    }
	} else {
	    size += COL_SIZE;
	    *vals = (unsigned *) realloc(*vals, sizeof(unsigned) * size);
	}
    }
    return cols;

}


/**
 * Merges all the rows in a (key,value) FFP
 *
 * All the frequencies in a (key,value) FFP row are added together
 * in a new hash.  The FFP in the file fp
 * must be a (Key,value) FFP.
 *
 * @param fp A file pointer to an FFP
 * @return None
 */


void mergeHash(FILE * fp)
{
    unsigned val;
    char ch[2];
    char s[Length];


    while (!feof(fp)) {
	fscanf(fp, "%s%[\t]%u%[\n\r]", s, ch, &val, ch);
	s[Length] = '\0';	/**<This line may be unnecessary**/
	hashAdd(s, val);	// increments hash by value;
    }
}

/**
 * Merges all the rows in a (key,value) FFP and orders by hash key
 *
 * All the frequencies in a (key,value) FFP row are added together
 * in a new hash.  The FFP in the file fp
 * must be a (Key,value) FFP.
 *
 * @param fp A file pointer to an FFP
 * @return None
 */


void mergeHashComplement(FILE * fp)
{
    unsigned val;
    char ch[2];
    char s[Length];
    char *r; //reverse complement

    while (!feof(fp)) {
	fscanf(fp, "%s%[\t]%u%[\n\r]", s, ch, &val, ch);
	s[Length] = '\0';	/**<This line may be unnecessary**/

	r = strdup(s, Length);
	rev(r,Length);
	complement(r,Length);

	// hash dna and rev complement, only store lower hash value for sorting
	int hashSVal = (*hashf)(s);
	int hashRVal = (*hashf)(r);

	if (hashSVal < hashRVal) {
    		hashAdd(s, val);	// increments hash by value;
	} else {
		hashAdd(r, val);
	}

	free(r);
    }
}
