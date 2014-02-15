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
/* UTILS.C */

#define _POSIX_C_SOURCE  200112L  // To use fdopen
#define _XOPEN_SOURCE  600  // to use mkstemp
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <limits.h>
#include "utils.h"

extern char PROG_NAME[];

/**@todo Error functions placed in their own object file */

char error_usage_str[] = "Usage: %s [OPTION]... [FILE]... \n\
Try `%s --help' for more information\n";


void printErrorUsageStr()
{
    fprintf(stderr, error_usage_str, PROG_NAME, PROG_NAME);
    exit(EXIT_FAILURE);
}


void fatal_msg(char *fmt, ...)
{
    va_list args;
    fprintf(stderr, "%s: ", PROG_NAME);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    exit(EXIT_FAILURE);
}

void fatal_at_line(char *fmt, ...)
{
    va_list args;	
    fprintf(stderr,"%s: %s:%d in function \'%s\': ",
	       PROG_NAME,__FILE__,__LINE__,__func__);
    va_start(args, fmt);
    vfprintf(stderr,fmt,args);
    va_end(args);
    exit(EXIT_FAILURE);
}


void warn_msg(char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "%s: ", PROG_NAME);
    vfprintf(stderr, fmt, args);
    va_end(args);
}


#define DIR_SEPARATOR '/'
#define IS_DIR_SEPARATOR(ch) ( (ch) == DIR_SEPARATOR ) 


char * basename (const char *name)
{
	const char *base;

	for (base = name; *name; name++) 
		if (IS_DIR_SEPARATOR (*name)) 
			base = name + 1;
	
	return (char *) base;
}




bool isDirectory(char * fname) {
	struct stat input;
	stat(fname,&input);
	return (S_ISDIR( input.st_mode) );
}



/**@todo possibly place isRegularFile and convertPipeToFile in a special file utilities object file */

/**
 *
 * Test whether a file is a regular file or FIFO/pipe
 * 
 * fstats the descriptor number of the file pointed to
 * by fp and returns result.  Exits fatally if the fstat
 * doesn't succeed.
 *
 * @param fp File pointer to be tested.
 * @return True or False
 * @retval 1 Is a regular file
 * @retval 0 Is not a regular file - a pipe or FIFO
 * 
 */

int isRegularFile(FILE * fp)
{
    struct stat fattr;
    if (fstat(fileno(fp), &fattr))
	fatal_msg("fstat error.\n");
    return (S_ISREG(fattr.st_mode));
}


int dirExists(char * name)
{
    struct stat fattr;
    if ((stat(name, &fattr))==0)
    	return (S_ISDIR(fattr.st_mode));
    return 0;
}




/**
 *
 * Test whether a file is a regular file or FIFO/pipe
 * 
 * Streams the contents of the FIFO/PIPE to a temporary
 * file.  The temp file is usually saved in /tmp. It is
 * immediately unlinked after creation so you will not
 * be able to see it in the file system.  See tempfile()
 * for specifics.  The file pointer of the
 * temporary file is returned.  
 *
 * @param fp File pointer to be tested.
 * @return FILE*
 * @retval FILE* Pointer to temporary file structure.
 *
 *
 */

/**@TODO Could fork and run this so file could be parsed immediately */

char * tempfile=NULL;

FILE *convertPipeToFile(FILE * fp)
{
    FILE *tmp;
    struct stat fattr;
    static size_t optimal_size;
    ssize_t nr, nw;
    static char *buf = NULL;
    char tmpdir[FILENAME_MAX];
    int fd;

    strcpy(tmpdir,"/tmp");

    if ( getenv("TMP") ) 
	strcpy(tmpdir,getenv("TMP"));

    if ( ! dirExists(tmpdir) )
	fatal_msg("%s: Directory not found.\n",tmpdir);

    if ((tempfile=(char*)malloc(FILENAME_MAX*sizeof(char)) ) == NULL ) 
	fatal_at_line("%s\n",strerror(errno));
	
    sprintf(tempfile,"%s/ffp.XXXXXX",tmpdir);

    //For some reason tempfile is not rewindable.
    if ((fd=mkstemp(tempfile)) == -1 )
	fatal_msg("%s: %s\n",tempfile,strerror(errno));


    if ((tmp=fdopen(fd, "w+")) == NULL )
	fatal_msg("%s: %s\n",tempfile,strerror(errno));

    if (fstat(fd, &fattr))
	fatal_msg("%s: %s\n",tempfile,strerror(errno));

    //Find optimal size for Disk IO 
    optimal_size = (fattr.st_blksize >= BUFSIZ) ? fattr.st_blksize : BUFSIZ;

    if ((buf = (char *) malloc(optimal_size)) == NULL)
	fatal_at_line("%s\n",strerror(errno));

    while ((nr = fread(buf, sizeof(char), optimal_size, fp)) != -1 && nr != 0)
	if ((nw = fwrite(buf, sizeof(char), nr, tmp)) == 0 || nw == -1)
	    fatal_msg("%s: %s\n",tempfile,strerror(errno));

    free(buf);
    rewind(tmp);
    return tmp;
}


unsigned long numFeatures(int l)
{
    if (l % 2 == 1)
	return pow(2, l - 1);
    else
	return (pow(2, l) + pow(2, l / 2)) / 2;
}


/**
 *
 * Tests whether the string is a valid nucleotide sequence
 * 
 * Valid nucleotide base characters: AaTtGgCcYyRr.  The
 * paramter s is converted to upper case if any lowercase
 * characters are encountered.  Calls macro isBase.
 *
 * @param s a string representing a sequence
 * @param len the string length
 * @return True or False
 * @retval 1 Is a valid nucleotide sequence
 * @retval 0 Is an invalid nucleotide sequence
 * 
 */


char isValid(register const char *s, int len)
{
    int i = 0;
    while (i < len) {
	if (!isBase(s[i++]))
	    return 0;
    }
    return 1;
}



/**
 *
 * Tests whether the string is a valid Amino Acid sequence
 * 
 * Valid amino acid characters: AaCcDdEeFfGgHhIiJjKkLlMmNn
 * PpQqRrSsTTVvWxYyZz.  The paramter s is converted to 
 * upper case if any lowercase characters are encountered.  
 * Calls macro isBase.
 *
 * @param s a string representing a sequence
 * @param len the string length
 * @return True or False
 * @retval 1 Is a valid AA sequence
 * @retval 0 Is an invalid AA sequence
 * 
 */


char isValidAA(register const char *s, int len)
{
    int i = 0;
    while (i < len) {
	if (!isAA(s[i++]))
	    return 0;
    }
    return 1;
}


/**
 *  
 * Find the number of columns in a columnar FFP
 *
 * @param fp a file pointer to an FFP
 * @return The number of columns
 *
 */

unsigned int numCols(FILE * fp)
{

    int ch;
    unsigned int cols = 0;

    do {
	ch = fgetc(fp);
	if (isspace(ch)) {
	    cols++;
	}
    }
    while (ch != '\n' && ch != '\r');
    rewind(fp);
    return (cols);
}


/**
 *
 * Find the number of rows in an FFP
 *
 * @param fp a file poitner to an FFP
 * @return The number of rows
 *
 */


unsigned numRows(FILE * fp)
{
    unsigned int rows = 0;
    char ch;

    while ((ch = getc(fp)) != EOF) {
	if ((ch) == '\n')
	    rows++;
    }
    rewind(fp);
    return rows;
}


/**
 *
 * Test whether this is a (key,value) FFP
 *
 * Perform file test to determine what kind of ffp this
 * is.  Either a columnar or key based FFP.
 *
 * @param fp a file pointer to an FFP file
 * @return true or false
 * @retval 1 Is a (key,value) based FFP
 * @retval 0 Is not (key,value) FFP
 * @todo this could be defined as a macro.
 *
 */


int isKeyBased(FILE * fp)
{
    int ch;
    if (isalpha(ch = fgetc(fp))) {
	ungetc(ch, fp);
	return 1;
    }
    ungetc(ch, fp);

    return 0;

}

/**
 *
 * Determines the key length (feature length) stored
 * in a (key,value) FFP
 *
 * @param fp A file pointer to an FFP file
 * @return Length of the key in nonwhitespace characters
 *
 */


int getKeyLength(FILE * fp)
{
    char s[40];  //@TODO change to a define.
    fscanf(fp,"%s",s);
    rewind(fp);
    return(strlen(s));
}



/**
 *
 * Test whether string contains valid text characters 
 *
 * @param s string to test
 * @param len Length of string
 * @retval 0 if invalid text string
 * @retval 1 if a valid text string
 */



char isValidTxt(register const char *s, int len)
{
    int i = 0;
    while (i < len) {
	if (!isalpha((int) s[i++]))
	    return 0;
    }
    return 1;
}


/**
*
* Creates a random amino acid word of length n
*
* Creates a random amino acid word in s 
* with length n.  Amino acid frequencies
* are the observed frequencies from
* vertebrate proteins (Dyer KF, 1971 J. Biol. Edu. 5, 15-24.)
*
*
* The new word is
* stored in s, space for which must already
* be allocated with enough space to store an
* array of length n.
*
*
*
* @param s a string representing an amino acid sequence
* @param n the length of the string
* @return none
* @todo Edit this function to use the new static hash arays and move into utils.h
*
***************************************/


void randaaword(char *s, int n)
{
    int i, j;
    float chance;
    float freq[] = { .074, .116, .160, .219,
	.252, .310, .347, .421,
	.450, .488, .564, .636,
	.654, .694, .744, .825,
	.887, .900, .933, 1.0001
    };


    s[n] = '\0';
    for (i = 0; i < n; i++) {
	chance = (float) rand() / INT_MAX;
	j = 0;
	while (chance > freq[j])
	    j++;
	switch (j) {
	case 0:
	    s[i] = 'A';
	    break;
	case 1:
	    s[i] = 'R';
	    break;
	case 2:
	    s[i] = 'N';
	    break;
	case 3:
	    s[i] = 'D';
	    break;
	case 4:
	    s[i] = 'C';
	    break;
	case 5:
	    s[i] = 'E';
	    break;
	case 6:
	    s[i] = 'Q';
	    break;
	case 7:
	    s[i] = 'G';
	    break;
	case 8:
	    s[i] = 'H';
	    break;
	case 9:
	    s[i] = 'I';
	    break;
	case 10:
	    s[i] = 'L';
	    break;
	case 11:
	    s[i] = 'K';
	    break;
	case 12:
	    s[i] = 'M';
	    break;
	case 13:
	    s[i] = 'F';
	    break;
	case 14:
	    s[i] = 'P';
	    break;
	case 15:
	    s[i] = 'S';
	    break;
	case 16:
	    s[i] = 'T';
	    break;
	case 17:
	    s[i] = 'W';
	    break;
	case 18:
	    s[i] = 'Y';
	    break;
	case 19:
	    s[i] = 'V';
	    break;
	}
    }
}


//allocates memory and clears to 0
void * chkcalloc(size_t size, size_t n) {
	void * new_block;
	if ((new_block = (void *) calloc(size,n)) == NULL )
		fatal_msg(" %d bytes: Error allocating memory.",size*n);
	return new_block;
}

//allocates memory w/o clearing.
void * chkmalloc(size_t size, size_t n) {
	void * new_block;
	if ((new_block = (void *)malloc(size*n)) == NULL )
		fatal_msg(" %d bytes: Error allocating memory.",size*n);
	return new_block;
}



/* UTILS.C */
