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
/* _UTILS_H_ */
#ifndef _UTILS_H_
#define _UTILS_H_
#include <stdbool.h>

/* prototypes */
unsigned int numCols(FILE * fp);
unsigned int numRows(FILE * fp);
int getKeyLength(FILE * fp);
int isKeyBased(FILE * fp);
char isValid(register const char *s, int len);
char isValidAA(register const char *s, int len);
char isValidTxt(register const char *s, int len);
unsigned long numFeatures(int l);
void fatal_msg(char *fmt, ...);
void fatal_at_line(char *fmt, ...);
void warn_msg(char *fmt, ...);
void printErrorUsageStr();
int isRegularFile(FILE * fp);
FILE *convertPipeToFile(FILE * fp);
int fileno(FILE * fp);
void randaaword(char *s, int n);
void * chkcalloc(size_t size, size_t n);
void * chkmalloc(size_t size, size_t n);
char * basename (const char *name);
int dirExists(char * name);
bool isDirectory(char * fname);

// For some reason this isn't recognized in the header file
// even though I can clearly see it in stdio.h
//int fileno(const FILE *stream);

/* Perform a check of the return value of malloc */

/*
#define smalloc(A) ({\
  void * __SAFEMALLOCPTR;  \
  if ( (__SAFEMALLOCPTR=malloc((A)))==NULL) \
  	fatal_msg("%s:%s: In function \'%s\':%s\n",__FILE__,__LINE__,__func__,strerror(errno)); \
  __SAFEMALLOCPTR ; \
  });
*/

#define printUsageStr() printf(usage_str,PROG_NAME, YEAR, AUTHORS, EMAIL);

#define isBase(c) base_values[(unsigned char)(c)]  /**<Macro: Tests if valid nucleotide character */
#define atgc_to_ry(c) base_rycoded_values[(unsigned char)(c)] /**<Macro: Returns RY class of Nucleotide */
#define ry(s,n) for(int i=0; i < (n); i++) s[i]=atgc_to_ry(s[i]) /**<Macro: Converts a str to ry in place */
#define randword(s,n) for(int i=0; i < n; i++){s[i]=bases[rand()%4]} /**<Macro: Creates random base sequence */
#define complement(s,n)  for (int i=0; i < (n); i++ ) flip((s)[i]) /**<Macro: Converts to complement of Nuct. seq */
#define flip(c) (c)=base_flip_values[(unsigned char)(c)] /**<Performs base flip for complement. */
#define rev(s, n) \
		{ \
		   for (int i=0; i < (n)/2; i++)   \
		   	{                   \
			s[(n)]=s[i];        \
			s[i]=s[(n)-i-1];    \
			s[(n)-i-1]=s[(n)];  \
			}                   \
		   s[(n)]='\0';             \
	  	   } /**<Reverses a sequence in place */


#define aa_to_class(c) aa_classcoded_values[(unsigned char)(c)]	/**<Macro: Converts AA to proper class character */
#define _class(s,n) for(int i=0; i < n; i++) s[i]=aa_to_class(s[i]) /**<Converts string to AA classes */
#define isAA(c) aa_values[(unsigned char)(c)] /**<Test whether character is an AA */

static const unsigned char bases[] = { 84, 67, 71, 84 };

static const unsigned char base_flip_values[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 84, 0, 71, 0, 0,	//A=65 C=67
    0, 67, 0, 0, 0, 0, 0, 0, 0, 0,	//G=71 
    0, 0, 89, 0, 65, 65, 0, 0, 0, 82,	//R=82 T=84 U=85 Y=89
    0, 0, 0, 0, 0, 0, 0, 84, 0, 71,
    0, 0, 0, 67, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 89, 0, 65, 65, 0, 0,
    0, 82, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
}; /**< ASCII values for base flips, equiv. to: tr [AaTtGgCcRrYy] [TTAACCGGYYRR] */


static const bool base_values[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 1, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 1, 1, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
}; /**< True if ASCII numeric subscript value is a base */


static const unsigned char base_rycoded_values[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 82, 0, 89, 0, 0,	//A=65 C=67
    0, 82, 0, 0, 0, 0, 0, 0, 0, 0,	//G=71 
    0, 0, 82, 0, 89, 89, 0, 0, 0, 89,	//R=82 T=84 U=85 Y=89
    0, 0, 0, 0, 0, 0, 0, 82, 0, 89,
    0, 0, 0, 82, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 82, 0, 89, 89, 0, 0,
    0, 89, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
}; /**< Supplies proper ascii values for RY conversions */

//(TS),(DE),(QKR),(VILM),(WFY),C,G,A,N,H,P

static const unsigned char aa_classcoded_values[] = {

    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 65, 0, 67, 68, 68,
    70, 71, 72, 73, 0, 75, 73, 73, 78, 0,
    80, 75, 75, 83, 83, 0, 73, 70, 0, 70,
    0, 0, 0, 0, 0, 0, 0, 65, 0, 67,
    68, 68, 70, 71, 72, 73, 0, 75, 73, 73,
    78, 0, 80, 75, 75, 83, 83, 0, 73, 70,
    0, 70, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
};  /**< Supplie proper ascii values for AA class conversions */


static const bool aa_values[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 1, 1, 1,	//A=65 C=67
    1, 1, 1, 1, 0, 1, 1, 1, 1, 0,	//G=71 
    1, 1, 1, 1, 1, 0, 1, 1, 0, 1,	//R=82 T=84 U=85 Y=89
    0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
    1, 1, 1, 1, 1, 1, 0, 1, 1, 1,
    1, 0, 1, 1, 1, 1, 1, 0, 1, 1,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
}; /**< True if array subscript is a valid AA character */




#endif				/* _UTILS_H_ */
