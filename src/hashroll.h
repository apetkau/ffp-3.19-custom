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
/* _ HASHROLL_H_ */
#ifndef _HASHROLL_H_
#define _HASHROLL_H_

#define BUCKETS 65536 /**< The Total number of buckets in feature hash table */
#define hashInc(X) hashAdd((X),(1)) /**< Macro for incrementing a key-value stored in the hash */
#define numKeys(void)  keyN /**< Macro for number of keys in hash */
#define MAX_WORD_SIZE 40
#include <stdbool.h>

/** Linked list for storing values in the hash */

typedef struct node {
    char *key;	    /**< Key value for storing features */
    unsigned value;	    /**< Positive integer value for hash */
    struct node *next;	    /**< Pointer to next node in linked list */
} NODE;


typedef struct hash {
    NODE *table[BUCKETS];     /**< The hash table consisting of an Array of NODES */
    unsigned int st_hash;     /**< A persistent hash value */
    unsigned int rt_hash;     /**< A persistent hash value */
    const unsigned char *hashi;
    bool isClass;
    bool reverse;
    int mode;
    char * r;     /**< last key stored in hash reverse complement*/
    char * s;     /**< last key stored in hash forward direction*/
    int k;		      /**< Length of the previous k-mer */
    int numChar;		/**< state: Was the last hash key valid?*/
			      /**<@todo we can achieve some object orientedness by performing the inithash on this variable and
			       * leaving the function pointers to various kinds of hashes here */
    int keyN;	  /**< The number of elements in the hash table */
    int (*strcmpf) (register const char *, register const char *);
} HASH;





/* prototypes */


void chkpushtxt(HASH * h, char *c, int n);
void pushtxt(HASH * h, char *c, int n);
int hashAddtxt(HASH * h, char *s, unsigned val);

unsigned int sumValues(HASH * h);
void pushaaw(HASH * h, char *c, int n,bool);
void chkpushaaw(HASH * h, char *c, int n,bool);
void chkpushaa(HASH * h, char *c, int n,bool);
void pushaa(HASH * h, char *c, int n,bool);
int hashAddwaa(HASH * h, char *s, unsigned val);
int hashAddaa(HASH * h, char *s, unsigned val);

void printFeatures(HASH * h);
int hashAdd(HASH * h, char *s, unsigned val);
int hashAddw(HASH * h, char *s, unsigned val);
void resetHash(HASH * h);
void init(HASH * h, int isMasked, int mode, bool isClass, bool reverse, int k);
void pushatgc(HASH *, char *, int,bool);
void pushatgcw(HASH *, char *, int,bool);
void chkpushatgc(HASH * h, char *c, int,bool);
void chkpushatgcw(HASH * h, char *c, int n,bool);
int new_strcmpw(register const char *s, register const char *t);
int new_strcmp(const char *s, const char *t);
unsigned int hashValNuc(HASH * h, register char *s);
void hashKeys(HASH *, char ***s);
void hashKeysAndValues(HASH * h, char ***s, unsigned **d);
int numKeys(void);
int freeHash(HASH *);
int hashAssign(HASH * h, register char *s, unsigned int val);
void hashValues(unsigned **values);
int hashDel(char *s);
int hashMax(char *s, unsigned val);
void hashValuesAndSet(HASH * h, unsigned **d);

enum hash_modes { nucleotide, amino, text }; /**< Type of hash to initialize */



// Base characters are mapped so that
// A is the complement of T
// G = 10111100
// C = 01000011

// consider adding RY and ry to here.


static const unsigned char base_hash_values[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,	// A C G  T
    0, 0, 0, 0, 0, 1, 0, 2, 0, 0,
    0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 4, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 2,
    0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 4, 0, 0, 0,
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
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
};




static const unsigned char base_ry_hash_values[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
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
}; /**< Character values for ry hash function */



static const unsigned char aac_values[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 2, 3, 3,
    4, 5, 6, 7, 0, 8, 7, 7, 9, 0,
    10, 8, 8, 11, 11, 0, 7, 4, 0, 4,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 2,
    3, 3, 4, 5, 6, 7, 0, 8, 7, 7,
    9, 0, 10, 8, 8, 11, 11, 0, 7, 4,
    0, 4, 0, 0, 0, 0, 0, 0, 0, 0,
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
}; /**< Character values for classed amino acids */



static const unsigned char aa_hash_values[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 2, 3, 4,
    5, 6, 7, 8, 0, 9, 10, 11, 12, 0,
    13, 14, 15, 16, 17, 0, 18, 19, 0, 20,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 2,
    3, 4, 5, 6, 7, 8, 0, 9, 10, 11,
    12, 0, 13, 14, 15, 16, 17, 0, 18, 19,
    0, 20, 0, 0, 0, 0, 0, 0, 0, 0,
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
}; /**< Character values for amino acids */


static const unsigned char txt_hash_values[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 2, 3, 4, 5,
    6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 0, 0, 0, 0, 0, 0, 1, 2, 3,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
    14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
    24, 25, 26, 0, 0, 0, 0, 0, 0, 0,
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
};

#endif				/* _HASHROLL_H_ */
