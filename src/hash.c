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
#include <ctype.h>
#include <math.h>
#include "hash.h"
#include "../config.h"



extern char *weightVector;
extern int Length;


/** @todo Mysteries to solve, why is abs() necessary in the hash functions? 
 *  Using unsigned int types should remove the need for this.
*/

/* Functions private to this file and therefore not part of hash.h */
/* This prevents anyone from accessing the hash improperly */


static int hashrywi(register const char *s);
static int hashryi(register const char *s);
static int hashatgcwi(register const char *s);
static int hashatgci(register const char *s);
static int new_strcmp(register const char *s, register const char *t);
static int new_strcmpw(register const char *s, register const char *t);
static int hashaawi(register const char *s);
static int hashaai(register const char *s);
static int hashaacwi(register const char *s);
static int hashaaci(register const char *s);
static int hashtxti(register const char *s);


/**
 *
 * Initializing feature that creates an empty hash table.
 *
 * This function initializes the hash table array by setting all
 * pointers in the the table array to NULL, setting the number of keys
 * keyN to zero and setting up the function pointer to the 
 * appropriate hash and string comparison functions.  The
 * function pointer is necessary to provide polymorphic
 * behavior depending on whether we are using a feature mask, RY
 * coded, ATGC coded or amino acid sequence
 *
 * Note however it does not free any memory, so freeHash should
 * be called before reinitializing the hash.
 *
 *
 * @param isMasked true if a mask is used
 * @param isNt true if using nucleic acid features
 * @param isClass true if using amino acid classes
 * @return None
 * @todo Since, everytime as hash is used it is first intialized, there may not be a need for a BUCKET variable instead the size of the hash can be dynamically allocated.
 */




void initHash(int isMasked, int mode, int isClass)
{

    // assign function pointers depending
    // on whether we are using feature
    // character masks or not
    if (mode == nucleotide) {
	if (isMasked)		// is a NA hash
	{
	    if (isClass)
		hashf = &hashrywi;
	    else
		hashf = &hashatgcwi;
	    strcmpf = &new_strcmpw;

	} else {
	    if (isClass)
		hashf = &hashryi;
	    else
		hashf = &hashatgci;
	    strcmpf = &new_strcmp;

	}
    }

    else if (mode == amino) {
	if (isMasked) {
	    if (isClass)
		hashf = &hashaacwi;
	    else
		hashf = &hashaawi;
	    strcmpf = &new_strcmpw;
	} else {
	    if (isClass)
		hashf = &hashaaci;
	    else
		hashf = &hashaai;

	    strcmpf = &new_strcmp;
	}
    } else			// is for text
    {
	hashf = &hashtxti;
	strcmpf = &new_strcmp;
    }

    int i;
    for (i = 0; i < BUCKETS; i++)
	table[i] = NULL;
    keyN = 0;
}







/**
 *
 * Returns an unmasked hash value for an RY coded feature
 * 
 * This is the hashing function used to store
 * features in the hash table.  The hashing function
 * has a direct affect on the collision rate of
 * hash stores.  The more hash collisions the less
 * effective the hash lookups are. 
 *
 *
 * @param s An RY-coded character string.
 * @retval >= 0 A valid hash index
 * @retval -1 s contains an invalid character
 */



inline static int hashryi(register const char *s)
{
    long unsigned hash = 0;
    register int index;

    while ( *s != '\0') {
	hash *= 3;
	index = base_ry_hash_values[(unsigned char) (*s)];
	s++;
	if (!index)
	    return -1;
	hash += index;
	hash %= BUCKETS;
    }
    return (int) hash;
}



/**
 *
 * Returns a masked hash value for an RY coded feature
 * 
 * This is the hashing function used to store
 * features in the hash table.  The hashing function
 * has a direct affect on the collision rate of
 * hash stores.  The more hash collisions the less
 * effective the hash lookups are. Uses a feature
 * mask, ignoring masked out features
 *
 *
 * @param s An RY-coded character string
 * @retval >= 0 A valid hash index
 * @retval -1 s contains an invalid character
 */



static int hashrywi(register const char *s)
{
    long unsigned hash = 0;
    int i = 0;
    register int index;

    while (s[i] != '\0') {
	if (weightVector[i] == '1') {
	    hash *= 3;
	    index = base_ry_hash_values[(unsigned char) s[i]];
	    if (!index)
		return -1;
	    hash += index;
	    hash %= BUCKETS;
	}
	i++;
    }
    return (int) hash;
}


static const int base_hash_values[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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

static int hashatgci(register const char *s)
{
    long unsigned hash = 0;
    int i = 0;
    register int index;

    while (s[i] != '\0') {
	hash *= 5;
	index = base_hash_values[(unsigned char) s[i]];
	if (!index)
	    return -1;
	hash += index;
	hash %= BUCKETS;
	i++;
    }
    return (int) hash;
}


static int hashatgcwi(register const char *s)
{
    long unsigned hash = 0;
    int i = 0;
    register int index;

    while (s[i] != '\0') {
	if (weightVector[i] == '1') {
	    hash *= 5;
	    index = base_hash_values[(unsigned char) s[i]];
	    if (!index)
		return -1;
	    hash += index;
	    hash %= BUCKETS;
	}
	i++;
    }
    return (int) hash;
}


/**
 *
 * Masked String comparison method used internally by the hash table functions.
 *
 * This function differs slightly from the library version of strcmp
 * It short circuits when it finds the first difference in strings
 * It is also declared inline so that the compiler has the option of
 * substituting this code in directly.  It is declared static so that
 * it can only be used by the hash functions.
 * Uses a feature mask, ignoring masked out features
 *
 * @param s a null terminated string
 * @param t a null terminated string
 * @retval 1 if s and t are equal
 * @retval 0 if s and t are different
 *
 */


static int new_strcmpw(register const char *s, register const char *t)
{
    int i = 0;
    while (s[i] != '\0') {
	if (weightVector[i] == '1')
	    if (s[i] != t[i])
		return 1;
	i++;
    }
    return 0;
}


/**
 *
 * String comparison method used internally by the hash table functions.
 *
 * This function differs slightly from the library version of strcmp
 * It short circuits when it finds the first difference in strings
 * It is also declared inline so that the compiler has the option of
 * substituting this code in directly.  It is declared static so that
 * it can only be used by the hash functions.
 *
 * @param s a null terminated string
 * @param t a null terminated string
 * @retval 1 if s and t are equal 
 * @ret val 0 if s and t are different
 *
 */



inline static int new_strcmp(const char *s, const char *t)
{
    while (*s != '\0') {
	if (*(s++) != *(t++))
	    return 1;
    }
    return 0;
}

/**
 *
 * Adds a key-value pairs to the hash table
 *
 *
 * Checks to see if s exists in the hash table, if it
 * doesn't add s to the hash table and store a value
 * of val.  If feature s exists in the hash table then
 * increment the frequency. The function has the 
 * following return values:
 * 0 or -1 if the feature isn't in the table
 * 1 if the feature is in the table
 * -1 more precisely indicates a hash collision has
 *  occurred.  i.e. two features hashed to the same
 *  key value, but are not the same.
 *
 *  Note that this function uses function pointers
 *  to a hashing function and a string comparison
 *  function.  This is to create polymorphic behavior
 *  depending upon whether we are hashing nucleotides
 *  or hashing proteins or using feature masking. 
 *
 * @param s a pointer to a string
 * @param val the integer value to add to the current value 
 * @retval 0 For feature not in hash 
 * @retval -1 For feature not in hash but hash collision,
 * @retval 1 For feature in hash
 * @see hashInc
 * @see strcmpf
 * @see hashf
 *
 *
 */


/*
int hashAdd(char *s, unsigned val)
{
    NODE *ptr;
    NODE *r;
    int index;

    if ((index = (*hashf) (s)) < 0)	
	return 0;


    ptr = table[index];

    while (ptr != NULL) {
	if ((*strcmpf) (ptr->key, s) == 0) {
	    ptr->value += val;
	    return 1;
	}

	ptr = ptr->next;
    }

    r = (NODE *) malloc(sizeof(NODE));
    r->key = (char *) malloc(sizeof(char) * (Length + 1));
    strcpy(r->key, s);
    r->value = val;
    r->next = table[index];
    table[index] = r;
    keyN++;
    return -1;

}
*/





/**
 *
 * The hash value associated with key s
 *
 * Returns the frequency stored in the hash which is associated with
 * feature s.
 *
 * @param s A pointer to a Key value s which is a string
 * @return an integer value associated the frequency of s
 * @retval >=0 if s is in hash
 * @retval 0 if s is not in hash
 */

int hashval(register char *s)
{
    NODE *ptr;

    ptr = table[(*hashf) (s)];


    while (ptr != NULL) {
	if ((*strcmpf) (ptr->key, s) == 0)
	    return ptr->value;
	ptr = ptr->next;
    }

    return 0;
}




/**
 *
 * Assign a new value to the key value s
 *
 * Returns the frequency stored in the hash which is associated with
 * feature s.
 *
 * Note: need to add checking for whether has key is stored at all!
 *
 * @param s A pointer to a Key value s which is a string
 * @param n An integer which is the new value
 * @return true or false
 * @retval 1 if found
 * @retval 0 if not found
 *
 */

/* Moved to Macro definition */

/*
void hashAssign(char *s, int n)
{
    NODE *ptr;

    ptr = table[(*hashf) (s)];

    while (ptr != NULL) {
	if ((*strcmpf) (ptr->key, s) == 0) {
		ptr->value = n;
	    return;
	}
	ptr = ptr->next;
    }

    return;
}
*/


/**
 *
 * returns a list of keys stored in the hash table
 *
 * This is function is used to return the values of the
 * keys hashed into the table.  No error checking is
 * provided to make sure that enough memory has been
 * allocated to store all the keys into the array of
 * strings.
 *
 * @param s a pointer to an array of character arrays.
 * @return None
 *
 */


void hashKeys(char ***s)
{
    int i, j;
    NODE *ptr;
    char ** sp;
    sp = (char **) malloc(sizeof(char *) * keyN);
    for (i = 0; i < keyN; i++)
	sp[i] = (char *) malloc(sizeof(char) * (Length + 1));
    i = 0;
    for (j = 0; j < BUCKETS; j++)
	if (table[j] != NULL) {
	    ptr = table[j];
	    while (ptr != NULL) {
		strcpy(sp[i++], ptr->key);
		ptr = ptr->next;
	    }
	}
   *s=sp;
}


/**
 *
 * Returns the stored hash values
 *
 * This function is used to return the values
 * stored in the hash table. The ordering is
 * not orderd by key in any way, but is retrieved
 * by increasing hash index.  This function can
 * be used to save the time involved in retrieving
 * all the keys and then looking up the value associated
 * with a key - one at a time.  This is especially useful
 * if the keys remain in the hash throughout and
 * ordering isn't essential.
 *
 * @param values
 * @return none
 *
 */




void hashValues(unsigned **values)
{
    int i, j;
    NODE *ptr;
    *values = (unsigned *) malloc(sizeof(unsigned *) * keyN);
    i = 0;
    for (j = 0; j < BUCKETS; j++)
	if (table[j] != NULL) {
	    ptr = table[j];
	    while (ptr != NULL) {
		*values[i++] = ptr->value;
		ptr = ptr->next;
	    }
	}

}




/**
 *
 * Frees the memory allocated to store the hash
 *
 * All memory allocated with malloc is freed from the
 * hash.  All values are deleted. returns 0 on success.
 * The number of hash keys is also returned to zero.
 *
 * @param None
 * @retval 1 on success
 *
 *
 */


int freeHash(void)
{
    int i;
    NODE *ptr;
    NODE *last;
    for (i = 0; i < BUCKETS; i++) {
	ptr = table[i];
	while (ptr != NULL) {
	    last = ptr;
	    ptr = ptr->next;
	    free(last->key);
	    free(last);
	}
	table[i] = NULL;
    }
    keyN = 0;
    return 1;
}

/**
 *
 *  Hash function for masked classed amino acid sequences
 *
 *  Returns a hash index, for char string s
 *  Positions which are masked out in w are
 *  not used to calculate index.
 *  String s must contain a valid character,
 *  one of the Amino acid classes specified in
 *  _class(), otherwise a hash of -1 is returned.
 *  
 *  Depends on global variable string weightVector
 *
 *  @param s hash key
 *  @return  returns the hash index
 *  @retval -1 invalid class character
 *  @retval >0 a valid hash index
 *  @see _class()
 *  @todo Consider the possibility that classing be built into the hash functions, so that conversion to classes is unnecessary.  This eliminates a step.
 *
 ************************************************/



static int hashaacwi(register const char *s)
{
    long unsigned hash = 0;
    int i = 0;
    register int index;

    while (s[i] != '\0') {
	if (weightVector[i] == '1') {
	    hash *= 12;
	    index = aac_values[(unsigned char) s[i]];
	    if (!index) {
		return -1;
	    }
	    hash += index;
	}
	i++;
    }
    return (int) abs(hash) % BUCKETS;
}



/**
 *
 *  Hash function for amino acid sequences
 *
 *  Returns a hash index, for char string s
 *  Positions which are masked out in weightVector are
 *  not used to calculate index.
 *  String s must contain a valid character,
 *  one of the Amino acids.
 *  
 *
 *  @param s hash key
 *  @return  returns the hash index
 *  @retval -1 Invalid Amino acid class character
 *  @retval >0 A valid hash index
 *
 ************************************************/




static int hashaaci(register const char *s)
{
    long unsigned hash = 0;
    int i = 0;
    register int index;

    while (s[i] != '\0') {
	hash *= 12;
	index = aac_values[(unsigned char) s[i]];
	if (!index) {
	    return -1;
	}
	hash += index;
	i++;
    }
    return (int) abs(hash) % BUCKETS;
}


/**
 *
 *  Hash function for amino acid sequences
 *
 *  Returns a hash index, for char string s
 *  Positions which are masked out in weightVector are
 *  not used to calculate index.
 *  String s must contain a valid character,
 *  one of the Amino acids.
 *  
 *
 *  @param s hash key
 *  @return  returns the hash index
 *  @retval -1 Invalid Amino acid class character
 *  @retval >0 A valid hash index
 *
 ************************************************/




static int hashaawi(register const char *s)
{
    long unsigned hash = 0;
    int i = 0;
    register int index;

    while (s[i] != '\0') {
	if (weightVector[i] == '1') {
	    hash *= 21;
	    index = aa_hash_values[(unsigned char) s[i]];
	    if (!index) {
		return -1;
	    }
	    hash += index;
	}
	i++;
    }
    return (int) abs(hash) % BUCKETS;
}



/**
 *
 *  Hash function for amino acid sequences
 *
 *  Returns a hash index, for char string s
 *  String s must contain a valid amino acid character,
 *  otherwise a hash of -1 is returned.
 *  
 *  @param s hash key
 *  @retval -1 invalid amino acid character
 *  @retval >0 a valid hash index
 *  @return  returns the hash index
 *
 **************************************************/




static int hashaai(register const char *s)
{
    long unsigned hash = 0;
    int i = 0;
    register int index;

    while (s[i] != '\0') {
	hash *= 21;
	index = aa_hash_values[(unsigned char) s[i]];
	if (!index) {
	    return -1;
	}
	hash += index;
	i++;
    }
    return (int) abs(hash) % BUCKETS;
}











/**
 *
 *  Hash function for text sequences
 *
 *  String s must contain a valid text character,
 *  
 *  @param s hash key
 *  @retval -1 invalid text character
 *  @retval >0 a valid hash index
 *  @return  returns the hash index
 *
 **************************************************/



static int hashtxti(register const char *s)
{
    long unsigned hash = 0;
    int i = 0;
    register int index;

    while (s[i] != '\0') {
	hash *= 27;
	index = txt_hash_values[(unsigned char) s[i]];
	if (!index) {
	    return -1;
	}
	hash += index;
	i++;
    }
    return (int) abs(hash) % BUCKETS;
}


/**
 *
 * Adds a key-value pair to the hash table if greate than
 * current value
 *
 *
 * Checks to see if s exists in the hash table, if it
 * doesn't add s to the hash table and store a value
 * of val.  If feature s exists in the hash table then
 * increment the frequency. The function has the 
 * following return values:
 * 0 or -1 if the feature isn't in the table
 * 1 if the feature is in the table
 * -1 more precisely indicates a hash collision has
 *  occurred.  i.e. two features hashed to the same
 *  key value, but are not the same.
 *
 *  Note that this function uses function pointers
 *  to a hashing function and a string comparison
 *  function.  This is to create polymorphic behavior
 *  depending upon whether we are hashing nucleotides
 *  or hashing proteins or using feature masking. 
 *
 * @param s a pointer to a string
 * @param val the integer value to add to the current value 
 * @retval 0 For feature not in hash 
 * @retval -1 For feature not in hash but hash collision,
 * @retval 1 For feature in hash
 * @see hashInc
 * @see strcmpf
 * @see hashf
 *
 *
 */


int hashMax(char *s, unsigned val)
{
    NODE *ptr;
    NODE *r;
    int index;

    if ((index = (*hashf) (s)) < 0)	/* If an invalid hash */
	return 0;


    ptr = table[index];

    while (ptr != NULL) {
	if ((*strcmpf) (ptr->key, s) == 0) {
	    ptr->value = ((val > ptr->value) ? val : ptr->value);
	    return 1;
	}

	ptr = ptr->next;
    }


    r = (NODE *) malloc(sizeof(NODE));
    r->key = (char *) malloc(sizeof(char) * (Length + 1));
    strcpy(r->key, s);
    r->value = val;
    r->next = table[index];
    table[index] = r;
    keyN++;
    return -1;

}




int hashDel(char *s)
{
    NODE *ptr;
    NODE *r;
    int index;

    if ((index = (*hashf) (s)) < 0)	/* If an invalid hash */
	return 0;


    ptr = table[index];
    r = table[index];

    while (ptr != NULL) {
	if ((*strcmpf) (ptr->key, s) == 0) {
	    if (r == table[index]) {
		table[index] = ptr->next;
	    } else
		r->next = ptr->next;
	    free(ptr);
	    keyN--;

	    break;
	}
	r = ptr;
	ptr = ptr->next;
    }


    return 0;
}
