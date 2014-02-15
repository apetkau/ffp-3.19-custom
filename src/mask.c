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
#include <stdlib.h>
#include "mask.h"
#include "../config.h"

extern char *weightVector;


/**
 * Create a random mismatch mask
 *
 * Creates a mismatch mask for allowing
 * mismatches in feature comparison.  A 1
 * at position i means that that feature must
 * match; a 0 at position i means that mismatch
 * is allowd.
 *
 * @param s a string representing the mismatch mask
 * @param n the length of the mask
 * @param mismatch the number of mismatches
 * @return none
 */


void randweight(char *s, int n, int mismatch)
{
    int i;

    for (i = 0; i < n; i++)
	weightVector[i] = '1';

    while (mismatch > 0) {
	weightVector[rand() % n] = '0';
	mismatch--;
    }

    weightVector[n] = '\0';
}
