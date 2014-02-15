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
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "hashroll.h"
#include "utils.h"
#include "../config.h"


extern bool zflag;
extern bool wflag;

void parseFeatureList(HASH *h, char * filename,int len,int mode) {
    char * s;
    FILE * fp;
    int items;
    unsigned lineno=1;
    long line_offset=0;
    int (*hashAddPtr) (HASH *, char *, unsigned) =NULL;

    if ( (s = (char *) malloc(sizeof(char) * (len + 1))) == NULL)
	fatal_at_line("%s\n",strerror(errno));
    if ((fp = fopen(filename, "r")) == NULL)
	fatal_msg("%s: %s\n", filename, strerror(errno));


    switch (mode) {
    case nucleotide:	    
    	if (zflag || wflag)
		hashAddPtr=hashAddw;
 	else
		hashAddPtr=hashAdd;
    break;
    case amino:
	if (zflag || wflag)
            hashAddPtr=hashAddwaa;
	else
	    hashAddPtr=hashAddaa;
    break;
    case text:
            hashAddPtr=hashAddtxt;
    break;
    }

    errno=0;
    while ( (items=fscanf(fp, "%s", s)) != EOF ) { 
	//fprintf(stderr,"%s %d\n",strerror(errno),items);
	if (errno) 
		fatal_msg("Parse error at line %u char %ld: %s.\n",
			lineno,ftell(fp)-line_offset, strerror(errno) );
	
	switch (items) {
		case 1:
			if (strlen(s) != len)
		    	fatal_msg("Feature lengths of %s differs from -l\n",
			      filename);
			hashAddPtr(h,s,0);
		break;	
		default:
			fatal_msg("Parse error at line %u char %ld.\n",
			lineno,ftell(fp)-line_offset);
		break;
	}
    }
    resetHash(h);
}
