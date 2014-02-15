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
#ifndef _VSTRING_H_
#define _VSTRING_H_
#include "../config.h"
#define VSTRING PACKAGE_VERSION
#define AUTHORS PACKAGE_AUTHORS
#define YEAR    COPY_YEAR
#define URL	PACKAGE_URL
#define EMAIL PACKAGE_BUGREPORT
#define printVersion() printf("%s %s\nCopyright (C) %s.\n%s\nWritten by %s.\n%s\n",PROG_NAME,VSTRING,YEAR,URL,AUTHORS,EMAIL);

#endif /*_VSTRING_H_*/
