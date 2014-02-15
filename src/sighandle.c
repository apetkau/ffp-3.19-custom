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
#include <signal.h>
#include <unistd.h>
#include "utils.h"
#include "sighandle.h"
#include "../config.h"

static void crash_handler(int sig_num);
static void cleanup();

/* Handle interrupt signals if they are defined
   * segfault, floating point exception, illegal instruction, bad pipe, bus error*/

void initSignalHandlers() 
{ 
//tempfile cleanup	
atexit(cleanup);

#ifdef SIGSEGV
  signal(SIGSEGV, crash_handler);
#endif /* SIGSEGV */
#ifdef SIGFPE
  signal(SIGFPE, crash_handler);
#endif /* SIGFPE */
#ifdef SIGILL
  signal(SIGILL, crash_handler);
#endif /* SIGILL */
#ifdef SIGPIPE
  signal(SIGPIPE, crash_handler);
#endif /* SIGPIPE */
#ifdef SIGBUS
  signal(SIGBUS, crash_handler);
#endif /* SIGBUS */
}

/* Error messages on receiving signal */
extern char * tempfile;

static void cleanup() {
  if (tempfile != NULL)  
  	unlink(tempfile);
}



static void crash_handler(int sig_num)
{ 
  cleanup();

  switch(sig_num) {
#ifdef SIGSEGV
    case SIGSEGV:
      fatal_msg("Segmentation fault.\n");
      break;
#endif /* SIGSEGV */
#ifdef SIGFPE
    case SIGFPE:
      fatal_msg("Floating Point Exception.\n");
      break;
#endif  /* SIGFPE */
#ifdef SIGILL
    case SIGILL:
      fatal_msg("Attempted an illegal instruction.\n");
      break;
#endif  /* SIGILL */
#ifdef SIGPIPE 
    case SIGPIPE:
      fatal_msg("Broken pipe.\n");
      break;
#endif  /* SIGPIPE */
#ifdef SIGBUS
    case SIGBUS:
      fatal_msg("Bus error.\n");
      break;
#endif /* SIGBUS */
  }   
  abort();
}



