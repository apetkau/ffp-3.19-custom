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
#ifndef _CDFMACROS_H_
#define _CDFMACROS_H_

#define EULER_MASCHERONI 0.5772156649015328606
#define PI_SQRT_6 1.282555555
#define SQRT_2 1.41421
#define evd_cdf(a,b,c) exp(-exp(-((c)-(a))/(b)))
#define norm_cdf(a)    0.5*(1+erf(a/sqrt(SQRT_2)))

/**< @todo Create a separate constants .h file */

#endif				/* _CDFMACROS_H_ */
