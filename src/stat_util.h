/*
 * stat_util.h
 *
 *  Created on: May 19, 2010
 *      Author: dnlennon
 */

#ifndef STAT_UTIL_H_
#define STAT_UTIL_H_

#include "gslvec.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


int multinomial_1sample(
		gsl_rng* rng,
		const_gslvec dist
		);


#endif /* STAT_UTIL_H_ */
