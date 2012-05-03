/*
 * stat_util.cc
 *
 *  Created on: May 19, 2010
 *      Author: dnlennon
 */

#include "stat_util.h"

#include "gslvec.h"

#include "vec.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


int multinomial_1sample(
		gsl_rng* rng,
		const_gslvec dist
		)
{
	int K = dist.length();

	const double* d = dist.getConstDoublePtr();
	vec<unsigned int> output(K);
	gsl_ran_multinomial (rng, K, 1, d, output);

	for(int i=0; i<K; i++)
		if(output[i] == 1)
			return i;

	return -1;
}
