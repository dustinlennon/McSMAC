/*
 * fastmap.h
 *
 *  Created on: Feb 9, 2010
 *      Author: dnlennon
 */

#ifndef FASTMAP_H_
#define FASTMAP_H_

#include "rootedbtree.h"
#include "gslmat.h"
#include "gslvec.h"

void phylo_multimap(
		rootedbtree& t,
		gslmat& rate,
		gslmat& label_mask,
		gslvec& dwell_mask,
		gslmat& site_data,  /* n.sites x n.tips */
		bool lhood_only,
		gslvec* piv,  /* if this is NULL, use stat_dist.  o/w use this as the root dist */
		gslmat& result  /* output */
		);

void phylo_onemap(
		rootedbtree& t,
		gslmat& rate,
		gslmat& label_mask,
		gslvec& dwell_mask,
		gslmat& state_data, /* n.states x n.tips */
		gslmat& result  /* output */
		);

void phylo_edgenames(
		rootedbtree& rbt,
		vec<string>& vs /* output */
		);

void sampletips(
		rootedbtree& rbt,
		gslmat& rate,
		const_gslvec stat_dist,
		int nsim,
		int seed,
		gslmat& result /* output */
		);


void sampletips_sterling_ascertainment(
		rootedbtree& rbt,
		gslmat& rate,
		int nsamp,
		int seed,
		gslmat& output
		);


#endif /* FASTMAP_H_ */
