/*
 * phylo_util.h
 *
 *  Created on: May 14, 2010
 *      Author: dnlennon
 */

#ifndef PHYLO_UTIL_H_
#define PHYLO_UTIL_H_

#include "rootedbtree.h"
#include "unrootedbtree.h"
#include "gslmat.h"
#include "gslvec.h"

#include <string>
using namespace std;

void treeFromNewickString(const string& s0, rootedbtree& rbt  /* output */);

void precompute(
		graph* t,
		gslmat& rate,
		gslmat& label_mask,
		const_gslvec dwell_mask,
		gslvec& stat_dist /* output */
		);

void precompute(
		graph* t,
		gslmat& rate,
		gslvec& stat_dist /* output */
		);

double fgval_vec_compute(
		graph* g,
		string root_id,
		const_gslvec root_dist,
		gslvec tip_vec /* a vector in {1..n_states}^N_tips */
		);

void fgval_vec_compute(
		graph* g,
		string root_id,
		const_gslvec root_dist,
		gslmat& tip_vec /* a vector in {1..n_states}^N_tips */
		);

void fgval_01mat_compute(
		graph* g,
		string root_id,
		const_gslvec root_dist,
		gslmat& tip_vec /* each column is a vector in {0,1}^n_states */
		);

void stochmapcompute(rootedbtree& rbt);

void smapextract(
		rootedbtree& rbt,
		gslvec result /* output */
		);

void print_phylotree(rootedbtree& t);

void sampletips(
		rootedbtree& rbt,
		gslmat& rate,
		const_gslvec stat_dist,
		int nsim,
		int seed,
		gslmat& result /* output */
		);

void unroot_phylorbt(rootedbtree& rbt, unrootedbtree& ut /* output */);

#endif /* PHYLO_UTIL_H_ */
