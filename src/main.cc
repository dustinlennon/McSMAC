/*
 * main.cc
 *
 *  Created on: May 4, 2010
 *      Author: dnlennon
 */

#include "test.h"
#include "graphalg.h"
#include "gslmat.h"
#include "gslvec.h"
#include "phylo_util.h"
#include "phylo_objs.h"

#include "rootedbtree.h"
#include "unrootedbtree.h"

#include "RecoverableAssertException.h"

#include "vec.h"
#include "graph.h"

#include "dataobject.h"

#include <gsl/gsl_blas.h>

#include <math.h>
#include <stdlib.h>

#include <string>
#include <fstream>
using namespace std;

void reopen(ifstream& f, string s);
void recycle(ifstream& f, string s);

int main(int argc, char** argv)
{
	run_container_tests(argc, argv);

	string tree_string;
	ifstream fin ("data/yang.tree");
	fin>>tree_string;

	rootedbtree rbt;
	treeFromNewickString(tree_string, rbt);

	unrootedbtree ut;
	unroot_phylorbt(rbt, ut);

	cout<<endl<<ut<<endl;

	reopen(fin, "data/rate.mat");
	gslmat rate(fin);

    reopen(fin, "data/label.mask.mat");
    gslmat label_mask( fin );

    reopen(fin, "data/dwell.mask.mat");
	gslvec dwell_mask( fin );

	reopen(fin, "data/data.yang");
	gslmat yang_data( fin );

	gslvec stat_dist;
	precompute(&ut, rate, stat_dist);

	fgval_vec_compute(&ut, "td1", stat_dist, yang_data);
	node* np;
	ut.findNode("td1", np);
	UphyloNodeObject* ndop = dynamic_cast<UphyloNodeObject*>(np->getDOP().ptr());
	double lhood;
	gsl_blas_ddot(ndop->fval, stat_dist, &lhood);
	cout<<lhood<<endl;

	fgval_vec_compute(&ut, "i1", stat_dist, yang_data);
	ut.findNode("i1", np);
	ndop = dynamic_cast<UphyloNodeObject*>(np->getDOP().ptr());
	gsl_blas_ddot(ndop->fval, stat_dist, &lhood);
	cout<<lhood<<endl;

	return 0;
}
