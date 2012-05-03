/*
 * fastmap.cc
 *
 *  Created on: Feb 9, 2010
 *      Author: dnlennon
 */

#include "fastmap.h"
#include "vec.h"
#include "graph.h"
#include "stat_util.h"

#include "rootedbtree.h"

#include "phylo_util.h"
#include "phylo_objs.h"
#include "phylo_rbtalg.h"

#include "RecoverableAssertException.h"

#ifdef ENABLE_MCWQ
#include "MultiCoreWorkQueue.h"
#else 
#define MAX_THREADS 1
#include "SingleCoreWorkQueue.h"
#endif

#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

void sampletips_noprecompute(
		rootedbtree& rbt,
		gslmat& rate,
		const_gslvec stat_dist,
		int nsim,
		int seed,
		gslmat& result /* output */
		);

int gvv_compare (const void *avp, const void *bvp)
{
	gsl_vector_view *a = (gsl_vector_view*)avp;
	gsl_vector_view *b = (gsl_vector_view*)bvp;
	int len = a->vector.size;
	for(int i=0; i<len; i++)
	{
		double A = gsl_vector_get(&(a->vector), i);
		double B = gsl_vector_get(&(b->vector), i);
		if(A < B)
			return -1;
		if(A > B)
			return 1;
	}
	return 0;
}

/* Work Queue Function Definitions
 * 		Our data structure for the work queue function evaluator
 */

typedef
struct smap_data
{
	rootedbtree* rbt;
	gsl_vector_view* gvv;
	gslvec* sdist;
	bool lhood_only;
	vec<size_t>* bundle;
} smap_data;

/*
 * 		Our function for the work queue evaluations
 */
void smap_fn(void* ans, void* data)
{
	gslmat& result = *((gslmat*)ans);
	smap_data* my_data = (smap_data*)data;
	vec<size_t>& vuip = *((vec<size_t>*) my_data->bundle);

	gsl_vector_view* gvv = my_data->gvv;

	for(int i=0; i<vuip.length(); i++)
	{
		gsl_vector_view row_view = gsl_matrix_row (result, vuip[i]);

		double lhood = fgval_vec_compute(my_data->rbt, my_data->rbt->getRootId(), *(my_data->sdist), gvv[ vuip[i] ]);
		if(my_data->lhood_only == false)
		{
			stochmapcompute(*(my_data->rbt));
			smapextract(*(my_data->rbt), row_view);
		}
		else
		{
			gsl_vector_set(&row_view.vector, 0, lhood);
		}
	}
}

/*
 * Monolithic function that computes the likelihood, edge labeled jumps,
 * and edge dwelling time site-by-site, conditional on the site data.
 */
void phylo_multimap(
		rootedbtree& rbt,
		gslmat& rate,
		gslmat& label_mask,
		gslvec& dwell_mask,
		gslmat& site_data,
		bool lhood_only,
		gslvec* piv,  /* if this is NULL, use stat_dist.  o/w use this as the root dist */
		gslmat& result  /* output */
		)
{
	vec<edge*> vep;
	rbt.edgeList(vep);
	int nr = site_data.nrow();
	int nc = 1 + 2* vep.length();

	if(lhood_only == false)
		result.setsize(nr , nc);
	else
		result.setsize(nr, 1);

	gslvec stat_dist;
	precompute(&rbt, rate, label_mask, dwell_mask, stat_dist);

	if(piv != NULL)
	{
		gsl_vector_memcpy(stat_dist, *piv);
	}

	/* Sort the site_data by dictionary order:
	 * 		If adjacent sites are identical, then the results can be
	 * 		copied instead of recomputed.
	 *
	gsl_vector_view* gvv = new gsl_vector_view[site_data.nrow()];
	size_t* idx = new size_t[site_data.nrow()];
	for(int r=0; r<nr; r++)
	{
		gvv[r] = gsl_matrix_row(site_data, r);
		idx[r] = r;
	}
	gsl_heapsort_index (idx, gvv, site_data.nrow(), sizeof(gsl_vector_view), gvv_compare);

	 * Pre-process data for use in a WorkQueue:
	 * 		loop 1..site_data.nrow:  get indices of unique tip vectors
	 *
	queue<size_t> idx_q;
	int num_unique_sites = 0;
	for(int r=0; r<site_data.nrow(); r++)
	{
		//gsl_vector_view row_view = gsl_matrix_row (result, idx[r]);
		if( r > 0 && gvv_compare(&gvv[idx[r]], &gvv[idx[r-1]]) == 0 )
		{
			continue;
		}

		idx_q.insert(idx[r]);
		num_unique_sites++;
	}
	vec<size_t> idx_v;
	idx_v.assign(idx_q);
	*

	// Split up workload:
	int ncomp = num_unique_sites;
	*/
	int ncomp = site_data.nrow();
	vec<size_t> idx_v(ncomp);
	gsl_vector_view* gvv = new gsl_vector_view[site_data.nrow()];

	for(int ii=0; ii<idx_v.length(); ii++)
		{
			gvv[ii] = gsl_matrix_row(site_data, ii);
			idx_v[ii] = ii;
		}

	int NumActiveThreads = ncomp < MAX_THREADS ? ncomp : MAX_THREADS;
	int MinPerBin = floor(ncomp / MAX_THREADS);
	int NumLeftOver;

	if(ncomp < MAX_THREADS)
		NumLeftOver = ncomp;
	else
		NumLeftOver = ncomp % (MinPerBin * MAX_THREADS);

	vec<size_t>* bundles = new vec<size_t>[NumActiveThreads];
	rootedbtree rbts[NumActiveThreads];

	int PerBin = MinPerBin;
	int curpos = 0;
	for(int i=0; i<NumActiveThreads; i++)
	{
		PerBin = MinPerBin;
		if(NumLeftOver > 0)
		{
			PerBin = MinPerBin + 1;
			NumLeftOver--;
		}

		rbts[i] = rbt;
		rbts[i].allDopMakeOwner();
		//vec<size_t>& vb = bundles[i];
		bundles[i].resize(PerBin);

		for(int j=0; j<PerBin; j++)
		{
			bundles[i][j] = idx_v[curpos+j];
		}
		curpos = curpos + PerBin;
	}

	// Create the work queue
#ifdef ENABLE_MCWQ
	WorkQueue* myWQ = new MultiCoreWorkQueue(smap_fn);
#else
	WorkQueue* myWQ = new SingleCoreWorkQueue(smap_fn);
#endif

	smap_data* my_data = new smap_data[NumActiveThreads];
	for(int i=0; i<NumActiveThreads; i++)
	{
		my_data[i].rbt = &rbts[i];
		my_data[i].gvv = gvv;
		my_data[i].sdist = &stat_dist;
		my_data[i].lhood_only = lhood_only;
		my_data[i].bundle = &bundles[i];
		myWQ->enqueue(&result, &my_data[i]);
	}
	myWQ->results();

	delete myWQ;
	delete [] my_data;
	delete [] bundles;

	/*
	for(int r=0; r<site_data.nrow(); r++)
	{
		gsl_vector_view row_view = gsl_matrix_row (result, idx[r]);
		if( r > 0 && gvv_compare(&gvv[idx[r]], &gvv[idx[r-1]]) == 0 )
		{
			gsl_vector_view prev_row_view = gsl_matrix_row (result, idx[r-1]);
			gsl_vector_memcpy(&row_view.vector, &prev_row_view.vector);
		}
	}
	delete[] idx;
	*/

	delete[] gvv;

	return;
}


void phylo_onemap(
		rootedbtree& rbt,
		gslmat& rate,
		gslmat& label_mask,
		gslvec& dwell_mask,
		gslmat& state_data, // n.states x n.tips
		gslmat& result  // output
		)
{
	vec<edge*> vep;
	rbt.edgeList(vep);
	int nc =1 + 2 * vep.length();
	result.setsize(1, nc);

	gsl_vector_view gvv = gsl_matrix_row(result, 0);

	// map the tips and edges of the tree into vectors
	gslvec stat_dist;
	precompute(&rbt, rate, label_mask, dwell_mask, stat_dist);

	fgval_01mat_compute(&rbt, rbt.getRootId(), stat_dist, state_data);
	stochmapcompute(rbt);
	smapextract(rbt, gvv);
}

//////////////////////////////////////////////////////////////////////////////////////

class phyloedgenames_Fn : public graphTraversalFunction
{
public:
	phyloedgenames_Fn(vec<string>& vs) : myvs(vs) {}

	virtual void visit(edge* e, node* prev)
	{
		string s = e->getId();
		myqs.insert(s);
	}

	virtual void finalize(void)
	{
		myvs.assign(myqs);
	}

	vec<string>& myvs;
	queue<string> myqs;
};


void phylo_edgenames(
		rootedbtree& rbt,
		vec<string>& vs /* output */
		)
{
	/*
	vec<edge*> vep;
	rbt.edgeList(vep);

	vs.resize(vep.length());

	for(int i=0; i<vep.length(); i++)
	{
		vs[i] = vep[i]->getId();
	}
	*/

	phyloedgenames_Fn phyloedgenames_Fn(vs);
	phylo_rbt_postorder(rbt, &phyloedgenames_Fn);
}


//////////////////////////////////////////////////////////////////////////////////////

class sampletips_Fn : public graphTraversalFunction
{
public:
	sampletips_Fn(gsl_rng * rng) : myrng(rng) {}
	virtual void visit(edge* e, node* prev)
	{
		phyloNodeObject* prev_ndop = dynamic_cast<phyloNodeObject*> ( prev->getDOP().ptr() );
		int state = prev_ndop->sim;

		phyloEdgeObject* edop = dynamic_cast<phyloEdgeObject*> ( e->getDOP().ptr() );
		gsl_vector_view gvv = gsl_matrix_row(edop->CP_mat, state);

		gslvec gvtmp(gvv);

		node* next = e->otherNode(prev);
		phyloNodeObject* next_ndop = dynamic_cast<phyloNodeObject*> ( next->getDOP().ptr() );

		int next_state = multinomial_1sample(myrng, gvv);
		next_ndop->sim = next_state;
	}

	gsl_rng * myrng;
};

void sampletips(
		rootedbtree& rbt,
		gslmat& rate,
		const_gslvec root_dist,
		int nsim,
		int seed,
		gslmat& result /* output */
		)
{
	gslvec tmp;
	precompute(&rbt, rate, tmp);
	sampletips_noprecompute(rbt, rate, root_dist, nsim, seed, result);
}


void sampletips_noprecompute(
		rootedbtree& rbt,
		gslmat& rate,
		const_gslvec stat_dist,
		int nsim,
		int seed,
		gslmat& result /* output */
		)
{
	// initialize random number generator with seed
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_ranlxd2);
	gsl_rng_set(rng, seed);

	// resize the result to match the output
	phyloGraphObject* gdop = dynamic_cast<phyloGraphObject*>(rbt.getDOP().ptr());
	int ntips = gdop->tip_strings.length();
	result.setsize(nsim, ntips);

	node* root_np;
	rbt.findNode( rbt.getRootId(), root_np );
	phyloNodeObject* ndop = dynamic_cast<phyloNodeObject*>( root_np->getDOP().ptr() );

	// generate the required number of tip vectors
	for(int i=0; i<nsim; i++)
	{
		gsl_vector_view gvv = gsl_matrix_row(result, i);

		int root_state = multinomial_1sample(rng, stat_dist);
		ndop->sim = root_state;

		sampletips_Fn sampletips_Fn(rng);
		bfs(rbt, rbt.getRootId(), &sampletips_Fn);

		for(int j=0; j<ntips; j++)
		{
			node* tip_np;
			rbt.findNode(gdop->tip_strings[j], tip_np);
			phyloNodeObject* tip_ndop = dynamic_cast<phyloNodeObject*>(tip_np->getDOP().ptr() );

			gsl_vector_set(&gvv.vector, j, tip_ndop->sim + 1);
		}
	}

	gsl_rng_free(rng);
}


///////////////////////////////////////////////////////////////////////////

void sampletips_sterling_ascertainment(
		rootedbtree& rbt,
		gslmat& rate,
		int nsamp,
		int seed,
		gslmat& output
		)
{
	gslvec stat_dist;
	precompute(&rbt, rate, stat_dist);

	graph g = rbt;
	g.allDopMakeOwner();

	node* np;
	if( g.findNode("hg18", np) == false )
	{
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	}

	vec<string> vs;
	np->getAdjNodeStrings(vs);
	if( vs.length() != 1)
	{
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	}

	edge* ep;
	if( g.findEdge(np->getId(), vs[0], ep ) == false )
	{
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	}

	node* new_root = ep->otherNode(np);

	phyloEdgeObject *edop = dynamic_cast<phyloEdgeObject*>(ep->getDOP().ptr());
	gsl_vector_view gvv = gsl_matrix_row(edop->CP_mat, 1);
	gslvec init_dist(2);
	gsl_vector_memcpy(init_dist, &gvv.vector);

	g.removeNode(np->getId());

	node* root_np;
	if( g.findNode( rbt.getRootId(), root_np) == false)
	{
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	}
	root_np->getAdjNodeStrings(vs);
	if( vs.length() != 2)
	{
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	}
	edge *ep1, *ep2;
	if( g.findEdge( root_np->getId(), vs[0], ep1 ) == false )
	{
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	}
	if( g.findEdge( root_np->getId(), vs[1], ep2 ) == false )
	{
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	}


	phyloEdgeObject* e1_dop = dynamic_cast<phyloEdgeObject*>(ep1->getDOP().ptr());
	phyloEdgeObject* e2_dop = dynamic_cast<phyloEdgeObject*>(ep2->getDOP().ptr());

	ep = g.addEdge(vs[0], vs[1], new phyloEdgeObject( edge::edgeString(vs[0], vs[1]), "", "", "", e1_dop->mylen + e2_dop->mylen ) );
	g.removeNode(rbt.getRootId());

	rootedbtree samp_tree(g, new_root->getId());
	phyloGraphObject* gdop = dynamic_cast<phyloGraphObject*>( rbt.getDOP().ptr() );
	vs.resize( gdop->tip_strings.length() - 1);
	for(int i=0; i<vs.length(); i++)
	{
		vs[i] = gdop->tip_strings[i+1];
	}
	samp_tree.setDOP( new phyloGraphObject( "sterling_tree", vs ) );

	/*
	cout<<"orig tree"<<endl;
	print_phylotree(rbt);
	*/

	precompute(&samp_tree, rate, stat_dist);
	/*
	cout<<"samp_tree"<<endl;
	print_phylotree(samp_tree);
	*/

	gslmat tmp;
	//cout<<init_dist<<endl;
	sampletips_noprecompute(samp_tree, rate, init_dist, nsamp, seed, tmp);

	output.setsize(nsamp, tmp.ncol() + 1);
	for(int r=0; r<nsamp; r++)
	{
		for(int c=0; c<tmp.ncol() + 1; c++)
		{
			double d;
			if(c == 0)
				d = 2.0;
			else
				d = gsl_matrix_get(tmp, r, c-1);

			gsl_matrix_set(output,r,c,d);
		}
	}

	return;
}
