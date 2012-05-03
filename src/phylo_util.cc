/*
 * phylo_util.cc
 *
 *  Created on: May 14, 2010
 *      Author: dnlennon
 */


#include "gslmat.h"
#include "gslvec.h"

#include "graph.h"
#include "graphalg.h"

#include "rootedbtree.h"
#include "unrootedbtree.h"

#include "phylo_objs.h"
#include "phylo_util.h"
#include "phylo_rbtalg.h"

#include "queue.h"

#include "RecoverableAssertException.h"

#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <string>
using namespace std;


#include <stdlib.h>

#include <string>
using namespace std;

void splitString(const string& s0, char split,
		string& s1, /* output */
		string& s2  /* output */
		)
{
	string s = s0;

	int sloc = s.find_first_of(split);

	if(sloc == -1)
	{
		s1 = s0;
		s2 = "";
	}
	else
	{
		s1 = s.substr(0, sloc);
		s2 = s.substr(sloc+1, s.length() - sloc);
	}
}

void newickPairSplice(const string& ps0,
		string& n1, /* output */
		double& d1, /* output */
		string& n2, /* output */
		double& d2  /* output */
		)
{
	string ps = ps0;

	string si1, si2;
	splitString(ps, ',', si1, si2);

	string sd1, sd2;
	splitString(si1, ':', n1, sd1);
	splitString(si2, ':', n2, sd2);

	d1 = atof(sd1.c_str());
	d2 = atof(sd2.c_str());
}

void tipsFromNewickString(const string& s0, vec<string>& vst)
{
	string s = s0;

	// int pctr = 0;
	string pstr;

	while(s.length() > 0)
	{
		int lp_idx = 0;
		while(s[lp_idx++] == '(');
		lp_idx--;

		if(lp_idx == 0)
			break;

		int rp_idx = lp_idx;
		while(s[1+rp_idx] != ')')
		{
			if(s[1+rp_idx] == '(')
				lp_idx = 2+rp_idx;
			rp_idx++;
		}

		int len = rp_idx - lp_idx + 1;

		string tbs = s.substr(lp_idx, len);

		string n1, n2;
		double d1, d2;
		newickPairSplice(tbs, n1, d1, n2, d2);

		pstr = n1 + "-" + n2;
		s.replace(lp_idx-1, len+2, pstr);
	}
	s = s.substr(0, s.length()-1);

	queue<string> sq;
	string tl, rest;
	while(s.length() > 0)
	{
		splitString(s, '-', tl, rest);
		s = rest;
		sq.insert(tl);
	}

	vst.assign(sq);

	return;
}

void treeFromNewickString(const string& s0, rootedbtree& rbt  /* output */)
{
	string s = s0;

	int pctr = 0;
	char pstr[10];

	graph g;

	while(s.length() > 0)
	{
		int lp_idx = 0;
		while(s[lp_idx++] == '(');
		lp_idx--;

		if(lp_idx == 0)
			break;

		int rp_idx = lp_idx;
		while(s[1+rp_idx] != ')')
		{
			if(s[1+rp_idx] == '(')
				lp_idx = 2+rp_idx;
			rp_idx++;
		}

		int len = rp_idx - lp_idx + 1;

		sprintf(pstr, "i%d", pctr++);

		string tbs = s.substr(lp_idx, len);

		//cout<<tbs<<" --> "<<pstr<<endl;
		//cout.flush();

		string n1, n2;
		double d1, d2;
		newickPairSplice(tbs, n1, d1, n2, d2);

		node *tmp;
		if( g.findNode(n1, tmp) == false )
			g.addNode(n1,   new phyloNodeObject(n1, "", "") );
		if( g.findNode(n2, tmp) == false )
			g.addNode(n2,   new phyloNodeObject(n2, "", "") );
		g.addNode(pstr, new phyloNodeObject(pstr, n1, n2) );

		g.addEdge(n1, pstr, new phyloEdgeObject( edge::edgeString(n1, pstr) , n1, n2, pstr, d1 ) );
		g.addEdge(n2, pstr, new phyloEdgeObject( edge::edgeString(n2, pstr) , n2, n1, pstr, d2 ) );

		s.replace(lp_idx-1, len+2, pstr);
	}
	//cout<<g<<endl;

	vec<string> tlv;
	tipsFromNewickString(s0, tlv);

/*
 *
 	cout<<"Note:  tip data must be indexed as follows:"<<endl;
	for(int i=0; i<tlv.length(); i++)
	{
		cout<<tlv[i];
		if(i < tlv.length() -1)
			cout<<", ";
		else
			cout<<endl;
	}
	*/

	g.setDOP( new phyloGraphObject("global data", tlv) );

	rbt.assign(g, pstr);
}
//////////////////////////////////////////////////////////////////////////////////////

/* Denote the rate matrix by Q.  We assume that the process is time reversible.
 * Then, in general, Q = B D where B is symmetric and D is diagonal.  See Yang,
 * "Computational Molecular Evolution," Section 1.5.
 *
 * In fact, we can say more.  Q is similar to a symmetric matrix.  Let Pi be the
 * stationary distribution.  Then diag(sqrt(Pi)) * Q * diag(1/sqrt(Pi)) is symmetric.
 * Compute the eigen-decomposition of this product (U Lambda U^T) and then transform
 * to obtain the eigen-decomposition of Q.
 *
 * In the above, we need to compute Pi, which satisfies Q^T Pi = 0.   This is done
 * via an SVD characterization of the null space of Q^T.
 */

void diagonalize(
		gslmat &Q  /* rate matrix (input) */,
		gslmat &U_result, /* output variables */
		gslvec &D_result,
		gslmat &Uinv_result,
		gslvec &Pi_result
		)
{
	int n = Q.nrow();

	/* clear output matrices and size them appropriately */
	U_result.setsize(n,n);
	D_result.setsize(n);
	Uinv_result.setsize(n,n);
	Pi_result.setsize(n);

	/* matrices for intermediate calculation */
	gslmat U(n,n);
	gslvec S(n);
	gslmat V(n,n);
	gslvec Pi(n);
	gslvec work(n);
	gslmat Q_trans(n,n);

	for(int r=0; r<n; r++)
	{
		for(int c=0; c<n; c++)
		{
			double Q_rc = gsl_matrix_get(Q, r, c);
			gsl_matrix_set(Q_trans, c, r, Q_rc);
		}
	}
	gsl_matrix_memcpy(U, Q_trans);

	/* compute the SVD and determine the stationary distribution */
	gsl_linalg_SV_decomp(U, V, S, work);
	gsl_matrix_get_col(Pi, V, n-1);
	double scale_factor =0 ;
	for(int i=0; i<n; i++)
		scale_factor += gsl_vector_get(Pi, i);
	gsl_vector_scale(Pi, 1.0/scale_factor);
	gsl_vector_memcpy(Pi_result, Pi);

	/* compute product:  diag(sqrt(Pi)) * Q * diag(1/sqrt(Pi)) */
	gslmat Ds_Q_Dsinv(n,n);
	for(int r=0; r<n; r++)
	{
		double Pi_r = gsl_vector_get(Pi, r);
		for(int c=0; c<n; c++)
		{
			double Pi_c = gsl_vector_get(Pi, c);
			double Q_rc = gsl_matrix_get(Q, r, c);
			gsl_matrix_set(Ds_Q_Dsinv, r, c, sqrt(Pi_r) * Q_rc / sqrt(Pi_c));
		}
	}

	/* compute symmetric eigen-decomposition of Ds_Q_Dsinv */
	gsl_eigen_symmv_workspace *eigen_work = gsl_eigen_symmv_alloc(n);
	gsl_eigen_symmv(Ds_Q_Dsinv, S, U, eigen_work);
	gsl_eigen_symmv_sort(S, U, GSL_EIGEN_SORT_VAL_DESC);

	gsl_vector_memcpy(D_result, S);
	for(int r=0; r<n; r++)
		{
		double Pi_r = gsl_vector_get(Pi, r);
		for(int c=0; c<n; c++)
			{
				double U_rc = gsl_matrix_get(U, r, c);
				gsl_matrix_set(U_result, r, c, U_rc * 1.0/sqrt(Pi_r));
				gsl_matrix_set(Uinv_result, c, r, U_rc * sqrt(Pi_r));
			}
		}

	gsl_eigen_symmv_free(eigen_work);
	return;
}

void Compute_I_tilde(
		gslmat& It,  /* output */
		gslvec& eval,
		double t
		)
{
	const double tol=1e-8;
	int n = It.nrow();
	for(int r=0; r<n; r++)
	{
		for(int c=0; c<n; c++)
		{
			double eval_r = gsl_vector_get(eval, r);
			double eval_c = gsl_vector_get(eval, c);
			if(r==c || fabs(eval_r - eval_c) < tol)
			{
				gsl_matrix_set(It, r, c, t * exp(eval_r * t) );
			}
			else
			{
				double dbl = ( exp( eval_r * t ) - exp( eval_c * t) ) / (eval_r - eval_c );
				gsl_matrix_set(It, r, c, dbl);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////

class precomputeFn : public graphTraversalFunction
{
public:

	precomputeFn(
			gslmat& Q_input,
			gslmat& U_input,
			gslvec& D_input,
			gslmat& Uinv_input,
			gslmat& Q_label_input,
			const_gslvec& dwell_mask_input
			) :
				Q(Q_input),
				U(U_input),
				D(D_input),
				Uinv(Uinv_input),
				Q_label(Q_label_input),
				dwell_mask(dwell_mask_input)
				{}


	virtual void visit(edge* e, node* prev)
	{
		int n = Q.nrow();

		phyloNodeObject* ndop;
		ndop = dynamic_cast<phyloNodeObject*>(prev->getDOP().ptr());
		if( ndop->fval.length() == 0 )
		{
			ndop->fval.setsize(n);
			ndop->gval.setsize(n);
		}

		node* np = e->otherNode(prev);
		ndop = dynamic_cast<phyloNodeObject*>(np->getDOP().ptr());
		ndop->fval.setsize(n);
		ndop->gval.setsize(n);

		UphyloEdgeObject* udop = dynamic_cast<UphyloEdgeObject*>(e->getDOP().ptr());

		/* compute conditional probability matrix (e->cprob) */
		udop->CP_mat.setsize(n,n);
		gslvec expDt(n);
		gslmat Ptmp(n,n);
		gsl_matrix_memcpy(Ptmp, U);
		for(int r=0; r<n; r++)
		{
			for(int c=0; c<n; c++)
			{
				double dbl = exp( gsl_vector_get(D, c) * udop->mylen );
				double U_rc = gsl_matrix_get(U, r, c);
				gsl_matrix_set(Ptmp, r, c, U_rc * dbl);
			}
		}
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Ptmp, Uinv, 0.0, udop->CP_mat);

		phyloEdgeObject* edop= dynamic_cast<phyloEdgeObject*>(e->getDOP().ptr());
		if(edop == NULL)
			return;

		edop->sval.setsize(n);

		/* compute conditional expected labeled jumps matrix (e->jump) */
		edop->EJ_mat.setsize(n,n);
		gslmat I_tilde(n,n);
		gslmat Etmp(n,n);
		gslmat Etmp2(n,n);
		Compute_I_tilde(I_tilde, D, edop->mylen);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Uinv, Q_label, 0.0, Etmp);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Etmp, U, 0.0, Etmp2);
		gsl_matrix_mul_elements(Etmp2, I_tilde);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, Etmp2, 0.0, Etmp);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Etmp, Uinv, 0.0, Etmp2);
		gsl_matrix_memcpy(edop->EJ_mat, Etmp2);

		/* compute conditional expected dwelling time (e->edwell) */
		edop->ED_mat.setsize(n,n);
		for(int r=0; r<n; r++)
		{
			for(int c=0; c<n; c++)
			{
				double dm_c = gsl_vector_get(dwell_mask, c);
				double Uinv_rc = gsl_matrix_get(Uinv, r, c);
				gsl_matrix_set(Etmp, r, c, Uinv_rc * dm_c);
			}
		}
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Etmp, U, 0.0, Etmp2);
		gsl_matrix_mul_elements(Etmp2, I_tilde);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, Etmp2, 0.0, Etmp);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Etmp, Uinv, 0.0, Etmp2);
		gsl_matrix_memcpy(edop->ED_mat, Etmp2);
	}

private:
	gslmat& Q;
	gslmat& U;
	gslvec& D;
	gslmat& Uinv;
	gslmat& Q_label;
	const_gslvec& dwell_mask;
};

void precompute(
		graph* g,
		gslmat& rate,
		gslvec& stat_dist /* output */
		)
{
	int nr = rate.nrow();
	int nc = rate.ncol();
	gslmat label_mask(nr, nc);
	gslvec dwell_mask(nr);

	precompute(g, rate, label_mask, dwell_mask, stat_dist);
}


void precompute(
		graph* g,
		gslmat& rate,
		gslmat& label_mask,
		const_gslvec dwell_mask,
		gslvec& stat_dist /* output */
		)
{
	int n = rate.nrow();
	gslmat U, Uinv;
	gslvec D;

	diagonalize(rate, U, D, Uinv, stat_dist);

	gslmat Q_label(n,n);
	gsl_matrix_memcpy(Q_label, rate);
	gsl_matrix_mul_elements(Q_label, label_mask);

	precomputeFn pf(rate, U, D, Uinv, Q_label, dwell_mask);

	vec<node*> nvp;
	g->nodeList(nvp);
	bfs(*g, nvp[0]->getId(), &pf);

	return;
}

//////////////////////////////////////////////////////////////////////////////////////

class fvalFn : public graphTraversalFunction
{
public:
	virtual void visit(edge* e, node* prev)
	{
		node* ch = e->otherNode(prev);

		edge *el = NULL;
		edge *er = NULL;
		vec<string> ins;
		ch->getAdjNodeStrings(ins);

		//UphyloNodeObject* pn_udop = dynamic_cast<UphyloNodeObject*>(prev->getDOP().ptr());
		UphyloNodeObject* cn_udop = dynamic_cast<UphyloNodeObject*>(ch->getDOP().ptr());


		if(ins.length() == 1)
			return;

		for(int i=0; i<ins.length(); i++)
		{
			edge* ep = ch->getep(ins[i]);
			if(ep == e)
			{
				continue;
			}
			else if(el == NULL)
			{
				el = ep;
			}
			else
			{
				er = ep;
			}
		}

		if(el == NULL || er == NULL)
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

		node* gcl = el->otherNode(ch);
		node* gcr = er->otherNode(ch);
		UphyloNodeObject* gcl_udop = dynamic_cast<UphyloNodeObject*>(gcl->getDOP().ptr());
		UphyloNodeObject* gcr_udop = dynamic_cast<UphyloNodeObject*>(gcr->getDOP().ptr());

		UphyloEdgeObject* el_udop = dynamic_cast<UphyloEdgeObject*>(el->getDOP().ptr());
		UphyloEdgeObject* er_udop = dynamic_cast<UphyloEdgeObject*>(er->getDOP().ptr());

		gslvec sv_sl(el_udop->CP_mat.ncol());
		gslvec sv_sr(el_udop->CP_mat.ncol());
		gsl_blas_dgemv(CblasNoTrans, 1.0, el_udop->CP_mat, gcl_udop->fval, 0.0, sv_sl);
		gsl_blas_dgemv(CblasNoTrans, 1.0, er_udop->CP_mat, gcr_udop->fval, 0.0, sv_sr);

		gsl_vector_memcpy(cn_udop->fval, sv_sl);
		gsl_vector_mul(cn_udop->fval, sv_sr);

		phyloEdgeObject* sl_pdop = dynamic_cast<phyloEdgeObject*>(el->getDOP().ptr());
		phyloEdgeObject* sr_pdop = dynamic_cast<phyloEdgeObject*>(er->getDOP().ptr());
		if(sl_pdop != NULL && sr_pdop != NULL)
		{
			gsl_vector_memcpy(sl_pdop->sval, sv_sl);
			gsl_vector_memcpy(sr_pdop->sval, sv_sr);
		}
	}

	virtual void root_patchup(node* r)
	{
		vec<string> vs;
		r->getAdjNodeStrings(vs);
		int nedge = vs.length();

		phyloNodeObject* rdop = dynamic_cast<phyloNodeObject*>(r->getDOP().ptr());
		if(rdop == NULL)
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

		edge** e = new edge*[nedge];
		node** ch = new node*[nedge];
		UphyloEdgeObject** e_dop = new UphyloEdgeObject*[nedge];
		UphyloNodeObject** n_dop = new UphyloNodeObject*[nedge];

		for(int i=0; i<nedge; i++)
		{
			e[i] = r->getep(vs[i]);
			e_dop[i] = dynamic_cast<UphyloEdgeObject*>(e[i]->getDOP().ptr());

			ch[i] = e[i]->otherNode(r);
			n_dop[i] = dynamic_cast<UphyloNodeObject*>(ch[i]->getDOP().ptr());

			gslvec rdop_orig_fval(e_dop[i]->CP_mat.ncol());
			gsl_vector_memcpy(rdop_orig_fval, rdop->fval);

			gslvec tmp(e_dop[i]->CP_mat.ncol());
			gsl_blas_dgemv(CblasNoTrans, 1.0, e_dop[i]->CP_mat, n_dop[i]->fval, 0.0, tmp);

			if(i == 0)
				gsl_vector_memcpy(rdop->fval, tmp);
			else
				gsl_vector_mul(rdop->fval, tmp);

			if(nedge == 1)
			{
				gsl_vector_mul(rdop->fval, rdop_orig_fval);
			}
			else if(nedge == 2)
			{
				phyloEdgeObject* pdop = dynamic_cast<phyloEdgeObject*>(e[i]->getDOP().ptr());
				gsl_vector_memcpy(pdop->sval, tmp);
			}
		}

		delete [] e;
		delete [] ch;
		delete [] e_dop;
		delete [] n_dop;
	}
};

class gvalFn : public graphTraversalFunction
{
public:
	virtual void visit(edge* e, node* prev)
	{
		phyloNodeObject* np_dop = dynamic_cast<phyloNodeObject*>(prev->getDOP().ptr());
		node* ch = e->otherNode(prev);
		phyloNodeObject* nc_dop = dynamic_cast<phyloNodeObject*>(ch->getDOP().ptr());

		phyloEdgeObject* el_dop = dynamic_cast<phyloEdgeObject*>(e->getDOP().ptr());
		edge* er = prev->getep( edge::edgeString( el_dop->sibstr, el_dop->parstr ) );
		phyloEdgeObject* er_dop = dynamic_cast<phyloEdgeObject*>(er->getDOP().ptr());

		//node* sib = er->otherNode(prev);
		// phyloNodeObject* ns_dop = dynamic_cast<phyloNodeObject*>(sib->getDOP().ptr());

		int n=np_dop->fval.length();
		gslvec Gtmp(n);

		gsl_vector_memcpy(Gtmp, np_dop->gval);
		gsl_vector_mul(Gtmp, er_dop->sval);
		gsl_blas_dgemv(CblasTrans, 1.0, el_dop->CP_mat, Gtmp, 0.0, nc_dop->gval);
	}
};

void fgval_vec_compute(
		graph* g,
		string root_id,
		const_gslvec root_dist,
		gslmat& tip_vec /* a vector in {1..n_states}^N_tips */
		)
{
	if(tip_vec.nrow() != 1)
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

	gsl_vector_view gvv = gsl_matrix_row(tip_vec, 0);
	fgval_vec_compute(g, root_id, root_dist, gvv);
}

double fgval_vec_compute(
		graph* g,
		string root_id,
		const_gslvec root_dist,
		gslvec tip_vec /* a vector in {1..n_states}^N_tips */
		)
{
	int nr = root_dist.size();
	int nc = tip_vec.size();

	gslmat tip_vecm(nr, nc);
	for(int c=0; c<nc; c++)
	{
		gsl_vector_view gvv = gsl_matrix_column(tip_vecm, c);
		gsl_vector_set_zero(&gvv.vector);
		int idx = gsl_vector_get(tip_vec, c);
		if(idx < 1 || idx > nc)
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

		idx = idx - 1;
		gsl_vector_set(&gvv.vector, idx, 1.0);
	}

	fgval_01mat_compute(g, root_id, root_dist, tip_vecm);

	node* np;
	g->findNode(root_id, np);
	phyloNodeObject* ndop = dynamic_cast<phyloNodeObject*> (np->getDOP().ptr());
	double lhood;
	gsl_blas_ddot(ndop->fval, root_dist, &lhood);
	return lhood;
}

void fgval_01mat_compute(
		graph* g,
		string root_id,
		const_gslvec root_dist,
		gslmat& tip_vec /* each column is a vector in {0,1}^n_states */
		)
{
	rootedbtree* rbtp = dynamic_cast<rootedbtree*>(g);

	phyloGraphObject* gdop = dynamic_cast<phyloGraphObject*>(g->getDOP().ptr());
	if( gdop == NULL)
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

	// initialize tip vectors
	for(int i=0; i<gdop->tip_strings.length(); i++)
	{
		node* np;
		if( g->findNode(gdop->tip_strings[i], np) == false )
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

		phyloNodeObject* ndop = dynamic_cast<phyloNodeObject*>(np->getDOP().ptr());
		if(ndop == NULL)
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

		gsl_vector_const_view gvv = gsl_matrix_const_column(tip_vec, i);
		gsl_vector_memcpy(ndop->fval, &gvv.vector);
	}

	// compute fval
	fvalFn fvalFn;
	rbfs(*g, root_id, &fvalFn);

	// compute gval
	if(rbtp != NULL)
	{
		node* np;
		if( rbtp->findNode(rbtp->getRootId(), np) == false)
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
		phyloNodeObject* npod = dynamic_cast<phyloNodeObject*>(np->getDOP().ptr());
		if(npod == NULL)
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

		gsl_vector_memcpy(npod->gval, root_dist);
		gvalFn gvalFn;
		bfs(*rbtp, root_id, &gvalFn);
	}
}

//////////////////////////////////////////////////////////////////////////////////////

class print_phyloFn : public graphTraversalFunction
{
public:
	virtual void visit(edge* e, node* prev)
	{
		node* np = e->otherNode(prev);

		cout<<"\n"<<np->getId()<<" :------------------------------"<<endl;
		phyloNodeObject* ndop = dynamic_cast<phyloNodeObject*>( np->getDOP().ptr() );

		double lhood;
		gsl_blas_ddot(ndop->fval, ndop->gval, &lhood);

		cout<<"\t fval = "<<ndop->fval<<endl;
		cout<<"\t gval = "<<ndop->gval<<endl;
		cout<<"\t lhood = "<<lhood<<endl;
		cout<<endl;

		cout<<"\n"<<e->getId()<<" :------------------------------"<<endl;
		phyloEdgeObject* edop = dynamic_cast<phyloEdgeObject*>( e->getDOP().ptr() );
		cout<<"\t len = "<<edop->mylen<<endl<<endl;
		cout<<"\t lhood_branch = "<<edop->lhood<<endl;
		cout<<"\t ej_branch = "<<edop->ej_branch<<endl;
		cout<<"\t ed_branch = "<<edop->ed_branch<<endl<<endl;
		cout<<"\t CP = \n"<<edop->CP_mat<<endl<<endl;
		cout<<"\t EJ = \n"<<edop->EJ_mat<<endl<<endl;
		cout<<"\t ED = \n"<<edop->ED_mat<<endl<<endl;
		cout<<endl;
	}
};

void print_phylotree(rootedbtree& rbt)
{
	print_phyloFn print_phyloFn;
	bfs(rbt, rbt.getRootId(), &print_phyloFn);

	string root = rbt.getRootId();
	node* np;
	rbt.findNode(root, np);

	cout<<"\n"<<np->getId()<<" :------------------------------"<<endl;
	phyloNodeObject* ndop = dynamic_cast<phyloNodeObject*>( np->getDOP().ptr() );

	double lhood;
	gsl_blas_ddot(ndop->fval, ndop->gval, &lhood);

	cout<<"\t fval = "<<ndop->fval<<endl;
	cout<<"\t gval = "<<ndop->gval<<endl;
	cout<<"\t lhood = "<<lhood<<endl;
	cout<<endl;

}

//////////////////////////////////////////////////////////////////////////////////////

class stochmap_Fn : public graphTraversalFunction
{
public:
	virtual void visit(edge* e, node* prev)
	{
		phyloNodeObject* par_dop = dynamic_cast<phyloNodeObject*>(prev->getDOP().ptr());
		node* ch = e->otherNode(prev);
		phyloNodeObject* ch_dop = dynamic_cast<phyloNodeObject*>(ch->getDOP().ptr());

		phyloEdgeObject* ecur_dop = dynamic_cast<phyloEdgeObject*>(e->getDOP().ptr());
		edge* er = prev->getep( edge::edgeString( ecur_dop->sibstr, ecur_dop->parstr ) );
		phyloEdgeObject* eopp_dop = dynamic_cast<phyloEdgeObject*>(er->getDOP().ptr());

		// node* sib = er->otherNode(prev);
		//phyloNodeObject* nsib_dop = dynamic_cast<phyloNodeObject*>(sib->getDOP().ptr());

		int n = par_dop->fval.size();
		gslvec GS(n);
		gslvec Vtmp(n);
		gsl_vector_memcpy(GS, par_dop->gval);
		gsl_vector_mul(GS, eopp_dop->sval);

		/*
		cout<<endl<<e->getId()<<endl;
		cout<<"fval:\t"<<par_dop->fval<<endl;
		cout<<"gval:\t"<<par_dop->gval<<endl;
		cout<<"sval:\t"<<eopp_dop->sval<<endl;
		*/

		gsl_blas_dgemv(CblasNoTrans, 1.0, ecur_dop->CP_mat, ch_dop->fval, 0.0, Vtmp);
		gsl_blas_ddot(GS, Vtmp, &(ecur_dop->lhood));

		gsl_blas_dgemv(CblasNoTrans, 1.0, ecur_dop->EJ_mat, ch_dop->fval, 0.0, Vtmp);
		gsl_blas_ddot(GS, Vtmp, &(ecur_dop->ej_branch));
		ecur_dop->ej_branch /= ecur_dop->lhood;

		gsl_blas_dgemv(CblasNoTrans, 1.0, ecur_dop->ED_mat, ch_dop->fval, 0.0, Vtmp);
		gsl_blas_ddot(GS, Vtmp, &(ecur_dop->ed_branch));
		ecur_dop->ed_branch /= ecur_dop->lhood;
	}
};

void stochmapcompute(rootedbtree& rbt)
{
	stochmap_Fn stochmap_Fn;
	bfs(rbt, rbt.getRootId(), &stochmap_Fn);
}

//////////////////////////////////////////////////////////////////////////////////////

class smapextract_Fn : public graphTraversalFunction
{
public:
	smapextract_Fn(gslvec& gv) : mygv(gv) { }

	virtual void visit(edge* e, node* prev)
	{
		phyloEdgeObject* edop = dynamic_cast<phyloEdgeObject*>(e->getDOP().ptr());

		lhood = edop->lhood;
		myqj.insert(edop->ej_branch);
		myqd.insert(edop->ed_branch);
	}

	virtual void finalize(void)
	{
		vec<double> vj;
		vec<double> vd;

		vj.assign(myqj);
		vd.assign(myqd);

		gsl_vector_set(mygv, 0, lhood);
		for(int i=1; i<=vj.length(); i++)
		{
			gsl_vector_set(mygv, i, vj[i-1]);
			gsl_vector_set(mygv, i + vd.length(), vd[i-1]);
		}
	}

	gslvec& mygv;
	double lhood;
	queue<double> myqj;
	queue<double> myqd;
};

void smapextract(
		rootedbtree& rbt,
		gslvec result /* output */
		)
{
	int ne = 1 + 2*rbt.numEdges();

	if( result.length() != ne )
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

	smapextract_Fn smapextract_Fn(result);
	phylo_rbt_postorder(rbt, &smapextract_Fn);
}

//////////////////////////////////////////////////////////////////////////////////////

void unroot_phylorbt(rootedbtree& rbt, unrootedbtree& ut /* output */)
{
	vec<graph> vg;
	rbt.connectedComponents(vg);

	string root_string = rbt.getRootId();
	node *np; //, *lp, *rp;
	if( vg[0].findNode(root_string, np) == false )
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

	vec<string> vs;
	np->getAdjNodeStrings(vs);

	if( vs.length() != 2)
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

	edge *ep1, *ep2;
	if( rbt.findEdge(vs[0], root_string, ep1) == false )
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	if( rbt.findEdge(vs[1], root_string, ep2) == false )
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

	phyloEdgeObject* e1_pod = dynamic_cast<phyloEdgeObject*>(ep1->getDOP().ptr());
	phyloEdgeObject* e2_pod = dynamic_cast<phyloEdgeObject*>(ep2->getDOP().ptr());
	double newlen = e1_pod->mylen + e2_pod->mylen;

	vg[0].removeNode(root_string);
	vg[0].addEdge(vs[0], vs[1], new UphyloEdgeObject( edge::edgeString(vs[0], vs[1]), vs[0], vs[1], newlen ));

	ut = vg[0];
	vec<edge*> vep;
	vg[0].edgeList(vep);
	for(int i=0; i<vep.length(); i++)
	{
		phyloEdgeObject* pdop = dynamic_cast<phyloEdgeObject*>(vep[i]->getDOP().ptr());
		if(pdop == NULL)
			continue;

		edge* ep;
		if( ut.findEdge(pdop->parstr, pdop->chstr, ep) == false)
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

		ep->setDOP( new UphyloEdgeObject(pdop->id(), pdop->parstr, pdop->chstr, pdop->mylen) );
	}

	ut.setDOP( rbt.getDOP() );
}
