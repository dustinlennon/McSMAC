/*
 * R_interface.cc
 *
 *  Created on: Feb 8, 2010
 *      Author: dnlennon
 */

#include "rootedbtree.h"
#include "phylo_util.h"
#include "fastmap.h"
#include "gslmat.h"
#include "gslvec.h"

#include "RecoverableAssertException.h"

#include <assert.h>

#include <Rinternals.h>
#include <gsl/gsl_matrix.h>

#include <math.h>

#include <iostream>
using namespace std;

void ASSIGN_COL_NAMES(SEXP Rmat, vec<string>& vnames)
{
	SEXP cnames, dimnames;

	int n = vnames.size();
	int nc = 1 + 2*n;
	PROTECT( cnames = Rf_allocVector( STRSXP , nc ) );
	SET_STRING_ELT(cnames, 0, mkChar(""));
	for(int i=1; i<=n; i++)
	{
		SET_STRING_ELT(cnames, i, mkChar( vnames[i-1].c_str() ) );
		SET_STRING_ELT(cnames, n+i, mkChar( vnames[i-1].c_str() ) );
	}

	PROTECT( dimnames  = allocVector(VECSXP, 2) );
	SET_VECTOR_ELT(dimnames, 0, R_NilValue);
	SET_VECTOR_ELT(dimnames, 1, cnames);

	setAttrib(Rmat, R_DimNamesSymbol, dimnames);
	UNPROTECT(2);
}

void RBLAS_to_GSL(gslmat& Gmat, SEXP Rmat)
{
	if(! isMatrix(Rmat) )
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	int nr = nrows(Rmat);
	int nc = ncols(Rmat);
	double* base = REAL(Rmat);

	Gmat.setsize(nr, nc);

	for(int r=0; r<nr; r++)
	{
		for(int c=0; c<nc; c++)
		{
			double dbl = base[r + nr*c];
			gsl_matrix_set(Gmat, r, c, dbl);
		}
	}
}

void GSL_to_RBLAS(SEXP Rmat, gslmat& Gmat)
 {
 	assert(isMatrix(Rmat));
 	int nr = Gmat.nrow();
 	int nc = Gmat.ncol();
 	double* base = REAL(Rmat);
 	int Rmat_nr = nrows(Rmat);
	int Rmat_nc = ncols(Rmat);

	if(! (nr==Rmat_nr && nc==Rmat_nc) )
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);


	for(int r=0; r<nr; r++)
	{
		for(int c=0; c<nc; c++)
		{
			double dbl = gsl_matrix_get(Gmat, r, c);
			base[r + nr*c] = dbl;
		}
 	}
 }


extern "C"
{

SEXP WQ_MAT(SEXP m)
{
	gslmat rate;
	RBLAS_to_GSL(rate, m);
	cout<<rate<<endl;

	SEXP ans;
	PROTECT( ans = Rf_allocMatrix(REALSXP, rate.nrow(), rate.ncol()) );
	GSL_to_RBLAS(ans, rate);
	UNPROTECT(1);
	return ans;
}

SEXP MultiMap_RAPI(
		SEXP newick_string_R,
		SEXP rate_R,
		SEXP label_mask_R,
		SEXP dwell_mask_R,
		SEXP site_data_R,
		SEXP root_dist_R,
		SEXP lhood_only_R
		)
{
	assert(isString(newick_string_R));

	string newick_string( CHAR(STRING_ELT(newick_string_R, 0)) );
	rootedbtree t;

	try
	{
		treeFromNewickString(newick_string, t);
	}
	catch(RecoverableAssertException &rae)
	{
		cout<<"FastMap internal error: "<<rae.what()<<endl;
	}

	gslmat rate, label_mask, site_data;
	RBLAS_to_GSL(rate, rate_R);
	RBLAS_to_GSL(label_mask, label_mask_R);
	RBLAS_to_GSL(site_data, site_data_R);

	assert(isNumeric(dwell_mask_R));
	gsl_vector_view dwell_mask_v = gsl_vector_view_array(REAL(dwell_mask_R), length(dwell_mask_R));
	gslvec dwell_mask( dwell_mask_v );

	bool lhood_only = Rf_asInteger(lhood_only_R) == 1 ? true : false;

	gslmat output;
	try
	{
		if(!Rf_isNull(root_dist_R))
		{
			gsl_vector_view root_dist_v = gsl_vector_view_array(REAL(root_dist_R), length(root_dist_R));
			gslvec root_dist ( root_dist_v );
			phylo_multimap(t, rate, label_mask, dwell_mask, site_data, lhood_only, &root_dist, output);
		}
		else
		{
			phylo_multimap(t, rate, label_mask, dwell_mask, site_data, lhood_only, NULL, output);
		}
	}
	catch(RecoverableAssertException &rae)
	{
		cout<<"FastMap internal error: "<<rae.what()<<endl;
	}

	SEXP ans;
	PROTECT( ans = Rf_allocMatrix(REALSXP, output.nrow(), output.ncol()) );
	GSL_to_RBLAS(ans, output);

	if(lhood_only == false)
	{
		vec<string> vnames;
		phylo_edgenames(t, vnames);
		ASSIGN_COL_NAMES(ans, vnames);
	}
	UNPROTECT(1);

	return ans;
}

SEXP OneMap_RAPI(
		SEXP newick_string_R,
		SEXP rate_R,
		SEXP label_mask_R,
		SEXP dwell_mask_R,
		SEXP state_data_R
		)
{
	assert(isString(newick_string_R));

	string newick_string( CHAR(STRING_ELT(newick_string_R, 0)) );
	rootedbtree t;

	try
	{
		treeFromNewickString(newick_string, t);
	}
	catch(RecoverableAssertException &rae)
	{
		cout<<"FastMap internal error: "<<rae.what()<<endl;
	}

	gslmat rate, label_mask, state_data;
	RBLAS_to_GSL(rate, rate_R);
	RBLAS_to_GSL(label_mask, label_mask_R);
	RBLAS_to_GSL(state_data, state_data_R);

	assert(isNumeric(dwell_mask_R));
	gsl_vector_view dwell_mask_v = gsl_vector_view_array(REAL(dwell_mask_R), length(dwell_mask_R));
	gslvec dwell_mask( dwell_mask_v );

	gslmat output;

	try
	{
		phylo_onemap(t, rate, label_mask, dwell_mask, state_data, output);
	}
	catch(RecoverableAssertException &rae)
	{
		cout<<"FastMap internal error: "<<rae.what()<<endl;
	}

	SEXP ans;
	PROTECT( ans = Rf_allocMatrix(REALSXP, output.nrow(), output.ncol()) );
	GSL_to_RBLAS(ans, output);
	UNPROTECT(1);

	vec<string> vnames;
	phylo_edgenames(t, vnames);
	ASSIGN_COL_NAMES(ans, vnames);

	return ans;
}


SEXP SimulateTips_RAPI(
		SEXP newick_string_R,
		SEXP rate_R,
		SEXP init_dist_R,
		SEXP nsamp_R,
		SEXP seed_R
		)
{
	assert(isString(newick_string_R));

	string newick_string( CHAR(STRING_ELT(newick_string_R, 0)) );
	rootedbtree t;

	try
	{
		treeFromNewickString(newick_string, t);
	}
	catch(RecoverableAssertException &rae)
	{
		cout<<"FastMap internal error: "<<rae.what()<<endl;
	}

	gslmat rate;
	RBLAS_to_GSL(rate, rate_R);

	assert(isNumeric(init_dist_R));
	gsl_vector_view init_dist_v = gsl_vector_view_array(REAL(init_dist_R), length(init_dist_R));
	gslvec init_dist( init_dist_v );

	int nsamp = asInteger(nsamp_R);
	int seed = asInteger(seed_R);

	gslmat output;

	try
	{
		sampletips(t, rate, init_dist, nsamp, seed, output);
	}
	catch(RecoverableAssertException &rae)
	{
		cout<<"FastMap internal error: "<<rae.what()<<endl;
	}

	SEXP ans;
	PROTECT( ans = Rf_allocMatrix(REALSXP, output.nrow(), output.ncol()) );
	GSL_to_RBLAS(ans, output);
	UNPROTECT(1);
	return ans;
}


SEXP SimulateTips_Sterling_RAPI(
		SEXP newick_string_R,
		SEXP rate_R,
		SEXP nsamp_R,
		SEXP seed_R
		)
{
	assert(isString(newick_string_R));
	string newick_string( CHAR(STRING_ELT(newick_string_R, 0)) );

	rootedbtree t;

	try
	{
		treeFromNewickString(newick_string, t);
	}
	catch(RecoverableAssertException &rae)
	{
		cout<<"FastMap internal error: "<<rae.what()<<endl;
	}

	gslmat rate;
	RBLAS_to_GSL(rate, rate_R);

	gslvec init_dist( 2 );
	gsl_vector_set(init_dist, 1, 1);

	int nsamp = asInteger(nsamp_R);
	int seed = asInteger(seed_R);

	gslmat output;

	try
	{
		sampletips_sterling_ascertainment(t, rate, nsamp, seed, output);
	}
	catch(RecoverableAssertException &rae)
	{
		cout<<"FastMap internal error: "<<rae.what()<<endl;
	}

	SEXP ans;
	PROTECT( ans = Rf_allocMatrix(REALSXP, output.nrow(), output.ncol()) );
	GSL_to_RBLAS(ans, output);
	UNPROTECT(1);
	return ans;
}


}
