/*
 * phylo_objs.h
 *
 *  Created on: May 14, 2010
 *      Author: dnlennon
 */

#ifndef PHYLO_OBJS_H_
#define PHYLO_OBJS_H_

#include "dataobject.h"
#include "gslvec.h"
#include "gslmat.h"
#include "vec.h"

#include <string>
using namespace std;

class UphyloNodeObject : public dataObject
{
	friend ostream& operator<<(ostream& o, const UphyloNodeObject& snop);
public:
	UphyloNodeObject(string id) : dataObject(id) { }
	gslvec fval;
	int sim;

protected:
	dataObject* deepCopy()
	{
		UphyloNodeObject* snop = new UphyloNodeObject(myid);
		snop->fval.setsize( fval.length() );
		if(fval.isGslNull() == false)
			gsl_vector_memcpy(snop->fval, fval);

		snop->sim = sim;
		return snop;
	}
};


class phyloNodeObject : public UphyloNodeObject
{
	friend ostream& operator<<(ostream& o, const phyloNodeObject& snop);
public:
	phyloNodeObject(string id, string l, string r) : UphyloNodeObject(id)
	{
		left = l;
		right = r;
	}
	gslvec gval;

	string left, right;
//	int sim;

private:
	dataObject* deepCopy()
	{
		phyloNodeObject* snop = new phyloNodeObject(myid, left, right);
		snop->fval.setsize( fval.length() );
		if(fval.isGslNull() == false)
			gsl_vector_memcpy(snop->fval, fval);

		snop->gval.setsize( gval.length() );
		if(gval.isGslNull() == false)
			gsl_vector_memcpy(snop->gval, gval);

//		snop->sim = sim;
		return snop;
	}
};

class UphyloEdgeObject : public dataObject
{
	friend ostream& operator<<(ostream& o, const UphyloEdgeObject& snop);
public:
	UphyloEdgeObject(
			string id,
			string first,
			string second,
			double len) : dataObject(id)
	{
		mylen = len;
		n1_str = first;
		n2_str = second;
	}

	double mylen;
	string n1_str, n2_str;
	gslmat CP_mat;

	double lhood;

protected:
	virtual dataObject* deepCopy()
	{
		UphyloEdgeObject* snop = new UphyloEdgeObject(myid, n1_str, n2_str, mylen);

		snop->CP_mat.setsize( CP_mat.nrow(), CP_mat.ncol() );
		if( CP_mat.isGslNull() == false )
			gsl_matrix_memcpy(snop->CP_mat, CP_mat);

		return snop;
	}
};


class phyloEdgeObject : public UphyloEdgeObject
{
	friend ostream& operator<<(ostream& o, const phyloEdgeObject& snop);
public:
	phyloEdgeObject(
			string id,
			string child,
			string sib,
			string parent,
			double len
			) :
				UphyloEdgeObject(id, child, parent, len),
				chstr(n1_str),
				parstr(n2_str)
	{
		chstr = n1_str;
		parstr = n2_str;
		sibstr = sib;
	}

	string &chstr, &parstr;
	string sibstr;
	gslmat EJ_mat, ED_mat;
	gslvec sval;

	double lhood, ej_branch, ed_branch;

protected:
	virtual dataObject* deepCopy()
	{
		phyloEdgeObject* snop = new phyloEdgeObject(myid, chstr, sibstr, parstr, mylen);

		snop->sval.setsize( sval.length() );

		if( sval.isGslNull() == false )
			gsl_vector_memcpy(snop->sval, sval);

		snop->CP_mat.setsize( CP_mat.nrow(), CP_mat.ncol() );
		if( CP_mat.isGslNull() == false )
			gsl_matrix_memcpy(snop->CP_mat, CP_mat);

		snop->EJ_mat.setsize( EJ_mat.nrow(), EJ_mat.ncol() );
		if( EJ_mat.isGslNull() == false )
			gsl_matrix_memcpy(snop->EJ_mat, EJ_mat);

		snop->ED_mat.setsize( ED_mat.nrow(), ED_mat.ncol() );
		if( ED_mat.isGslNull() == false )
			gsl_matrix_memcpy(snop->ED_mat, ED_mat);

		return snop;
	}
};



class phyloGraphObject : public dataObject
{
public:
	phyloGraphObject(string myid, vec<string>& sv) : dataObject(myid), tip_strings(sv) { }
	vec<string> tip_strings;

private:
	dataObject* deepCopy()
	{
		phyloGraphObject* geop = new phyloGraphObject(myid, tip_strings);
		return geop;
	}
};


#endif /* PHYLO_OBJS_H_ */
