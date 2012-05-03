/*
 * gslvec.cc
 *
 *  Created on: Feb 4, 2010
 *      Author: dnlennon
 */

#include "gslvec.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex_double.h>

#include <string>
#include <fstream>

#include <assert.h>

using namespace std;

gslvec::gslvec(istream& fin)
{
	int n, ncol, nrow;
	double dbl;

	fin.ignore();
	fin>>nrow>>ncol;
	if(nrow == 1)
		n = ncol;
	else if(ncol == 1)
		n = nrow;
	else
		assert(0);

	mygsl = NULL;
	setsize(n);

	for(int r=0; r<n; r++)
	{
		assert(fin.eof() == false);
		fin>>dbl;
		fin.ignore();
		gsl_vector_set(mygsl, r, dbl);
	}
	myview = false;
}

gslvec::gslvec(gsl_vector_view& gvv)
{
	myview = true;
	mygsl = &gvv.vector;
}

gslvec::gslvec(gsl_vector& gv)
{
	myview = true;
	mygsl = &gv;
}

/*
gslvec::gslvec(gslvec& gv)
{
	myview = true;
	mygsl = gv.mygsl;
}
*/

gslvec::~gslvec()
{
	if(mygsl != NULL && myview == false)
		gsl_vector_free(mygsl);
}

void gslvec::setsize(int n)
{
	if(mygsl != NULL && myview == false)
		gsl_vector_free(mygsl);

	if( n == 0 )
		return;

	mygsl = gsl_vector_alloc(n);
	gsl_vector_set_zero(mygsl);
}

/////////////////////////////////////////////////////////////////////////////


ostream& operator<<(ostream& o, const gslvec& m)
{
	int len = m.length();

	double dbl;
	o.setf(ios_base::fixed, ios_base::floatfield);
	for(int r=0; r<len; r++)
	{
		dbl = gsl_vector_get(m, r);
		o.width(12);
		o.precision(8);
		o<<dbl;
	}
	o<<endl;
	o.setf(std::_Ios_Fmtflags(0), ios_base::floatfield);


	return o;
}
