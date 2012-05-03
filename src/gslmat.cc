/*
 * gslmat.cc
 *
 *  Created on: Feb 4, 2010
 *      Author: dnlennon
 */

#include "gslmat.h"

#include <string>
#include <fstream>

#include <assert.h>

using namespace std;

gslmat::gslmat(istream& fin)
{
	int nrow, ncol;
	double dbl;

	fin.ignore();
	fin>>nrow>>ncol;
	mygsl = gsl_matrix_alloc(nrow, ncol);
	assert(mygsl != NULL);
	for(int r=0; r<nrow; r++)
	{
		for(int c=0; c<ncol; c++)
		{
			assert(fin.eof() == false);
			fin>>dbl;
			fin.ignore();
			gsl_matrix_set(mygsl, r, c, dbl);
		}
	}
	myview = false;
}

gslmat::gslmat(gsl_matrix_view& gmv)
{
	mygsl = &gmv.matrix;
	myview = true;
}

gslmat::~gslmat()
{
	if(mygsl != NULL && myview == false)
		gsl_matrix_free(mygsl);
}

void gslmat::setsize(int n, int c)
{
	if(mygsl != NULL && myview == false)
		gsl_matrix_free(mygsl);

	if(n == 0 || c == 0)
		return;

	mygsl = gsl_matrix_alloc(n, c);
	gsl_matrix_set_zero(mygsl);
}

ostream& operator<<(ostream& o, const gslmat& m)
{
	int nrow = m.nrow();
	int ncol = m.ncol();

	double dbl;
	o.setf(ios_base::fixed, ios_base::floatfield);
	for(int r=0; r<nrow; r++)
	{
		for(int c=0; c<ncol; c++)
		{
			dbl = gsl_matrix_get(m, r, c);
			o.width(12);
			o.precision(8);
			o<<dbl;
		}
		o<<endl;
	}
	o.setf(std::_Ios_Fmtflags(0), ios_base::floatfield);

	return o;
}

