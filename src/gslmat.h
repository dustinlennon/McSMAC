/*
 * gslmat.h
 *
 *  Created on: Feb 4, 2010
 *      Author: dnlennon
 */

#ifndef GSLMAT_H_
#define GSLMAT_H_

#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_matrix_complex_double.h>

#include <string>
#include <fstream>

using namespace std;

typedef class gslmat gslmat;
typedef class gslmat_z gslmat_z;

class gslmat
{
public:
	gslmat()
	{
		mygsl = NULL;
		myview = false;
	}

	gslmat(int n, int c)
	{
		mygsl=NULL;
		myview = false;
		setsize(n,c);
	}

	gslmat(istream& fin);
	gslmat(gsl_matrix_view& gmv);
	~gslmat();

	void setsize(int n, int c);

	operator gsl_matrix*() { return mygsl; }
	operator const gsl_matrix*() const { return mygsl; }

	int nrow(void) const
	{
		if(isGslNull())
			return 0;
		else
			return mygsl->size1;
	}

	int ncol(void) const
	{
		if(isGslNull())
			return 0;
		else
			return mygsl->size2;
	}

	bool isGslNull(void) const  { return mygsl == NULL; }

	gslmat& operator=(const gslmat& g)
	{
		if(this == &g)
			return *this;

		if(myview == false && mygsl != NULL)
			gsl_matrix_free(mygsl);

		myview = true;
		mygsl = g.mygsl;
	}


private:
	gslmat(const gslmat&) {}
	gsl_matrix* mygsl;
	bool myview;
};
ostream& operator<<(ostream& o, const gslmat& m);


#endif /* GSLMAT_H_ */
