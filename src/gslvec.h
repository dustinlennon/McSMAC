/*
 * gslvec.h
 *
 *  Created on: Feb 4, 2010
 *      Author: dnlennon
 */

#ifndef GSLVEC_H_
#define GSLVEC_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex_double.h>

#include <string>
#include <fstream>

using namespace std;

typedef class gslvec gslvec;
typedef class const_gslvec const_gslvec;

class gslvec
{
	friend class const_gslvec;
public:
	gslvec()
	{
		mygsl = NULL;
		myview = false;
	}

	gslvec(int n)
	{
		mygsl=NULL;
		myview = false;
		setsize(n);
	}

	gslvec(const gslvec& gv) : mygsl(gv.mygsl) { myview = true; }
	gslvec(gsl_vector_view& gvv);
	gslvec(gsl_vector& gv);

	gslvec(istream& fin);

	~gslvec();

	operator gsl_vector*() { return mygsl; }
	operator const gsl_vector*() const { return mygsl; }

	int length(void) const
	{
		if(mygsl == NULL)
			return 0;
		return mygsl->size;
	}

	int size(void) const
	{
		return length();
	}

	void setsize(int n);

	double* getDoublePtr(void) { return mygsl->data; }
	const double* getConstDoublePtr(void) const { return mygsl->data; }

	bool isGslNull(void) const { return mygsl == NULL; }

private:
	gsl_vector* mygsl;
	bool myview;
};


class const_gslvec
{
public:
	const_gslvec(gsl_vector_const_view& gvv) : mygsl(&gvv.vector) { myview = true; }
	const_gslvec(gsl_vector_view& gvv) : mygsl(&gvv.vector) { myview = true; }
	const_gslvec(const gsl_vector& gv) : mygsl(&gv) { myview = true; }
	const_gslvec(gsl_vector& gv) : mygsl(&gv) { myview = true; }

	const_gslvec(const const_gslvec& gv) : mygsl(gv) { myview = true; }
	const_gslvec(const gslvec& gv) : mygsl(gv) { myview = true; }

	operator const gsl_vector*() const { return mygsl; }

	~const_gslvec()
	{
		if(myview == false)
			delete mygsl;
	}

	const_gslvec(istream& fin)
	{
		myview = false;
		gslvec gv(fin);
		mygsl = gv.mygsl;
		gv.myview = true;
	}

	int length(void) const
	{
		if(mygsl == NULL)
			return 0;
		return mygsl->size;
	}

	int size(void) const
	{
		return length();
	}

	const double* getConstDoublePtr(void) const { return mygsl->data; }

private:
	const gsl_vector* mygsl;
	bool myview;
};

ostream& operator<<(ostream& o, const gslvec& v);

#endif /* GSLVEC_H_ */
