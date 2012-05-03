/*
 * rootedtree.h
 *
 *  Created on: May 13, 2010
 *      Author: dnlennon
 */

#ifndef GENERALTREE_H_
#define GENERALTREE_H_

#include "graph.h"

class generaltree : public graph
{
public:
	generaltree() {}

	generaltree(const graph& g) : graph(g) { sanityCheckWithException(true); }
	generaltree(const generaltree& t) : graph(t) { sanityCheckWithException(true); }

	generaltree& operator=(const graph& t)
	{
		if(&t == this)
			return *this;

		this->graph::operator=(t);
		sanityCheckWithException(false);
		return *this;
	}

	generaltree& operator=(const generaltree& t)
	{
		if(&t == this)
			return *this;

		this->graph::operator=(t);
		sanityCheckWithException(false);
		return *this;
	}

	virtual bool checkSanity(bool incremental);
};

#endif /* GENERALTREE_H_ */
