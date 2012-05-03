/*
 * unrootedbinarytree.h
 *
 *  Created on: May 14, 2010
 *      Author: dnlennon
 */

#ifndef UNROOTEDBTREE_H_
#define UNROOTEDBTREE_H_

#include "generaltree.h"

/* unrootedbtree must satisfy
	 1.  the criteria for a generaltree  AND
	 2.  each vertex must be an internal vertex of degree 3 or a leaf vertex
*/

class unrootedbtree : public generaltree
{
public:
	unrootedbtree() {}

	unrootedbtree(const graph& g) : generaltree(g) { sanityCheckWithException(true); }
	unrootedbtree(const generaltree& gt) : generaltree(gt) { sanityCheckWithException(true); }
	unrootedbtree(const unrootedbtree& ubt) : generaltree(ubt) { sanityCheckWithException(true); }

	unrootedbtree& operator=(const graph& g)
	{
		this->graph::operator=(g);
		sanityCheckWithException(false);
		return *this;
	}

	unrootedbtree& operator=(const generaltree& gt)
	{
		this->graph::operator=(gt);
		sanityCheckWithException(false);
		return *this;
	}

	unrootedbtree& operator=(const unrootedbtree& ubt)
	{
		this->graph::operator=(ubt);
		sanityCheckWithException(false);
		return *this;
	}

	virtual bool checkSanity(bool incremental);
};

#endif /* UNROOTEDBINARYTREE_H_ */
