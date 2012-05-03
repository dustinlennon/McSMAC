/*
 * rootedbinarytree.h
 *
 *  Created on: May 14, 2010
 *      Author: dnlennon
 */

#ifndef ROOTEDBTREE_H_
#define ROOTEDBTREE_H_

#include "graph.h"
#include "generaltree.h"
#include "graphalg.h"

/* rootedbtree must satisfy
	 1.  the criteria for a generaltree  AND
	 2.  a root with degree 2 must be specified
	 3.  every other vertex must be an internal vertex of degree 3 or a leaf vertex
*/

class rootedbtree : public generaltree
{
public:
	rootedbtree()
	{
		myroot = NULL;
	}

	rootedbtree(const graph& g, const string& rootname) : generaltree(g)
	{
		assignRootWithException(rootname);
		sanityCheckWithException(true);
	}

	rootedbtree(const generaltree& gt, const string& rootname) : generaltree(gt)
	{
		assignRootWithException(rootname);
		sanityCheckWithException(true);
	}

	rootedbtree(const rootedbtree& rbt, const string& rootname) : generaltree(rbt)
	{
		assignRootWithException(rootname);
		sanityCheckWithException(true);
	}

	rootedbtree& operator=(const rootedbtree& rbt)
	{
		if(&rbt == this)
			return *this;

		this->graph::operator=(rbt);

		assignRootWithException(rbt.myroot->getId());
		sanityCheckWithException(false);

		return *this;
	}

	void assign(const graph& g, const string& rootname);
	void assign(const generaltree& gt, const string& rootname);

	virtual bool checkSanity(bool incremental);
	virtual void deepClean(void);

	string getRootId(void)
	{
		checkRootNullityWithException();
		return myroot->getId();
	}

private:
	void checkRootNullityWithException(void);
	void assignRootWithException(const string& rootname);

	node* myroot;
};

//bool bfs(rootedbtree& rbt, graphTraversalFunction* gtf);
//bool rbfs(rootedbtree& rbt, graphTraversalFunction* gtf);


#endif /* ROOTEDBINARYTREE_H_ */
