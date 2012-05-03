/*
 * unrootedbtree.cc
 *
 *  Created on: May 14, 2010
 *      Author: dnlennon
 */

#include "rootedbtree.h"
#include "RecoverableAssertException.h"

bool rootedbtree::checkSanity(bool incremental)
{
	checkRootNullityWithException();

	if(incremental == false)
	{
		generaltree::checkSanity(false);
	}

	vec<node*> vnp;
	np_set.getVec(vnp);
	for(int i=0; i<vnp.length(); i++)
	{
		int ne = vnp[i]->numEdges();

		if( vnp[i] == myroot )
		{
			if(ne == 2)
				continue;
			else
				throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
		}

		if( ne == 3 || ne == 1)
			continue;
		else
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	}

	return true;
}

void rootedbtree::assignRootWithException(const string& rootname)
{
	if( findNode(rootname, myroot) == false )
	{
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	}
}

void rootedbtree::checkRootNullityWithException(void)
{
	if(myroot == NULL)
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
}


void rootedbtree::assign(const graph& g, const string& rootname)
{
	deepClean();

	this->graph::operator=(g);
	assignRootWithException(rootname);
	sanityCheckWithException(false);
}

void rootedbtree::assign(const generaltree& gt, const string& rootname)
{
	deepClean();

	this->generaltree::operator=(gt);
	assignRootWithException(rootname);
	sanityCheckWithException(false);
}

void rootedbtree::deepClean(void)
{
	graph::deepClean();
	myroot = NULL;
}

/*
bool bfs(rootedbtree& rbt, graphTraversalFunction* gtf)
{
	return bfs(rbt, rbt.getRootId(), gtf);
}

bool rbfs(rootedbtree& rbt, graphTraversalFunction* gtf)
{
	return rbfs(rbt, rbt.getRootId(), gtf);
}
*/
