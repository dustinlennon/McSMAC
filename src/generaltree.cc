/*
 * generaltree.cc
 *
 *  Created on: May 13, 2010
 *      Author: dnlennon
 */

#include "graph.h"
#include "generaltree.h"
#include "search_set.h"

#include "RecoverableAssertException.h"

bool generaltree::checkSanity(bool incremental)
{
	if(incremental == false)
	{
		graph::checkSanity(false);
	}

	/* 1.  there should be a single component. */
	vec< graph > vcc;
	connectedComponents(vcc);
	if(vcc.length() > 1)
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

	/* 2.  there should be no cycles. */
	vec<node*> nvec;
	np_set.getVec(nvec);
	node* np_root = nvec[0];

	dictionary<node*, int> d;
	search_set<edge*> ss_ep;

	queue<node*> qnp;
	qnp.insert(np_root);
	while( qnp.empty() == false )
	{
		int ct;
		node* np = qnp.front();
		if( d.lookup(np, ct) == false )
			d.insert(np, 1);
		else
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

		vec<string> vs;
		np->getAdjNodeStrings(vs);
		for(int i=0; i<vs.length(); i++)
		{
			edge* ep;
			if( findEdge(np->getId(), vs[i], ep) == false )
				assert(0);

			if(ss_ep.inSet(ep) == true)
				continue;

			ss_ep.insert(ep);

			node* tmp = ep->otherNode(np);
			qnp.insert(tmp);
		}

		qnp.remove();
	}

	return true;
}
