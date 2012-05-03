/*
 * graphalg.cc
 *
 *  Created on: May 16, 2010
 *      Author: dnlennon
 */

#include "graph.h"
#include "queue.h"
#include "stack.h"
#include "graphalg.h"

class enpair
{
public:
	enpair(edge* ep0, node* np0) { np = np0;  ep = ep0; }
	node* np;
	edge* ep;
};

bool bfs(graph& g, const string& root, graphTraversalFunction* gtf, bool forward);

bool bfs(graph& g, const string& root, graphTraversalFunction* gtf)
{
	return bfs(g, root, gtf, true);
}

bool rbfs(graph& g, const string& root, graphTraversalFunction* gtf)
{
	return bfs(g, root, gtf, false);
}

bool bfs(graph& g, const string& root, graphTraversalFunction* gtf, bool forward)
{
	node* np_root;
	if( g.findNode(root, np_root) == false)
		return false;

	queuestackBase<enpair> *visit_order;
	if(forward == true)
		visit_order = new queue<enpair>;
	else
		visit_order = new stack<enpair>;

	search_set<node*> ssnp;
	ssnp.insert(np_root);

	queue<node*> tip_nodes;

	queue<node*> nq;
	nq.insert(np_root);
	while(nq.empty() == false)
	{
		node* np = nq.next();
		vec<string> ves;
		np->getAdjNodeStrings(ves);
		for(int i=0; i<ves.length(); i++)
		{
			edge* ep;
			if( g.findEdge(np->getId(), ves[i], ep) == false )
				return false;

			node* next = ep->otherNode(np);
			if( ssnp.inSet(next) == true )
				continue;

			enpair enp(ep, np);
			visit_order->insert( enp );

			ssnp.insert(next);
			nq.insert(next);
		}
		if(ves.length() == 1)
		{
			tip_nodes.insert(np);
		}

		nq.remove();
	}

	while(visit_order->empty() == false)
	{
		enpair tmp = visit_order->next();
		gtf->visit(tmp.ep, tmp.np);

		visit_order->remove();
	}

	if(forward == false)
	{
		gtf->root_patchup(np_root);
	}
	else
	{
		vec<node*> vnp;
		vnp.assign(tip_nodes);
		for(int i=0; i<vnp.length(); i++)
			gtf->tip_patchup(vnp[i]);
	}

	gtf->finalize();

	delete visit_order;
	return true;
}
