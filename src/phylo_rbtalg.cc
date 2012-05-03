/*
 * phylo_rbtalg.cc
 *
 *  Created on: May 31, 2010
 *      Author: dnlennon
 */

#include "rootedbtree.h"
#include "graphalg.h"
#include "graph.h"

#include "RecoverableAssertException.h"

#include "phylo_objs.h"
#include "phylo_rbtalg.h"

void phylo_rbt_postorder_recurse(node* np, graphTraversalFunction* gtf)
{
	phyloNodeObject* npod = dynamic_cast<phyloNodeObject*>(np->getDOP().ptr());

	edge* ep = np->getep(npod->left);
	if(ep != NULL)
	{
		phylo_rbt_postorder_recurse(ep->otherNode(np), gtf);
		gtf->visit(ep, np);
	}

	ep = np->getep(npod->right);
	if(ep != NULL)
	{
		phylo_rbt_postorder_recurse(ep->otherNode(np), gtf);
		gtf->visit(ep, np);
	}
}

void phylo_rbt_postorder(rootedbtree& rbt, graphTraversalFunction* gtf)
{
	node* np;
	if( rbt.findNode(rbt.getRootId(), np) == false)
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

	if( np->numEdges() != 2)
		throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);

	phylo_rbt_postorder_recurse(np, gtf);

	gtf->finalize();
}
