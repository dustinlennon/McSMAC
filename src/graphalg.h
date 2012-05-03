/*
 * graphalg.h
 *
 *  Created on: May 16, 2010
 *      Author: dnlennon
 */

#ifndef GRAPHALG_H_
#define GRAPHALG_H_

#include "graph.h"
#include "dataobject.h"

class graphTraversalFunction : public DOPuser
{
public:
	virtual void visit(edge* ep, node* prev) = 0;
	virtual void root_patchup(node* root) {}  /* called for root node after a reverse bfs */
	virtual void tip_patchup(node* tip) {}    /* called for each tip node after a forward bfs */
	virtual void finalize(void) {}

	virtual ~graphTraversalFunction() {} ;
};

bool bfs(graph& g, const string& root, graphTraversalFunction* gtf);
bool rbfs(graph& g, const string& root, graphTraversalFunction* gtf);

#endif /* GRAPHALG_H_ */
