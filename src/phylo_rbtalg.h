/*
 * phylo_rbtalg.h
 *
 *  Created on: May 31, 2010
 *      Author: dnlennon
 */

#ifndef PHYLO_RBTALG_H_
#define PHYLO_RBTALG_H_

#include "rootedbtree.h"
#include "graphalg.h"

void phylo_rbt_postorder(rootedbtree& rbt, graphTraversalFunction* gtf);

#endif /* RBTALG_H_ */
