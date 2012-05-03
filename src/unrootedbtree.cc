/*
 * unrootedbtree.cc
 *
 *  Created on: May 14, 2010
 *      Author: dnlennon
 */

#include "unrootedbtree.h"
#include "RecoverableAssertException.h"

bool unrootedbtree::checkSanity(bool incremental)
{
	if(incremental == false)
	{
		generaltree::checkSanity(false);
	}

	vec<node*> vnp;
	np_set.getVec(vnp);
	for(int i=0; i<vnp.length(); i++)
	{
		int ne = vnp[i]->numEdges();
		if( ne == 3 || ne == 1)
			continue;
		else
			throw RecoverableAssertException(__FUNCTION__, __FILE__, __LINE__);
	}

	return true;
}
