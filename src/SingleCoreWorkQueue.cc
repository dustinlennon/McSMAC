/*
 * SingleCoreWorkQueue.cc
 *
 *  Created on: Oct 16, 2010
 *      Author: dnlennon
 */

#include "SingleCoreWorkQueue.h"

void SingleCoreWorkQueue::results(void)
{
	while(_head != NULL)
	{
		WQfunc fn = _head->fn;
		fn(_head->res, _head->data);

		wqnode* tmp = _head;
		_head = _head->next;
		delete tmp;

		if(_head == NULL)
		{
			_tail = NULL;
		}
	}
}
