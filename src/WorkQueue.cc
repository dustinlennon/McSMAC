/*
 * WorkQueue.cc
 *
 *  Created on: Oct 16, 2010
 *      Author: dnlennon
 */

#include "WorkQueue.h"

void WorkQueue::enqueue(void* res, void* data)
{
	if(_head == NULL)
	{
		_tail = new wqnode;
		_tail->fn = _fn;
		_tail->res = res;
		_tail->data = data;
		_tail->next = NULL;
		_head = _tail;
	}
	else
	{
		_tail->next = new wqnode;
		_tail = _tail->next;
		_tail->fn = _fn;
		_tail->res = res;
		_tail->data = data;
		_tail->next = NULL;
	}
	return;
}
