/*
 * MultiCoreWorkQueue.h
 *
 *  Created on: Oct 16, 2010
 *      Author: dnlennon
 */

#include "MultiCoreWorkQueue.h"
extern  "C"
{
#include "workq.h"
}
#include "errors.h"

MultiCoreWorkQueue::MultiCoreWorkQueue(WQfunc fn) : WorkQueue(fn)
{
	int status = workq_init (&_workq, MAX_THREADS, _fn);
	if (status != 0)
		err_abort (status, "Init work queue");
}

void MultiCoreWorkQueue::results(void)
{
	while(_head != NULL)
	{
		workq_add(&_workq, _head->res, _head->data);

		wqnode* tmp = _head;
		_head = _head->next;
		delete tmp;

		if(_head == NULL)
		{
			_tail = NULL;
		}
	}

	int status = workq_destroy (&_workq);
	if (status != 0)
		err_abort (status, "Destroy work queue");

}
