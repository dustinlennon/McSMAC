/*
 * WorkQueue.h
 *
 *  Created on: Oct 16, 2010
 *      Author: dnlennon
 *
 *      WorkQueues
 */

#ifndef WORKQUEUE_H_
#define WORKQUEUE_H_

// const int MAX_THREADS=8;

#include <stdlib.h>

typedef void (*WQfunc)(void* result, void* data);

class WorkQueue
{
protected:
	typedef struct wqnode
	{
		wqnode* next;
		WQfunc fn;
		void* res;
		void* data;
	} wqnode;


public:
	WorkQueue(WQfunc fn)
	{
		_head = _tail = NULL;
		_fn = fn;
	}

	virtual void enqueue(void* res, void* data);
	virtual void results(void) = 0;
	virtual ~WorkQueue() {}

protected:
	WQfunc _fn;

	wqnode* _head;
	wqnode* _tail;
};


#endif /* WORKQUEUE_H_ */
