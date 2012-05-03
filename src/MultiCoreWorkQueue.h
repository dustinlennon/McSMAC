/*
 * MultiCoreWorkQueue.h
 *
 *  Created on: Oct 16, 2010
 *      Author: dnlennon
 */

#ifndef MULTICOREWORKQUEUE_H_
#define MULTICOREWORKQUEUE_H_

#include "WorkQueue.h"
extern "C"
{
#include "workq.h"
}

class MultiCoreWorkQueue : public WorkQueue
{
public:
	MultiCoreWorkQueue(WQfunc fn);
	virtual void results(void);

	virtual ~MultiCoreWorkQueue() {}

protected:
	workq_t _workq;
};


#endif /* MULTICOREWORKQUEUE_H_ */
