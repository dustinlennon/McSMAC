/*
 * SingleCoreWorkQueue.h
 *
 *  Created on: Oct 16, 2010
 *      Author: dnlennon
 */

#ifndef SINGLECOREWORKQUEUE_H_
#define SINGLECOREWORKQUEUE_H_

#include "WorkQueue.h"

class SingleCoreWorkQueue : public WorkQueue
{
public:
	SingleCoreWorkQueue(WQfunc fn) : WorkQueue(fn) {}
	virtual void results(void);
};


#endif /* SINGLECOREWORKQUEUE_H_ */
