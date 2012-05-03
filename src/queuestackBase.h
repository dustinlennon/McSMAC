/*
 * queuestackbase.h
 *
 *  Created on: May 16, 2010
 *      Author: dnlennon
 */

#ifndef QUEUESTACKBASE_H_
#define QUEUESTACKBASE_H_

template <class T>
class queuestackBase
{
public:
	virtual void insert(T& t) = 0;
	virtual void remove(void) = 0;
	virtual T& next(void) = 0;
	virtual bool empty(void) = 0;
};

#endif /* QUEUESTACKBASE_H_ */
