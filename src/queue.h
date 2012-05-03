/*
 * queue.h
 *
 *  Created on: Feb 5, 2010
 *      Author: dnlennon
 */

#ifndef QUEUE_H_
#define QUEUE_H_

#include "queuestackBase.h"

#include <stdlib.h>
#include <assert.h>

#include <fstream>
#include <iostream>
using namespace std;

template <class T> class queue;
template <class T> class vec;

/*  The queue class maintains a copy of each object inserted into the queue.
	Calling the object's copy constructor may be expensive.  */

template <class T>
class queue_node
{
public:
	// This will call T's copy constructor.
	queue_node(T& t, queue_node<T>* n) : mydata(t), mynext(n) { }

	queue_node* next(void) const
	{
		return mynext;
	}

	const T& data(void) const
	{
		// This will call T's copy constructor.
		return mydata;
	}

private:
	friend class queue<T>;
	T mydata;
	queue_node<T>* mynext;
};

template <class T>
class queue : public queuestackBase<T>
{
public:
	queue(void)
	{
		head = NULL;
		tail = NULL;
		mysize = 0;
	}

	~queue()
	{
		while(!empty())
			remove();
	}

	virtual void insert(T& t)
	{
		if(head == NULL && tail == NULL)
		{
			// This will call T's copy constructor.
			tail = new queue_node<T>(t, NULL);
			head = tail;
		}
		else
		{
			// This will call T's copy constructor.
			tail->mynext = new queue_node<T>(t, NULL);
			tail = tail->mynext;
		}
		mysize++;
	}

	virtual void remove(void)
	{
		queue_node<T>* tmp;
		if(head == NULL)
		{
			return;
		}
		else if(head == tail)
		{
			delete head;
			head = NULL;
			tail = NULL;
		}
		else
		{
			tmp = head;
			head = head->mynext;
			delete tmp;
		}
		mysize--;
		return;
	}

	virtual T& next(void)
	{
		return front();
	}

	T& front(void)
	{
		assert(head != NULL);
		return head->mydata;
	}

	bool empty(void)
	{
		return head == NULL;
	}

	int size(void) const
	{
		return mysize;
	}

	template <class S> friend ostream& operator<<(ostream& o, const queue<S>& q);
	friend class vec<T>;

private:
	queue_node<T> *head, *tail;
	int mysize;
};

template <class T>
ostream& operator<<(ostream& o, const queue<T>& q)
{
	queue_node<T>* tmp = q.head;
	o<<"("<<q.size()<<") ";
	while(tmp != NULL)
	{
		o<<tmp->data();
		tmp = tmp->next();
		if(tmp != NULL)
			o<<" ";
	}
	return o;
}


#endif /* QUEUE_H_ */
