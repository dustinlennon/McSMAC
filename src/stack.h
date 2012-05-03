/*
 * stack.h
 *
 *  Created on: Apr 13, 2010
 *      Author: dnlennon
 */

#ifndef STACK_H_
#define STACK_H_

#include "queuestackBase.h"

#include <stdlib.h>
#include <assert.h>

#include <iostream>
using namespace std;

template <class T> class stack;
template <class T> class vec;

/*  The stack class maintains a copy of each object pushed on to the stack.
	Calling the object's copy constructor may be expensive.  */

template <class T>
class stack_node
{
public:
	// This calls the object's copy constructor.
	stack_node(T& t, stack_node<T>* n) : mydata(t), mynext(n) { }

	stack_node* next(void) const
	{
		return mynext;
	}

	const T& data(void) const
	{
		// This will call T's copy constructor.
		return mydata;
	}

private:
	friend class stack<T>;
	T mydata;
	stack_node<T>* mynext;
};

template <class T>
class stack : public queuestackBase<T>
{
public:
	stack(void)
	{
		head = NULL;
		mysize = 0;
	}

	~stack()
	{
		while(!empty())
			pop();
	}

	void push(T& t)
	{
		// This calls the object's copy constructor.
		stack_node<T>* tmp = new stack_node<T>(t, head);

		mysize++;
		head = tmp;
	}

	void pop(void)
	{
		assert(head != NULL);

		stack_node<T>* tmp = head;
		head = head->mynext;

		mysize--;
		delete tmp;
	}

	T& top(void)
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

	virtual T& next(void)
	{
		return top();
	}

	virtual void insert(T& t)
	{
		push(t);
	}

	virtual void remove(void)
	{
		pop();
	}


	template <class S> friend ostream& operator<<(ostream& o, const stack<S>& q);
	friend class vec<T>;
private:
	stack<T>& operator=(const stack<T>&) { }
	stack(const stack<T>&) { }

	stack_node<T> *head;
	int mysize;
};

template <class T>
ostream& operator<<(ostream& o, const stack<T>& s)
{
	stack_node<T>* tmp = s.head;
	o<<"("<<s.size()<<") ";
	while(tmp != NULL)
	{
		o<<tmp->data();
		tmp = tmp->next();
		if(tmp != NULL)
			o<<" ";
	}
	return o;
}

#endif /* STACK_H_ */
