/*
 * vec.h
 *
 *  Created on: Feb 6, 2010
 *      Author: dnlennon
 */

#ifndef VEC_H_
#define VEC_H_

#include "queue.h"
#include "stack.h"

#include <assert.h>
#include <iostream>
using namespace std;

template <class T>
class vec
{
public:
	vec()
	{
		myLen = 0;
		myMaxIndex = 0;
		myv = NULL;
	}

	vec(const int sz)
	{
		myLen = sz;
		myv = new T[length()];
		myMaxIndex = length();
	}

	vec(const vec<T>& v0)
	{
		if(this == &v0)
			return;

		myLen = v0.length();
		myMaxIndex = length();
		myv = new T[length()];
		for(int i=0; i<length(); i++)
			myv[i] = v0[i];
	}

	void assign(const vec<T>& v0)
	{
		if(this == &v0)
			return;

		resize(v0.length());

		for(int i=0; i<length(); i++)
			myv[i] = v0[i];
	}

	void assign(const queue<T>& q)
	{
		resize(q.size());

		queue_node<T>* tmp = q.head;
		int i=0;
		while(tmp != NULL)
		{
			assert(i < length());

			// Object's operator= will be called.
			myv[i++] = tmp->data();
			tmp = tmp->next();
		}
	}

	void assign(const stack<T>& s)
	{
		resize(s.size());

		stack_node<T>* tmp = s.head;
		int i=0;
		while(tmp != NULL)
		{
			assert(i < length());

			// Object's operator= will be called.
			myv[i++] = tmp->data();
			tmp = tmp->next();
		}
	}

	void resize(int i)
	{
		if(myv != NULL)
			delete [] myv;

		myLen = i;
		myMaxIndex = length();
		myv = new T[length()];
	}

	~vec()
	{
		if(myv != NULL)
			delete [] myv;
		myv = NULL;
	}

	operator T*() { return myv; }

	T& operator[] (int i)
	{
		assert(i>=0 && i < length());
		return myv[i];
	}

	const T& operator[] (int i) const
	{
		assert(i>=0 && i < length());
		return myv[i];
	}

	int length(void) const
	{
		return myLen;
	}

	int size(void) const
	{
		return length();
	}

	void sort(vec<int>& order,  /* output */
			  int (*cmpfn)(T&, T&)
			 )
	{
		order.resize(length());
		for(int i=0; i<length(); i++)
			order[i] = i;

		for(int i=0; i<length(); i++)
		{
			for(int j=i+1; j<length(); j++)
			{
				if( cmpfn(myv[order[i]], myv[order[j]]) > 0 )
				{
					int k = order[i];
					order[i] = order[j];
					order[j] = k;
				}
			}
		}
	}

	void insertAtEnd(T& o)
	{
		if(length() + 1 >= myMaxIndex)
		{
			vec<T> tmp;
			tmp.assign(*this);
			resize(myMaxIndex * 2);
			assign(tmp);
		}
		else
		{
			operator[](myLen) = o;
			myLen++;
		}
	}

	template <class S> friend ostream& operator<<(ostream& o, const vec<S>& v);

private:
	vec<T> operator=(const vec<T>& v0) {}
	T* myv;
	int myLen;
	int myMaxIndex;
};

template <class T>
ostream& operator<<(ostream& o, const vec<T>& v)
{
	for(int i=0; i<v.length(); i++)
	{
		o<<v[i];
		if(i < v.length() - 1)
			o<<" ";
	}
	return o;
}


#endif /* VEC_H_ */
