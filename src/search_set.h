/*
 * rb_tree.h
 *
 *  Created on: May 4, 2010
 *      Author: dnlennon
 */

#ifndef SEARCH_SET_H_
#define SEARCH_SET_H_

#include "vec.h"
#include <map>

template <class T>
class search_set
{
public:
	search_set() { }

	void insert(T t)
	{
		my_map[t] = true;
	}

	void remove(T t)
	{
		my_map.erase(t);
	}

	bool inSet(T t)
	{
		return my_map.find(t) != my_map.end();
	}

	void getVec(vec<T>& v  /* output */)  const
	{
		v.resize(my_map.size());
		typename map<T, bool>::const_iterator it = my_map.begin();
		for(int i=0; it != my_map.end(); it++, i++)
		{
			v[i] = it->first;
		}
	}

	int size(void) const
	{
		return my_map.size();
	}

private:
	search_set& operator=(const search_set& st) { }
	search_set(search_set& st) { }

	map<T, bool> my_map;
};


#endif /* SEARCH_SET_H_ */
