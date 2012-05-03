/*
 * dictionary.h
 *
 *  Created on: May 4, 2010
 *      Author: dnlennon
 */

#ifndef DICTIONARY_H_
#define DICTIONARY_H_

#include "vec.h"

#include <map>
using namespace std;

template <class K, class V>
class dictionary
{
public:
	dictionary() {}

	void insert(K k, V v)
	{
		my_map[k] = v;
	}

	void remove(K k)
	{
		my_map.erase(k);
	}

	bool lookup(K k, V& vout  /* output */)
	{
		typename std::map<K,V>::iterator iter = my_map.find(k);
		if(iter == my_map.end())
			return false;

		vout = iter->second;
		return true;
	}

	int size(void) const
	{
		return my_map.size();
	}

	void getKeysVec( vec<K>& kv )
	{
		kv.resize( size() );
		int i=0;
		typename std::map<K,V>::iterator iter = my_map.begin();
		while(iter != my_map.end())
		{
			assert(i < kv.length() );
			kv[i++] = iter->first;
			iter++;
		}
	}

	void getValsVec( vec<V>& vv )
	{
		vv.resize( size() );
		int i=0;
		typename std::map<K,V>::iterator iter = my_map.begin();
		while(iter != my_map.end())
		{
			assert(i < vv.length() );
			vv[i++] = iter->second;
			iter++;
		}
	}

private:
	map<K,V> my_map;
};

#endif /* DICTIONARY_H_ */
