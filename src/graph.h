/*
 * Graph.h
 *
 *  Created on: May 4, 2010
 *      Author: dnlennon
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <string>
using namespace std;

#include "search_set.h"
#include "vec.h"
#include "dictionary.h"
#include "dataobject.h"


typedef class edge edge;
typedef class node node;

class node : public DOPuser
{
public:
	node(string id, DOP dop) : DOPuser(dop), myid(id) { }

	int numEdges(void) const
	{
		return deb.size();
	}

	edge* getep(string s);

	void setep(string& s, edge* ep)
	{
		deb.insert(s, ep);
	}

	string getId(void) const
	{
		return myid;
	}

	void addEdge(edge* ep);
	void removeEdge(edge* ep);

	void getAdjNodeStrings(vec<string>& vds /* output */);

private:
	node(const node&) {}
	node& operator=(const node&) { return *this; }

	dictionary<string, edge*> deb;
	string myid;
};


class edge: public DOPuser
{
public:
	//edge(node* n1, node* n2, double len, DOP dop);
	edge(node* n1, node* n2, DOP dop);
	~edge();

	node* getnp(int i);
	node* getnp(const string& s);

	void setnp(int i, node* np);
	void setnp(const string& s, node* np);

	node* otherNode(node* n);

	string getId(void) const
	{
		return myid;
	}

	static string edgeString(const string& s1, const string& s2);

private:
	edge(const edge&) {}
	edge& operator=(const edge&) { return *this; }

	vec<node*> vnp;
	string myid;
};


class graph: public DOPuser
{
	friend ostream& operator<<(ostream&o, const graph&);

public:
	graph() {}
	virtual ~graph();

	/* Graph copying means the following:
	 * 1.  (sub)connectivity strucutres are deep copied
	 * 2.  associated DOP's are shallow copied (with respective ref. counter increments)
	 */
	graph(const graph& g);
	graph& operator=(const graph& g);

	node* addNode(string id)
	{
		return addNode(id, NULL);
	}

	edge* addEdge(const string& s1, const string& s2)
	{
		return addEdge(s1, s2, NULL);
	}

	node* addNode(const string& id, DOP o);
	edge* addEdge(const string& s1, const string& s2, DOP dop);

	/* allows operator=() to addEdges incrementally and check sanity conditions posthoc */
	edge* addEdgeSaneCheckOpt(const string& s1, const string& s2, DOP dop, bool doSaneCheck);

	void removeNode(const string& s1);
	void removeEdge(const string& s1, const string& s2);

	bool findEdge(const string& s1, const string& s2, edge*& epout /* output */);
	bool findNode(const string& s1, node*& nout /* output */);

	/* Returns a vector containing graph objects.  Each of these graphs is connected. */
	void connectedComponents(vec< graph >& vcc /* output */ );

	/* Useful for derived classes:  initialize with graph constructor, sanity check
	 * with derived version of sanity.  incremental == false specifies that the
	 * sanity checks of all classes in the hierarchy are performed.  incremental == true
	 * specifies that only the local sanity check in the current derived class is performed. */
	virtual bool checkSanity(bool incremental) { return true; }

	virtual void deepClean(void);

	void nodeList(vec<node*>& vnp);
	void edgeList(vec<edge*>& vep);

	int numNodes(void) const
	{
		return np_set.size();
	}

	int numEdges(void) const
	{
		return ep_set.size();
	}

	void allDopMakeOwner(void);

protected:

	/* This may throw RecoverableAssertException.  corresponding catch code should
	 * clean up the object.  no guarantees are made on internal state after the throw.
	 */
	void sanityCheckWithException(bool incremental);

	dictionary<string, node*> node_dict;

	search_set<edge*> ep_set;
	search_set<node*> np_set;
};


#endif /* GRAPH_H_ */
