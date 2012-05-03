/*
 * graph.cc
 *
 *  Created on: May 4, 2010
 *      Author: dnlennon
 */

#include "graph.h"
#include "search_set.h"
#include "dictionary.h"
#include "dataobject.h"
#include "RecoverableAssertException.h"

////////////////////////////////////////////////////////////////////////////////

void node::addEdge(edge* ep)
{
	deb.insert(ep->getId(), ep);
}

void node::removeEdge(edge* ep)
{
	deb.remove(ep->getId());
}

void node::getAdjNodeStrings(vec<string>& vds /* output */)
{
	vec<edge*> tmp;
	deb.getValsVec( tmp );

	vds.resize( tmp.length() );
	for(int i=0; i<vds.length(); i++)
	{
		vds[i] = tmp[i]->otherNode(this)->getId();
	}
}

edge* node::getep(string s)
{
	edge* ep = NULL;

	/*
	vec<string> vs;
	deb.getKeysVec(vs);
	cout<<getId()<<":\t";
	for(int i=0; i<vs.length(); i++)
		cout<<vs[i]<<",";
	cout<<endl;
	*/

	if( deb.lookup(s, ep) != false )
		return ep;

	if( deb.lookup( edge::edgeString(getId(), s) , ep) != false)
		return ep;

	if( deb.lookup( edge::edgeString(s, getId()) , ep) != false)
		return ep;

	return NULL;
}


////////////////////////////////////////////////////////////////////////////////

edge::edge(node* n1, node* n2, DOP dop) : DOPuser(dop)
{
	vnp.resize(2);
	vnp[0] = n1;
	vnp[1] = n2;

	mydop = dop;
	myid = edgeString( n1->getId(), n2->getId() );	

	n1->addEdge(this);
	n2->addEdge(this);
}

edge::~edge()
{
	vnp[0]->removeEdge(this);
	vnp[1]->removeEdge(this);
}

string edge::edgeString(const string& s1, const string& s2)
{
	return s1 + " -- " + s2;
}

node* edge::getnp(int i)
{
	assert(i == 0 || i == 1);
	return vnp[i];
}

node* edge::getnp(const string& s)
{
	if(vnp[0]->getId() == s)
		return vnp[0];
	else if(vnp[1]->getId() == s)
		return vnp[1];
	else
		assert(0);

	return NULL;
}

void edge::setnp(int i, node* np)
{
	vnp[i] = np;
}

void edge::setnp(const string& s, node* np)
{
	if(vnp[0]->getId() == s)
		vnp[0] = np;
	else if(vnp[1]->getId() == s)
		vnp[1] = np;
	else
		assert(0);

	return;
}

node* edge::otherNode(node* n)
{
	if(vnp[0] == n)
		return vnp[1];
	else if(vnp[1] == n)
		return vnp[0];
	else
		return NULL;
}

////////////////////////////////////////////////////////////////////////////////

node* graph::addNode(const string& id, DOP dop)
{
	node* n;
	if( node_dict.lookup(id, n) == false )
	{
		n = new node(id, dop);
		np_set.insert(n);
		node_dict.insert(id, n);
	}
	else
	{
		n->setDOP(dop);
	}

	return n;
}

edge* graph::addEdge(const string& s1, const string& s2, DOP dop)
{
	return addEdgeSaneCheckOpt(s1, s2, dop, true);
}


edge* graph::addEdgeSaneCheckOpt(const string& s1, const string& s2, DOP dop, bool doSaneCheck)
{
	assert(s1 != s2);

	edge* ep;
	if ( findEdge(s1, s2, ep) == false )
	{
		node *n1, *n2;
		if( node_dict.lookup(s1, n1) == false )
			n1 = addNode(s1);

		if( node_dict.lookup(s2, n2) == false )
			n2 = addNode(s2);

		ep = new edge(n1, n2, dop);
		ep_set.insert(ep);
	}
	else
	{
		ep->setDOP(dop);
	}

	if(doSaneCheck == true)
	{
		sanityCheckWithException(false);
	}

	return ep;
}

bool graph::findEdge(const string& s1, const string& s2, edge*& epout /* output */)
{
	epout = NULL;
	node *n1;
	if( node_dict.lookup(s1, n1) == false )
		return false;

	node *n2;
	if( node_dict.lookup(s2, n2) == false )
		return false;

	string schk1 = edge::edgeString(s1 , s2);
	string schk2 = edge::edgeString(s2,  s1);

	epout = n1->getep(schk1);
	if(epout == NULL)
		epout = n1->getep(schk2);

	if(epout == NULL)
		return false;

	return true;
}

bool graph::findNode(const string& s1, node*& nout /* output */)
{
	nout = NULL;
	node* n;
	if( node_dict.lookup(s1, n) == false)
		return false;

	nout = n;
	return true;
}

void graph::removeNode(const string& n1)
{
	node* np;
	if( findNode(n1, np) == false)
		return;

	vec<string> vds;
	np->getAdjNodeStrings(vds);
	for(int i=0; i<vds.length(); i++)
	{
		removeEdge(np->getId(), vds[i]);
	}

	node_dict.remove( np->getId() );
	np_set.remove(np);

	delete np;
}

void graph::removeEdge(const string& s1, const string& s2)
{
	edge* ep;
	if ( findEdge(s1, s2, ep) == false )
		return;

	ep_set.remove(ep);
	delete ep;
}

void graph::deepClean(void)
{
	vec<edge*> vep;
	ep_set.getVec(vep);
	for(int i=0; i<vep.length(); i++)
	{
		edge* ep = vep[i];
		removeEdge( ep->getnp(0)->getId(), ep->getnp(1)->getId() );
	}

	vec<node*> vnp;
	np_set.getVec(vnp);
	for(int i=0; i<vnp.length(); i++)
	{
		removeNode( vnp[i]->getId() );
	}
}

graph::~graph()
{
	deepClean();
	/*
	vec<edge*> vep;
	ep_set.getVec(vep);
	for(int i=0; i<vep.length(); i++)
	{
		delete vep[i];
	}

	vec<node*> vnp;
	np_set.getVec(vnp);
	for(int i=0; i<vnp.length(); i++)
	{
		delete vnp[i];
	}
	*/
}

void graph::connectedComponents(vec< graph >& vcc /* output */ )
{
	queue< graph > gq;
	search_set<node*> np_global;

	vec<node*> vnp;
	np_set.getVec(vnp);

	for(int i=0; i<vnp.length(); i++)
	{
		if( np_global.inSet( vnp[i] ) )
			continue;

		/* BFS to find the nodes connected to vnp[i] */
		graph cc;
		cc.addNode(vnp[i]->getId(), vnp[i]->getDOP());
		
		search_set<node*> np_local;
		np_local.insert(vnp[i]);
		np_global.insert(vnp[i]);

		queue<node*> npq;
		npq.insert(vnp[i]);
		while( npq.empty() == false )
		{
			node* np_src = npq.front();

			vec<string> vns;
			np_src->getAdjNodeStrings(vns);
			for(int j=0; j<vns.length(); j++)
			{
				edge* ep;
				if( findEdge(np_src->getId(), vns[j], ep) == false )
					assert(0);

				node* np_dest = ep->otherNode(np_src);

				if( np_local.inSet(np_dest) == true )
					continue;

				cc.addNode(np_dest->getId(), np_dest->getDOP());

				np_local.insert(np_dest);
				np_global.insert(np_dest);
				npq.insert(np_dest);
			}
			npq.remove();
		}

		/* cc has the correct nodes connected to vnp[i].  now add the edges */
		vec<node*> cc_vnp;
		cc.np_set.getVec(cc_vnp);
		for(int j=0; j<cc_vnp.length(); j++)
		{
			node* np;
			if( findNode(cc_vnp[j]->getId(), np) == false )
				assert(0);

			vec<string> cc_vs;
			np->getAdjNodeStrings(cc_vs);
			for(int k=0; k<cc_vs.length(); k++)
			{
				edge* ep;
				if( findEdge(cc_vnp[j]->getId(), cc_vs[k], ep) == false )
					assert(0);

				cc.addEdge(cc_vnp[j]->getId(), cc_vs[k], ep->getDOP() );
			}
		}

		gq.insert(cc);
	}
	vcc.assign(gq);
}

graph::graph(const graph& g)
{
	vec<node*> vnp;
	g.np_set.getVec(vnp);
	for(int i=0; i<vnp.length(); i++)
	{
		node* np = vnp[i];
		addNode(np->getId(), np->getDOP());
	}

	vec<edge*> vep;
	g.ep_set.getVec(vep);
	for(int i=0; i<vep.length(); i++)
	{
		edge* ep = vep[i];
		addEdgeSaneCheckOpt(ep->getnp(0)->getId(),
					ep->getnp(1)->getId(),
					ep->getDOP(),
					false);
	}
}

void graph::sanityCheckWithException(bool incremental)
{
	checkSanity(incremental);
}

graph& graph::operator=(const graph& g)
{
	if(this == &g)
		return *this;

	deepClean();

	vec<node*> vnp;
	g.np_set.getVec(vnp);
	for(int i=0; i<vnp.length(); i++)
	{
		node* np = vnp[i];
		addNode(np->getId(), np->getDOP());
	}

	vec<edge*> vep;
	g.ep_set.getVec(vep);
	for(int i=0; i<vep.length(); i++)
	{
		edge* ep = vep[i];
		string ns1 = ep->getnp(0)->getId();
		string ns2 = ep->getnp(1)->getId();
		addEdgeSaneCheckOpt( ns1, ns2, ep->getDOP(), false );
	}

	setDOP( g.getDOP() );
	//sanityCheckWithException(false);

	return *this;
}


void graph::nodeList(vec<node*>& vnp)
{
	np_set.getVec(vnp);
}

void graph::edgeList(vec<edge*>& vep)
{
	ep_set.getVec(vep);
}

////////////////////////////////////////////////////////////////////////////////

ostream& operator<<(ostream& o, const graph& g)
{
	o<<"\n----------------------"<<endl;
	o<<"\ngraph"<<endl;
	o<<hex<<"0x"<<(unsigned long)(g.getDOP().ptr())<<dec<<endl;

	o<<"\nedges"<<endl;
	vec<edge*> vep;
	g.ep_set.getVec(vep);
	for(int i=0; i<vep.length(); i++)
	{
		o<<vep[i]->getnp(0)->getId()<<" to "<<vep[i]->getnp(1)->getId()<<" (aka "<<vep[i]->getId()<<" )\t"
		<<hex<<"0x"<<(unsigned long)(vep[i]->getDOP().ptr())<<dec<<endl;
	}

	o<<"\nnodes"<<endl;
	vec<node*> vnp;
	g.np_set.getVec(vnp);
	for(int i=0; i<vnp.length(); i++)
	{
		o<<vnp[i]->getId()<<"\t"<<hex<<"0x"<<(unsigned long)(vnp[i]->getDOP().ptr())<<dec<<endl;
	}
	return o;
}

void graph::allDopMakeOwner(void)
{
	DOP dop = getDOP();
	dop.makeOwner();
	setDOP(dop);

	vec<node*> vnp;
	nodeList(vnp);
	for(int i=0; i<vnp.length(); i++)
	{
		dop = vnp[i]->getDOP();
		dop.makeOwner();
		vnp[i]->setDOP(dop);
	}

	vec<edge*> vep;
	edgeList(vep);
	for(int i=0; i<vep.length(); i++)
	{
		dop = vep[i]->getDOP();
		dop.makeOwner();
		vep[i]->setDOP( dop );
	}
}
