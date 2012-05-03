/*
 * test.cc
 *
 *  Created on: May 4, 2010
 *      Author: dnlennon
 */

#include <string>
using namespace std;

#include "queue.h"
#include "stack.h"
#include "vec.h"
#include "dummy.h"
#include "dictionary.h"
#include "search_set.h"
#include "dataobject.h"
#include "graph.h"
#include "generaltree.h"
#include "unrootedbtree.h"
#include "rootedbtree.h"
#include "graphalg.h"
#include "RecoverableAssertException.h"
#include "gslvec.h"
#include "gslmat.h"
#include "test.h"

#include "fastmap.h"
#include "phylo_util.h"
#include "phylo_objs.h"

#include <math.h>

#include <gsl/gsl_matrix.h>

#include <fstream>
#include <iostream>
using namespace std;

void reopen(ifstream& f, string s)
{
	f.close();
	f.clear();
	f.open(s.c_str());
}

void recycle(ofstream& f, string s)
{
	f.close();
	f.clear();
	f.open(s.c_str());
}


void incfn(gslvec v)
{
	for(int i=0; i<v.length(); i++)
	{
		double d = gsl_vector_get(v, i);
		d += 1.0;
		gsl_vector_set(v, i, d);
	}
}


int gsl_test(int argc, char** argv)
{
	ifstream fin("data/sim.data.mat");
	gslmat sdm(fin);

	gsl_vector_view gvv = gsl_matrix_row(sdm, 0);
	gslvec v(gvv);

//	gsl_vector_const_view gvcv = gsl_matrix_const_row(sdm, 1);

	cout<<"gvv: "<<gvv<<endl;
	incfn(gvv);
	cout<<"gvv: "<<gvv<<endl;
	incfn(v);
	cout<<"gvv: "<<gvv<<endl;
	cout<<"v: "<<v<<endl;

	return 0;
}

int run_container_tests(int argc, char** argv)
{
	/*
	simple_container_test(argc, argv);
	dop_test(argc, argv);
	graph_test(argc, argv);
	generaltree_test(argc, argv);
	ubt_test(argc, argv);
	rbt_test(argc, argv);
	bfs_test(argc, argv);
	gsl_test(argc, argv);
	*/
	rooted_phylo(argc, argv);

	return 0;
}

int simple_container_test(int argc, char** argv)
{
	cout<<"simple_container_test:"<<endl;

	dummy d1(1);
	dummy d2(2);
	dummy d3(3);

	queue<dummy> q;
	stack<dummy> s;
	q.insert(d1);
	q.insert(d2);
	q.insert(d3);
	cout<<q.front()<<endl;

	cout<<q<<endl;

	s.push(d1);
	s.push(d2);
	s.push(d3);
	cout<<s.top()<<endl;

	cout<<s<<endl;

	vec<dummy> vdq;
	vdq.assign(q);

	vec<dummy> vds;
	vds.assign(s);

	cout<<vdq<<endl;
	cout<<vds<<endl;

	cout<<vdq[0]<<endl;
	cout<<vds[0]<<endl;

	dictionary<char, int> d;
	d.insert('a', 1);
	d.insert('b', 2);
	d.insert('c', 3);

	cout<<"dictionary d has size "<<d.size()<<endl;

	int i;
	d.remove('b');
	cout<<"dictionary d has size "<<d.size()<<endl;
	for(char ch='a'; ch<='d'; ch++)
	{
		if( d.lookup(ch, i) == false )
		{
			cout<<ch<<" not found in dictionary"<<endl;
		}
		else
		{
			cout<<ch<<" found in dictionary with value "<<i<<endl;
		}
	}

	search_set<char> ss;
	ss.insert('a');
	ss.insert('b');
	ss.insert('c');

	ss.remove('b');
	ss.insert('d');
	for(char ch='a'; ch<='d'; ch++)
	{
		if( ss.inSet(ch) == true )
		{
			cout<<ch<<" is in the set"<<endl;
		}
		else
		{
			cout<<ch<<" in not in the set"<<endl;
		}
	}

	vec<char> vc;
	ss.getVec(vc);
	for(i=0; i<vc.length(); i++)
	{
		cout<<i<<" "<<vc[i]<<endl;
	}


	return 0;
}


int dop_test(int argc, char** argv)
{
	cout<<"dop_test:"<<endl;

	dataObject* dop = new dataObject("first dop");
	DOP d1(dop);
	DOP d4 = d1;
	DOP d5;
	d5 = d1;

	if(true)
	{
		DOP d2 = d1;
		DOP d3;
		d3 = d2;
		d3.makeOwner();
		d2 = d3;
		d2 = NULL;
		cout<<*d3<<endl;
	}

	return 0;
}

class specialNodeObject : public dataObject
{
	friend ostream& operator<<(ostream& o, const specialNodeObject& snop);
public:
	specialNodeObject(string id, int x) : dataObject(id) { myx = x; }

private:
	dataObject* deepCopy()
	{
		specialNodeObject* snop= new specialNodeObject(myid, myx);
		return snop;
	}

	int myx;
};

ostream& operator<<(ostream& o, const specialNodeObject& snop)
{
	o<<snop.myid<<" "<<snop.myx;
	return o;
}

int graph_test(int argc, char** argv)
{
	cout<<"graph_test:"<<endl;

	graph g;

	string s1("n1");
	g.addNode("n1");
	g.addNode("n2", new specialNodeObject("two", 2));
	g.addEdge("n1", "n2", new dataObject("edge"));

	cout<<g<<endl;

	node* n;
	specialNodeObject* snop;
	g.findNode("n2", n);
	if(n != NULL)
	{
		DOP d = n->getDOP();

		snop = dynamic_cast<specialNodeObject*>( d.ptr() );
		cout<<*snop<<endl;

		g.removeNode("n2");
		cout<<g<<endl;

		// when we leave the scope of the if {}, snop will point to nowhere
		// since we removed "n2" from graph!!
	}
	// THIS WILL DEREFERENCE FREE'D MEMORY:  cout<<*snop<<endl;
	snop = NULL;

	g.addNode("n1", new specialNodeObject("one", 1));
	cout<<g<<endl;

	g.addNode("n1", new specialNodeObject("zero", 0));
	cout<<g<<endl;

	edge* e;
	g.findEdge("n1", "n2", e);
	g.removeEdge("n1", "n2");
	cout<<g<<endl;

	g.addEdge("n1", "n2", new dataObject("edge2"));
	cout<<g<<endl;

	g.removeEdge("n2", "n1");
	cout<<g<<endl;

	g.addEdge("n1", "n2");
	g.addEdge("n2", "n3");
	g.addEdge("n3", "n4");
	g.addEdge("n4", "n1");
	g.addEdge("n2", "n4");
	g.addEdge("n1", "n3");

	g.addEdge("n6", "n7", new dataObject("uhoh"));
	g.addEdge("n6", "n8");
	g.addEdge("n7", "n8");

	g.addEdge("n9", "n10");

	g.addNode("nxx", new dataObject("fp9"));
	g.addNode("mxx", new dataObject("bw2"));
	g.addEdge("nxx", "mxx");

	cout<<g;

	vec<graph> cc;
	g.connectedComponents(cc);	

	for(int i=0; i<cc.length(); i++)
	{
		cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		cout<<"component "<<i<<endl;
		cout<<cc[i];
	}

	return 0;
}

int generaltree_test(int argc, char** argv)
{
	cout<<"generaltree_test:"<<endl;
	graph g;
	g.addEdge("1", "2");
	g.addEdge("2", "3");
	g.addEdge("2", "4");
	g.addEdge("2", "5");
	g.addEdge("2", "5");
	g.addEdge("6", "5");
	g.addEdge("7", "5");
	g.addEdge("7", "8");

	graph g2(g);
	g2.addEdge("8", "6");

	cout<<g<<endl;
	cout<<g2<<endl;

	generaltree gt1(g);
	generaltree gt2;
	try
	{
		gt2 = g2;
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RecoverableAssertException:  "<<e.what()<<endl;
	}
	gt2 = gt1;
	generaltree gt3(gt2);

	try
	{
		gt3.addEdge("6", "8");
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RecoverableAssertException:  "<<e.what()<<endl;
	}

	return 0;
}

int ubt_test(int argc, char** argv)
{
	cout<<"ubt_test:"<<endl;
	graph g0;
	g0.addEdge("1", "2");
	g0.addEdge("2", "3");
	g0.addEdge("2", "5");
	g0.addEdge("2", "5");
	g0.addEdge("6", "5");
	g0.addEdge("7", "5");
	g0.addEdge("7", "8");

	unrootedbtree ubt0;
	try
	{
		ubt0 = g0;
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RAE: "<<e.what()<<"\tubt1 initialization failed"<<endl;
		ubt0.deepClean();
	}

	cout<<"what is the state of ubt0 after throwing exception?"<<endl;
	cout<<ubt0<<endl;

	g0.addEdge("7", "9");

	generaltree gt0(g0);
	ubt0 = g0;

	generaltree gtx(gt0);
	gtx.addEdge("2", "4");

	graph gx(g0);
	gx.addEdge("8", "6");

	unrootedbtree b1(g0);
	unrootedbtree b2(gt0);
	unrootedbtree b3(ubt0);

	unrootedbtree b4;
	unrootedbtree b5;

	try
	{
		b4 = gx;
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RAE:  "<<e.what()<<"\tb4 initialization failed"<<endl;
	}

	try
	{
		b5 = gtx;
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RAE:  "<<e.what()<<"\tb5 initialization failed"<<endl;
	}

	try
	{
		unrootedbtree b6(gx);
		cout<<"b6 initialized successfully."<<endl;
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RAE:  "<<e.what()<<"\tb6 initialization failed"<<endl;
	}

	try
	{
		unrootedbtree b7(gtx);
		cout<<"b7 initialized successfully."<<endl;
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RAE:  "<<e.what()<<"\tb7 initialization failed"<<endl;
	}

	return 0;
}

int rbt_test(int argc, char** argv)
{
	cout<<"rbt_test:"<<endl;
	graph g0;
	g0.addEdge("1", "2");
	g0.addEdge("2", "3");
	g0.addEdge("2", "5");
	g0.addEdge("2", "5");
	g0.addEdge("6", "5");
	g0.addEdge("7", "5");
	g0.addEdge("7", "8");

	generaltree gt0(g0);

	rootedbtree rbt0(g0, "7");
	rootedbtree rbt1(gt0, "7");

	rootedbtree rbt9;
	rbt9 = rbt0;

	try
	{
		rootedbtree rbt2(g0, "5");
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RAE:  "<<e.what()<<"\trbt2 initialization failure"<<endl;
	}

	try
	{
		rootedbtree rbt3(g0, "5");
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RAE:  "<<e.what()<<"\trbt3 initialization failure"<<endl;
	}

	graph gx0(g0);
	gx0.addEdge("4", "2");
	try
	{
		rootedbtree rbt4(gx0, "7");
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RAE:  "<<e.what()<<"\trbt4 initialization failure"<<endl;
	}

	graph gx1(g0);
	gx1.addEdge("6", "8");
	try
	{
		rootedbtree rbt5(gx1, "7");
	}
	catch (RecoverableAssertException &e)
	{
		cout<<"RAE:  "<<e.what()<<"\trbt5 initialization failure"<<endl;
	}

	return 0;
}

class outGtf : public graphTraversalFunction
{
public:
	virtual void visit(edge* ep, node* np)
	{
		cout<<ep->getId()<<"\t("<<np->getId()<<")"<<endl;
	}
};

int bfs_test(int argc, char** argv)
{
	cout<<"bfs_test:"<<endl;
	graph g0;
	g0.addEdge("1", "2");
	g0.addEdge("2", "3");
	g0.addEdge("2", "4");
	g0.addEdge("2", "5");
	g0.addEdge("6", "5");
	g0.addEdge("7", "5");
	g0.addEdge("7", "8");
	g0.addEdge("6", "8");

	outGtf outGtf;
	bfs(g0, "5", &outGtf);
	cout<<endl;
	rbfs(g0, "5", &outGtf);
	cout<<endl;

	cout<<g0<<endl;
	g0.removeEdge("6","8");
	g0.removeNode("4");
	cout<<g0<<endl;

	rootedbtree rbt(g0, "7");
	bfs(rbt, rbt.getRootId(), &outGtf);
	cout<<endl;
	rbfs(rbt, rbt.getRootId(), &outGtf);
	cout<<endl;

	return 0;
}

int rooted_phylo(int argc, char** argv)
{
	string yang;
	ifstream fin;

	cout<<"opening yang.tree"<<endl;
	reopen(fin, "data/yang.tree");
	fin>>yang;

	rootedbtree rbt;
	treeFromNewickString(yang, rbt);

	reopen(fin, "data/rate.mat");
    gslmat rate( fin );

    reopen(fin, "data/label.mask.mat");
    gslmat label_mask( fin );

    reopen(fin, "data/dwell.mask.mat");
	gslvec dwell_mask( fin );

	reopen(fin, "data/data.yang");
	gslmat yang_data( fin );

	gslvec stat_dist;
	precompute(&rbt, rate, label_mask, dwell_mask, stat_dist);

	fgval_vec_compute(&rbt, rbt.getRootId(), stat_dist, yang_data);
	stochmapcompute(rbt);

	cout<<stat_dist<<endl;
	cout<<"--------------------------------"<<endl;

	gslvec result_v;
	result_v.setsize(17);
	smapextract(rbt, result_v);
	cout<<result_v<<endl;

	cout<<"--------------------------------"<<endl;
	gslmat result_m;

    gslmat state_data(4, 5);
	for(int c=0; c<5; c++)
	{
		gsl_vector_view gv = gsl_matrix_column(state_data, c);
		int pos = gsl_matrix_get(yang_data, 0, c) - 1;
		gsl_vector_set(&gv.vector, pos, 1.0);
	}
	cout<<state_data<<endl;
	phylo_onemap(rbt, rate, label_mask, dwell_mask, state_data, result_m);
	cout<<result_m<<endl;

	cout<<"--------------------------------"<<endl;

	//Simulation from root w/ stationary distribution
    gslmat result;
    ofstream fout("data/sim4.data.mat");
    sampletips(rbt, rate, stat_dist, 10, 0, result);
    fout<<"# "<<result.nrow()<<" "<<result.ncol()<<endl;
    for(int r=0; r<result.nrow(); r++)
    {
    	for(int c=0; c<result.ncol(); c++)
    	{
    		double val = gsl_matrix_get(result, r, c);
    		int ival = floor(val + 1e-6);

    		fout<<ival;
    		if(c<result.ncol() - 1)
    			fout<<",";
    	}
    	if(r<result.nrow() - 1)
    		fout<<endl;
    }
    fout.close();

	cout<<"--------------------------------"<<endl;

	/* Compute stochastic mapping for multiple sites */
	reopen(fin, "data/sim.data2.mat");
	gslmat sim_data( fin );
	phylo_multimap(rbt, rate, label_mask, dwell_mask, sim_data, false, NULL, result_m);
	cout<<result_m<<endl;

	cout<<"--------------------------------"<<endl;

	reopen(fin, "data/sterling_rate.mat");
	gslmat sterling_rate(fin);

	string sterling_treestring;
	reopen(fin, "data/sterling.tree");
	fin >> sterling_treestring;
	rootedbtree sterling_rbt;
	treeFromNewickString(sterling_treestring, sterling_rbt);
	sampletips_sterling_ascertainment(sterling_rbt, sterling_rate, 10, 1, result_m);

    recycle(fout, "data/sim2.data.mat");
    fout<<"# "<<result_m.nrow()<<" "<<result_m.ncol()<<endl;
    for(int r=0; r<result_m.nrow(); r++)
    {
    	for(int c=0; c<result_m.ncol(); c++)
    	{
    		double val = gsl_matrix_get(result_m, r, c);
    		int ival = floor(val + 1e-6);

    		fout<<ival;
    		if(c<result_m.ncol() - 1)
    			fout<<",";
    	}
    	if(r<result_m.nrow() - 1)
    		fout<<endl;
    }
    fout.close();

	return 0;
}
