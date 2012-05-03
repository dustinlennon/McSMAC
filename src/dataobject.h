/*
 * dataobject.h
 *
 *  Created on: May 4, 2010
 *      Author: dnlennon
 */

#ifndef DATAOBJECT_H_
#define DATAOBJECT_H_

#include <iostream>
using namespace std;

// This is a simple dataObject with corresponding auto pointer.

typedef class dataObject dataObject;
typedef class DOP DOP;
typedef class DOPuser DOPuser;

class dataObject
{
	friend class DOP;
	friend ostream& operator<<(ostream& o, const dataObject& dor);

public:

	dataObject(string id) : myid(id) { }
	virtual ~dataObject() {}

	string id(void) const { return myid; }

private:
	virtual dataObject* deepCopy()
	{
		dataObject* copy = new dataObject(myid);
		return copy;
	};

protected:
	string myid;
};

ostream& operator<<(ostream& o, const dataObject& dor);

/* The DOP class is an auto pointer for the dataObject and derived classes.  It allows a user
   to associate data with nodes and edges without worrying about cleaning up that data.
   */
class DOP
{
public:
	DOP()
	{
		refct = NULL;
		mydop = NULL;
	}

	DOP(dataObject* dop)
	{
		if(dop == NULL)
		{
			mydop = NULL;
			refct = NULL;
			return;
		}

		mydop = dop;
		refct = new int;
		*refct = 1;
	}

	DOP(const DOP& dop)
	{
		if(dop.mydop == NULL)
		{
			mydop = NULL;
			refct = NULL;
			return;
		}
		mydop = dop.mydop;
		refct = dop.refct;
		(*refct)++;
	}

	const DOP operator=(const DOP& dop)
	{
		if(mydop == dop.mydop)
			return dop;

		if(mydop != NULL)
		{
			(*refct)--;
			if(*refct == 0)
			{
				delete refct;
				delete mydop;
			}
		}

		if(dop.mydop == NULL)
		{
			mydop = NULL;
			refct = NULL;
			return NULL;
		}

		refct = dop.refct;
		mydop = dop.mydop;
		(*refct)++;

		return *this;
	}

	~DOP()
	{
		if(mydop == NULL)
			return;

		(*refct)--;
		if(*refct == 0)
		{
			delete refct;
			delete mydop;
		}
	}

	void makeOwner(void)
	{
		if(mydop == NULL)
			return;

		(*refct)--;
		if(*refct == 0)
		{
			(*refct)++;
			return;
		}

		mydop = mydop->deepCopy();
		refct = new int;
		*refct = 1;
	}

	dataObject& operator*()
	{
		return *mydop;
	}

	dataObject* ptr()
	{
		return mydop;
	}

private:
	dataObject* mydop;
	int *refct;
};

class DOPuser
{
public:
	DOPuser() { }
	DOPuser(DOP d) : mydop(d) { }

	DOP getDOP(void) const
	{
		return mydop;
	}

	void setDOP(DOP d)
	{
		mydop = d;
	}

	void makeOwnerDOP(void)
	{
		mydop.makeOwner();
	}

protected:
	DOP mydop;
};

#endif /* DATAOBJECT_H_ */
