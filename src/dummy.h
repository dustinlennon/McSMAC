/*
 * dummy.h
 *
 *  Created on: May 4, 2010
 *      Author: dnlennon
 */

#ifndef DUMMY_H_
#define DUMMY_H_

#include <iostream>

class dummy
{
friend ostream& operator<<(ostream& o, const dummy& d);

public:
	dummy(void) { myval = 0; }
	dummy(int val) : myval(val) { }

	dummy(const dummy& d)
	{
		cout<<"dummy copy constructor called"<<endl;
		myval = d.myval;
	}

	dummy& operator=(const dummy& d)
	{
		cout<<"dummy operator= called"<<endl;
		myval = d.myval;
		return *this;
	}

private:
	int myval;
};

ostream& operator<<(ostream& o, const dummy& d)
{
	o<<d.myval;
	return o;
}


#endif /* DUMMY_H_ */
