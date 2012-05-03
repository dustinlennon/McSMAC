/*
 * dataobject.cc
 *
 *  Created on: May 5, 2010
 *      Author: dnlennon
 */

#include "dataobject.h"

#include <ostream>
using namespace std;

ostream& operator<<(ostream& o, const dataObject& dor)
{
	o<<dor.myid;
	return o;
}
