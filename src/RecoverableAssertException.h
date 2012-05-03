/*
 * AssertException.h
 *
 *  Created on: May 13, 2010
 *      Author: dnlennon
 */

#ifndef ASSERTEXCEPTION_H_
#define ASSERTEXCEPTION_H_

#include <string.h>

#include <string>
using namespace std;

class RecoverableAssertException
{
public:
	RecoverableAssertException(const char* fn, const char* file, int line)
	{
		char linestr[10];
		sprintf(linestr, "%d", line);

		pMessage = string(file) + ":" + string(linestr) + "\t" + string(fn) + "()";
	}

	const char* what() const { return pMessage.c_str(); }

private:
	string pMessage;
};

#endif /* ASSERTEXCEPTION_H_ */
