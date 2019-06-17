/**
 * \file    main.cpp
 * \author  Salvador García <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Main file for the MCC matching algorithm.
 */

#include <iostream>
#include "MCC.h"
#include <string.h>

using namespace std;

int main(int argc, char * argv[])
{
	MCC a1;
	MCC a2;

	if(argc < 3){
			cout << "Usage: MCC <fingerprint1> <fingerprint2> -N {8|16} -C {LSS|LSSR|LSA|LSAR|LGS|NHS} [-H] [-B]" << endl;
			return 0;
	}

	MCC::configureAlgorithm(argc, argv);

	if(a1.readFile(argv[1])!=0){
			cout << "Error opening fingerprint files: " << argv[1] << endl;
			return 0;
	}

	if(a2.readFile(argv[2])!=0){
			cout << "Error opening fingerprint files: " << argv[2] << endl;
			return 0;
	}

	a1.initialize();
	a2.initialize();
	
	cout << "First fingerprint: " << a1 << endl;
	cout << "Second fingerprint: " << a2 << endl;

	cout << a1.match(a2) << endl;

	return 0;
}
