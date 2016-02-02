/**
 * \file    main.cpp
 * \author  Salvador García <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Main file for the Jiang matching algorithm.
 */

#include <iostream>
#include "FingerprintJiang.h"

using namespace std;

int main(int argc, char * argv[]){

    FingerprintJiang a1;
    FingerprintJiang a2;
    int state;

    Fingerprint::configureAlgorithm(argc, argv);

    if(argc != 3){
        cout << "Usage: jiangMatching <fingerprint1> <fingerprint2>" << endl;
        return 0;
    }

    state=a1.readFile(argv[1]);

    if(state!=0){
        cout << "Error opening fingerprint files: " << argv[1];
        return 0;
    }

    state=a2.readFile(argv[2]);

    if(state!=0){
        cout << "Error opening fingerprint files: " << argv[2];
        return 0;
    }

    a1.initialize();
    a2.initialize();

    cout << a1.match(a2) << endl;

    return 0;
}
