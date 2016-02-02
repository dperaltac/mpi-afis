/**
 * \file    genericMatching.cpp
 * \author  Daniel Peralta <dperalta@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Main file for a parallel application of any matching algorithm.
 */

#include <iostream>
#include <fstream>
#include "mpi.h"
#include <string.h>
#include <cstdlib>
#include <cmath>

#include "ParallelHandler.h"
#include "ParallelSlave.h"
#include "ParallelMaster.h"
#include "IOHandler.h"
#include "Score.h"

using namespace std;

int main(int argc, char * argv[])
{
	// Variables
	double begin_time = MPI::Wtime();
	vector<double> times(3);
	vector<string> cinput;
	double avgmatchingtime;

	// Structures
	ParallelHandler *parallel_handler;
	vector< Score > matches;
	IOHandler iohandler(argc, argv);
	vector< ParallelHandler::fingerprint_t > fp_types = iohandler.getFpTypes();
	vector< Matrix<double> > feature_matrices;

	// MPI initialization
	parallel_handler = ParallelHandler::getHandler(iohandler.getThreads(), iohandler.getClassification(), MPI_THREAD_MULTIPLE);

	parallel_handler->setFusion(iohandler.getFusionType());
	parallel_handler->setFingerprintType(fp_types, argc, argv);

	parallel_handler->loadDBFile(iohandler.getTemplateFiles(), iohandler.getClassFiles(), iohandler.getQuality());

	if(parallel_handler->isMaster())
		iohandler.printInitialOutput();

	parallel_handler->setStopMode(iohandler.getStopType());
	parallel_handler->setRanking(iohandler.getRanking());
	parallel_handler->setThreshold(iohandler.getThreshold());

	// Time to load the database
	MPI::COMM_WORLD.Barrier();
	times[0] = MPI::Wtime() - begin_time;

	if(!parallel_handler->isMaster())
			((ParallelSlave *) parallel_handler)->initializeDatabase();

	// Time to initialize the database
	MPI::COMM_WORLD.Barrier();
	times[1] = MPI::Wtime() - begin_time;

	iohandler.readFileName(cinput, parallel_handler);

	while(cinput[0] != "exit")
	{
		MPI::COMM_WORLD.Barrier();
		begin_time = MPI::Wtime();

		vector<Fingerprint *> newF(parallel_handler->getNumFingers(), 0);

		for(unsigned int j = 2; j < times.size(); j++)
			times[j] = 0;

		// Read fingerprint
		if(!parallel_handler->isMaster())
		{
			for(unsigned int i = 0; i < newF.size(); ++i)
			{
				newF[i] = parallel_handler->newFingerprint(i);

				if(cinput.size() > 1 && newF[i]->readFile(cinput[i], iohandler.getQuality()))
					cerr << "Error when reading file " << cinput[i] << endl;
				else if(cinput.size() == 1 && newF[i]->readFile(cinput[0], iohandler.getQuality()))
					cerr << "Error when reading file " << cinput[0] << endl;
			}
		}

		if(!parallel_handler->isMaster())
			for(unsigned int i = 0; i < newF.size(); ++i)
				newF[i]->initialize();

		parallel_handler->run(newF, matches);

		avgmatchingtime = parallel_handler->getAvgMatchingTime();

		times[2] += MPI::Wtime() - begin_time;
		begin_time = MPI::Wtime();

		if(!parallel_handler->isMaster())
			for(unsigned int i = 0; i < newF.size(); ++i)
				delete newF[i];

		if(parallel_handler->isMaster())
			iohandler.printIdentificationOutput(cinput[0], matches, avgmatchingtime, times, parallel_handler->getIteration(), 1, parallel_handler->getPenetrationRate());

		iohandler.readFileName(cinput, parallel_handler);
	}

	if(parallel_handler->isMaster())
		iohandler.printFinalOutput(parallel_handler);

	MPI::COMM_WORLD.Barrier();

	delete parallel_handler;

	MPI::Finalize();

	return 0;
}
