/**
 * \file    DPDDFF.cpp
 * \author  Daniel Peralta <dperalta@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Main file for the dual phase identification system DPD-DFF
 */

#include <iostream>
#include <fstream>
#include "mpi.h"
#include <string.h>
#include <cstdlib>
#include <cmath>
#include <unistd.h>

#include "ParallelHandler.h"
#include "ParallelSlave.h"
#include "ParallelMaster.h"
#include "IOHandler.h"
#include "Score.h"

using namespace std;

int main(int argc, char * argv[])
{
	// Variables
	unsigned int num_candidates=0, total_candidates;
	double begin_time = MPI::Wtime();
	vector<double> times(3);
	string cinput[2];
	int numfps_phase1, numfps_phase2;

	// Structures
	ParallelHandler *parallel_handler1;
	ParallelHandler *parallel_handler2;
	vector< Score > candidates, matches;
	vector< vector< float > > scoresmatrix;
	
	IOHandler iohandler(argc, argv);
	ParallelHandler::variant_t variant = iohandler.getVariant();

	// MPI initialization
	parallel_handler1 = ParallelHandler::getHandler(iohandler.getThreads(), false, MPI_THREAD_MULTIPLE);
	parallel_handler2 = ParallelHandler::getHandler(iohandler.getThreads(), false, MPI_THREAD_MULTIPLE);

	// Set matching algorithm
	parallel_handler1->setFusion(iohandler.getFusionType());
	parallel_handler2->setFusion(iohandler.getFusionType());
	parallel_handler1->setFingerprintType(iohandler.getFpType(0), argc, argv);
	parallel_handler2->setFingerprintType(iohandler.getFpType(1), argc, argv);

	// Set auxiliary variables
	if(variant == ParallelHandler::DPDDFF_SS || variant == ParallelHandler::DPDDFF_SD)
		numfps_phase1 = 1;
	else
		numfps_phase1 = 2;
	
	if(variant == ParallelHandler::DPDDFF_SS || variant == ParallelHandler::DPDDFF_DS)
		numfps_phase2 = 1;
	else
		numfps_phase2 = 2;
	

	// Load databases
	if(numfps_phase1 == 1)
		parallel_handler1->loadInitializeDBFile(iohandler.getTemplateFile(0), iohandler.getQuality());
	else
		parallel_handler1->loadInitializeDBFile(iohandler.getTemplateFiles(), iohandler.getQuality());
	
	if(numfps_phase2 == 1)
		parallel_handler2->loadInitializeDBFile(iohandler.getTemplateFile(1), iohandler.getQuality());
	else
		parallel_handler2->loadInitializeDBFile(iohandler.getTemplateFiles(), iohandler.getQuality());

	
	// Initial output
	if(parallel_handler1->isMaster())
		iohandler.printInitialOutput();

	// Set parameters
	parallel_handler1->setStopMode(iohandler.getStopType());
	parallel_handler1->setRanking(iohandler.getRanking());
	parallel_handler1->setThreshold(iohandler.getThreshold(0));
	parallel_handler2->setStopMode(ParallelHandler::STOP_MAX);
	parallel_handler2->setThreshold(iohandler.getThreshold(1));

	// Measure initialization time
	MPI::COMM_WORLD.Barrier();
	times[0] = MPI::Wtime() - begin_time;

	// Read input fingerprints
	iohandler.readFileName(0, cinput[0], parallel_handler1);
	iohandler.readFileName(1, cinput[1], parallel_handler2);

	while(cinput[0] != "exit")
	{
		for(unsigned int j = 1; j < times.size(); j++)
			times[j] = 0.0;

		vector<Fingerprint *> newF1(numfps_phase1, NULL);
		vector<Fingerprint *> newF2(numfps_phase2, NULL);

		// Read input fingerprints
		if(!parallel_handler1->isMaster())
		{
			for(unsigned int i = 0; i < newF1.size(); ++i)
			{
				newF1[i] = parallel_handler1->newFingerprint();

				if(newF1[i]->readFile(cinput[i]))
					parallel_handler1->exit(string("Error when reading file ") + cinput[i], -1);
			}
			
			for(unsigned int i = 0; i < newF2.size(); ++i)
			{
				// The index ensures that if only one fingeprint is used in the second phase,
				// it will be the second one
				unsigned int pos = i - newF2.size() + 2;
				
				newF2[i] = parallel_handler2->newFingerprint();

				if(newF2[i]->readFile(cinput[pos]))
					parallel_handler2->exit(string("Error when reading file ") + cinput[pos], -1);
			}
		}

		// time[1]: measuring the candidate selection with the first fingerprint
		MPI::COMM_WORLD.Barrier();
		begin_time = MPI::Wtime();

		if(!parallel_handler1->isMaster())
		{
			for(unsigned int i = 0; i < newF1.size(); ++i)
				newF1[i]->initialize();
			((ParallelSlave *)parallel_handler1)->runFirst(newF1, candidates);
		}

// 		time[2]: measuring the final identification with the second fingerprint
		MPI::COMM_WORLD.Barrier();
		times[1] += MPI::Wtime() - begin_time;
		begin_time = MPI::Wtime();

		if(!parallel_handler2->isMaster())
			for(unsigned int i = 0; i < newF2.size(); ++i)
				newF2[i]->initialize();

// 		if(!parallel_handler1->isMaster())
// 			iohandler.printIdentificationOutput(cinput[0], candidates, -1, times, parallel_handler2->getIteration(), candidates.size());
                
		parallel_handler2->runSecond(newF2, candidates, matches);

		times[2] += MPI::Wtime() - begin_time;
		begin_time = MPI::Wtime();

		num_candidates = candidates.size();

		MPI::COMM_WORLD.Reduce(&num_candidates, &total_candidates, 1, MPI::INT, MPI::SUM, parallel_handler2->MASTER_ID);

		// Clear memory
		for(unsigned int i = 0; i < newF1.size(); ++i)
			delete newF1[i];
		for(unsigned int i = 0; i < newF2.size(); ++i)
			delete newF2[i];

		// Output
		if(parallel_handler1->isMaster())
			iohandler.printIdentificationOutput(cinput[0], matches, -1, times, parallel_handler2->getIteration(), total_candidates, parallel_handler1->getPenetrationRate());

		// Read next fingeprint pair
		iohandler.readFileName(0, cinput[0], parallel_handler1);
		iohandler.readFileName(1, cinput[1], parallel_handler2);
	}

	if(parallel_handler2->isMaster())
		iohandler.printFinalOutput(parallel_handler2);

	MPI::COMM_WORLD.Barrier();

	delete parallel_handler1;
	delete parallel_handler2;

	MPI::Finalize();

	return 0;
}
