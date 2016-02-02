/**
 * \file    ParallelHandler.cpp
 * \author  Daniel Peralta <daninep@correo.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Implementation file for the ParallelHandler class.
 */

#include <cstdlib>
#include <cmath>
#include <sstream>
#include "mpi.h"

#include "Fingerprint.h"
#include "MCC.h"
#include "FingerprintJiang.h"
#include "ParallelHandler.h"
#include "ParallelMaster.h"
#include "ParallelSlave.h"
#include "Score.h"
#include "Functions.h"

using namespace std;

ParallelHandler::ParallelHandler()
{
	stop = STOP_MAX;
	ranking = 1;
	threshold = 0;
	iteration = 0;
	explored = 0;
	totalsize = 0;
	template_fp = vector<fingerprint_t>(1,FP_UNKNOWN);
	mpi_pair_t = Score::getDatatype();
	fusion = FUSION_SUM;
}

ParallelHandler::~ParallelHandler()
{
	freeRequests();
	(mpi_pair_t.Free)();
}


ParallelHandler * ParallelHandler::getHandler(int threads, bool classification, int level)
{
	if(!MPI::Is_initialized())
	{
		if(MPI::Init_thread(level) < level)
			exit("The parallelism level is not high enough. At least MPI_THREAD_MULTIPLE is required.", -1);
	}

	omp_set_num_threads(threads);

	if(getProcessID() == MASTER_ID)
		return new ParallelMaster();
	else
		return new ParallelSlave();
}


void ParallelHandler::exit(const string & message, int status)
{
	cerr << message << endl;

	if(MPI::Is_initialized())
		MPI::Finalize();

	std::exit(status);
}


void ParallelHandler::freeRequests()
{
	for(vector<MPI::Request>::iterator i = stop_requests.begin(); i != stop_requests.end(); ++i)
		if(*i != MPI::REQUEST_NULL)
		{
			i->Cancel();
			(i->Free)();
		}

	stop_requests.clear();
}


// void ParallelHandler::calculatePenetration()
// {
// 	int penetration[2] = {explored, totalsize};
// 	int totalpenetration[2*getProcesses()];
//
// 	// Gathering the number of matches of each process
// 	MPI::COMM_WORLD.Gather(penetration, 2, MPI::INT, totalpenetration, 2, MPI::INT, MASTER_ID);
//
// 	if(isMaster())
// 	{
// 		explored = 0;
// 		totalsize = 0;
//
// 		for(int i = 2; i < 2*getProcesses(); i+=2)
// 		{
// 			explored += totalpenetration[i];
// 			totalsize += totalpenetration[i+1];
// 		}
// 	}
// }




Fingerprint * ParallelHandler::newFingerprint(int finger) const
{
	finger = min<int>(finger, template_fp.size()-1);

	if(template_fp[finger] == FP_UNKNOWN)
		exit("newFingerprint(): the fingerprint type has not been set yet (use function setFingerprintType())", -1);
	else if(template_fp[finger] == FP_JIANG)
		return new FingerprintJiang();
	else if(template_fp[finger] == FP_MCC)
		return new MCC();

	return 0;
}

void ParallelHandler::setFingerprintType(const vector<ParallelHandler::fingerprint_t> &fp, int argc, char *argv[])
{
	template_fp.resize(fp.size());

	for(unsigned int i = 0; i < fp.size(); ++i)
		setFingerprintType(fp[i], argc, argv, i);
}

void ParallelHandler::setFingerprintType(ParallelHandler::fingerprint_t fp, int argc, char *argv[], int finger)
{
	template_fp[finger] = fp;

	if(fp == FP_UNKNOWN)
		exit("setFingerprintType(): the fingerprint type is invalid", -1);
	else if(fp == FP_JIANG)
		FingerprintJiang::configureAlgorithm(argc, argv);
	else if(fp == FP_MCC)
		MCC::configureAlgorithm(argc, argv);
}

void ParallelHandler::calculateDBSize()
{
	int sum = 0;

	MPI::COMM_WORLD.Reduce(&totalsize, &sum, 1, MPI::INT, MPI::SUM, MASTER_ID);

	if(isMaster())
		totalsize = sum;
}

void ParallelHandler::calculateExplored()
{
	int sum = 0;

	MPI::COMM_WORLD.Reduce(&explored, &sum, 1, MPI::INT, MPI::SUM, MASTER_ID);

	if(isMaster())
		explored = sum;
}


float ParallelHandler::aggregateScores(const vector<float> &fscore) const
{
	float acum = 0.0;
	vector<float> valid_scores;

	for(vector<float>::const_iterator i = fscore.begin(); i != fscore.end(); ++i)
		if(*i != -1)
			valid_scores.push_back(*i);

	if(valid_scores.empty())
		return -1;
	else if(valid_scores.size() == 1)
		return valid_scores[0];

	switch(fusion)
	{
		case FUSION_SUM:
			acum = 0.0;
			for(vector<float>::const_iterator i = valid_scores.begin(); i != valid_scores.end(); ++i)
				acum += *i;
			acum /= fscore.size();

			break;
		case FUSION_PROD:
			acum = 1.0;
			for(vector<float>::const_iterator i = valid_scores.begin(); i != valid_scores.end(); ++i)
				acum *= *i;

			break;
		case FUSION_MAX:
			acum = -1;
			for(vector<float>::const_iterator i = valid_scores.begin(); i != valid_scores.end(); ++i)
				if(*i > acum)
					acum = *i;

			break;
		case FUSION_MIN:
			acum = 100000000000;
			for(vector<float>::const_iterator i = valid_scores.begin(); i != valid_scores.end(); ++i)
				if(*i < acum)
					acum = *i;

			break;
		default:
			acum = -1;
	}

	return acum;
}


void ParallelHandler::loadDBFile(const vector<string> &filelists, const vector<string> &classfiles, unsigned int quality)
{
	if(filelists.size() == 1 && classfiles.size() == 1)
		loadDBFile(filelists[0], classfiles[0], quality);
	else if(classfiles.empty())
		loadDBFile(filelists, quality);
	else
	{
		stringstream ss;
		ss << "ERROR: received " << filelists.size() << " template files and " << classfiles.size() << " class files (Process " << getProcessID() << ")";
		exit(ss.str(), -1);
	}
}


void ParallelHandler::loadInitializeDBFile(const vector<string> &filelists, const vector<string> &classfiles, unsigned int quality)
{
	if(filelists.size() == 1 && classfiles.size() == 1)
		loadInitializeDBFile(filelists[0], classfiles[0], quality);
	else if(classfiles.empty())
		loadInitializeDBFile(filelists, quality);
	else
	{
		stringstream ss;
		ss << "ERROR: received " << filelists.size() << " template files and " << classfiles.size() << " class files (Process " << getProcessID() << ")";
		exit(ss.str(), -1);
	}
}
