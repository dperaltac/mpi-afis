/**
 * \file    ParallelMaster.cpp
 * \author  Daniel Peralta <daninep@correo.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Implementation file for the ParallelMaster class.
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "ParallelMaster.h"

using namespace std;

ParallelMaster::ParallelMaster() : ParallelHandler()
{
	found_others = 0;
}

ParallelMaster::~ParallelMaster()
{
	delete [] found_others;
}

void ParallelMaster::setBoundaries(unsigned int size) const
{
	int processes = getProcesses();
	unsigned int boundaries[2];
	int chunk_size = ceil(size / (double)processes);

	boundaries[1] = 0;

	for(int i = 0; i < processes; ++i)
	{
		boundaries[0] = boundaries[1];
		boundaries[1] = boundaries[0] + chunk_size;

		if(boundaries[1] > size)
			boundaries[1] = size;

		MPI::COMM_WORLD.Isend(boundaries, 2, MPI::INT, i, TAG_OFFSET*iteration + TAG_BOUNDARIES);
	}
}


void ParallelMaster::waitForStop()
{
	// Clear the structures from possible preceding runs
	freeRequests();

	for(int i = 0; i < getProcesses(); ++i)
		if(i != getProcessID())
			stop_requests.push_back(MPI::COMM_WORLD.Irecv(&found_others[i], 1, MPI::BOOL, i, TAG_OFFSET*iteration + TAG_FOUND));
}

void ParallelMaster::run(vector< Score > &v, unsigned int numfps)
{
	v.clear();
	iteration++;

	found = false;
	
	// Calculate the average matching time
	avgmatchingtime = 0.0;
	double tmp = 0.0;
	
	MPI::COMM_WORLD.Reduce(&tmp, &avgmatchingtime, 1, MPI::DOUBLE, MPI::SUM, MASTER_ID);
	avgmatchingtime /= numfps;
	
	if(stop == STOP_MAX)
		retrieveBestMatch(v);
	else if(stop == STOP_YES)
		retrieveMatchesStop(v);
	else if(stop == STOP_NO || stop == STOP_ALL)
		retrieveAllMatches(v);
	else if(stop == STOP_RANKING)
	{
		retrieveAllMatches(v);
		sort(v.begin(), v.end(), Score::better);
		v.resize(ranking);
	}
}


void ParallelMaster::retrieveMatchesStop(vector< Score > &v)
{
	Score match, best_match;
	int count = 0, source;
	int received;
	int index[getProcesses()-1];
	MPI::Status status[getProcesses()-1];
	
	v.clear();

	if(found_others == 0)
		found_others = new bool[getProcesses()];

	for(int i = 0; i < getProcesses(); ++i)
		found_others[i] = false;

	found = false;

	waitForStop();

	do
	{
		// Wait for a stop request
		received = MPI::Request::Waitsome(stop_requests.size(), &stop_requests[0], index, status);

		// Sums the processes that have concluded their search
		count += received;

		// Delete the obsolete requests
		for(int i = 0; i < received; i++)
		{
			source = status[i].Get_source();

			// If a slave process has found a match, receive it
			if(found_others[source])
			{
				found = true;

				MPI::COMM_WORLD.Recv(&match, 1, mpi_pair_t, source, TAG_OFFSET*iteration + TAG_MATCH);

				if(Fingerprint::better(match.getScore(), best_match.getScore()))
					best_match = match;
			}
			else
			{
				stop_requests.erase(stop_requests.begin() + index[i]);
				for(int j = i+1; j < received; j++)
					if(index[j] > index[i])
						index[j]--;
			}
		}

	} while(!found && count < getProcesses()-1);

	// The requests are not needed any more
	freeRequests();

	v.push_back(best_match);
}



void ParallelMaster::retrieveAllMatches(vector< Score > &v) const
{
	int processes = getProcesses();
	int displs[processes], sizes[processes];
	int self_size = 0;
	int total_size = 0;

// 	cout << "Master: gathering results" << endl;
	// Gathering the number of matches of each process
	MPI::COMM_WORLD.Gather(&self_size, 1, MPI::INT, sizes, 1, MPI::INT, MASTER_ID);

	// The master process calculates the displacements
	displs[0] = 0;

	for(int i = 1; i < processes; ++i)
		displs[i] = displs[i-1] + sizes[i-1];

	total_size = sizes[processes-1] + displs[processes-1];

	v.resize(total_size);

	// Gathering all the matches
	MPI::COMM_WORLD.Gatherv(&v[0], self_size, mpi_pair_t, &v[0], sizes, displs, mpi_pair_t, MASTER_ID);
}


void ParallelMaster::retrieveBestMatch(vector<Score> &v) const
{
	vector<Score> received(getProcesses());
	
	v.resize(1);

	// Gathering all the matches
	MPI::COMM_WORLD.Gather(&v[0], 1, mpi_pair_t, &received[0], 1, mpi_pair_t, MASTER_ID);

	// Eliminate the master position, which is not valid (anyway it's the null score)
	received.front() = received.back();
	received.pop_back();

	for(vector< Score >::const_iterator i = received.begin(); i != received.end(); ++i)
		if(i->getScore() > v[0].getScore())
			v[0] = *i;
}



void ParallelMaster::loadDBFile(const vector<string> &filenames, unsigned int quality)
{
	numfingers = max(filenames.size(), template_fp.size());
	calculateDBSize();
}

