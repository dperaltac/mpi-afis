/**
 * \file    ParallelSlave.cpp
 * \author  Daniel Peralta <dperalta@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Implementation file for the ParallelSlave class.
 */

#include <cstdlib>
#include <cmath>
#include <sstream>
#include <fstream>

#include "ParallelSlave.h"
#include "Fingerprint.h"
#include "Functions.h"

using namespace std;

ParallelSlave::ParallelSlave() : ParallelHandler()
{
	begin = 0;
	end = 0;
}

ParallelSlave::~ParallelSlave()
{
	for(int i = 0; i < fp_list.rows(); ++i)
		for(int j = 0; j < fp_list.cols(); ++j)
			delete fp_list[i][j];
}

void ParallelSlave::loadInitializeDBFile(const vector<string> &filelists, unsigned int quality)
{
	loadDBFile(filelists, quality);
	initializeDatabase();
}

void ParallelSlave::loadDBFile(const vector<string> &filelists, unsigned int quality)
{
	vector<string> filenames;
	vector<string>::iterator fname;
	
	// There should be either:
	// - A single algorithm and any number of fingers
	// - A single finger and any number of algorithms
	// - One finger per algorithm
	if(filelists.size() > 1 && template_fp.size() > 1 &&
		filelists.size() != template_fp.size())
	{
		stringstream ss;
		ss << "ERROR: received " << template_fp.size() << " fingerprint types and " << filelists.size() << " template files (Process " << getProcessID() << ")";
		exit(ss.str(), -1);
	}
	
	numfingers = max(filelists.size(), template_fp.size());
	fp_list.clear();

	// Read the file list of each finger
	for(unsigned int j = 0; j < filelists.size(); ++j)
	{
		filenames = readFileByLines(filelists[j]);

		if(filenames.empty())
		{
			stringstream ss;
			ss << "ERROR: file " << filelists[j] << " is empty or could not be read (Process " << getProcessID() << ")";
			exit(ss.str(), -1);
		}
	
		// This is executed only once
		if(fp_list.empty())
			setListLimits(filenames.size(), filelists.size());

		// The process reads its chunk of files
		for(unsigned int i = 0; i < end-begin; ++i)
		{
			// Use the corresponding matcher
			fp_list[i][j] = newFingerprint(j);

			if(fp_list[i][j]->readFile(filenames[begin+i], quality))
				cerr << "WARNING: Process " << getProcessID() << ": Error when reading file " + *fname << endl;
		}
	}

	totalsize = fp_list.rows();
	calculateDBSize();
}




void ParallelSlave::sendStopSignal(const Score &match) const
{
	for(int k = 0; k < getProcesses(); ++k)
		if(k != getProcessID())
			MPI::COMM_WORLD.Isend(&found, 1, MPI::BOOL, k, TAG_OFFSET*iteration + TAG_FOUND);

	MPI::COMM_WORLD.Isend(&match, 1, mpi_pair_t, MASTER_ID, TAG_OFFSET*iteration + TAG_MATCH);
}


void ParallelSlave::sendNotFoundSignal() const
{
	MPI::COMM_WORLD.Isend(&found, 1, MPI::BOOL, MASTER_ID, TAG_OFFSET*iteration + TAG_FOUND);
}



void ParallelSlave::setListLimits(unsigned int num_files, unsigned int num_fingers)
{
	// Calculates the number of files that will be read by each process
	int chunk_size = ceil(num_files / (getProcesses() - 1.0));

	begin = (getProcessID()-1) * chunk_size;
	end = begin + chunk_size;

	if(end > num_files)
		end = num_files;

	// This is useful when there are more processes than files
	if(begin >= end)
		return;

	fp_list.resize(end-begin, num_fingers);
}


void ParallelSlave::initializeDatabase()
{
	// Fingerprint preprocessing
	#pragma omp parallel for schedule(dynamic) num_threads(omp_get_num_procs())
	for(int i = 0; i < fp_list.rows(); ++i)
		for(int j = 0; j < fp_list.cols(); ++j)
			fp_list[i][j]->initialize();
}


void ParallelSlave::run(vector<Fingerprint *> &vinput, vector< Score > &v)
{
	Score score;
	vector<float> fscore(getNumFingers());

	if(vinput.size() != getNumFingers())
	{
		cerr << "ERROR: the number of input fingers (" << vinput.size() << ") is different from the number of template fingers (" << getNumFingers() << ")" << endl;
		sendResultsToMaster();
		return;
	}

	iteration++;

	found_matches.clear();
	freeRequests();

	if(stop == STOP_YES)
	{
		waitForStop();
		found_matches.push_back(Score());
	}
	else if(stop == STOP_MAX)
		found_matches.push_back(Score());
	else if(stop == STOP_ALL)
		found_matches.resize(fp_list.rows());
	else if(stop == STOP_RANKING)
		found_matches.resize(ranking);

	avgmatchingtime = 0.0;

	#pragma omp parallel for schedule(dynamic) firstprivate(score, fscore)
	for(int i = 0; i < fp_list.rows(); ++i)
		if(!(stop == STOP_YES && (found || found_other)))
		{
			score.setId( fp_list[i][0]->getId() );
			score.setRelIndex(i);
			score.setAbsIndex(i+begin);
			score.setProcess(getProcessID());

			double tinicio = MPI::Wtime();

			for(int j = 0; j < fp_list.cols(); ++j)
				fscore[j] = fp_list[i][j]->match(*vinput[j]);

			score.setScore(aggregateScores(fscore));

			tinicio = MPI::Wtime() - tinicio;
			
			#pragma omp atomic
			avgmatchingtime += tinicio;

			// If the stop criterion is MAX, check if the found score is better than the best found score of the process
			if(stop == STOP_MAX && Fingerprint::betterOrEqual(score.getScore(), threshold))
			{
				#pragma omp critical
				if(Fingerprint::better(score, found_matches[0]))
					found_matches[0] = score;
			}

			// If the stop criteria is YES, check if another process has sent a stop signal, and if the score is better than the threshold send the stop signal
			else if(stop == STOP_YES && Fingerprint::betterOrEqual(score.getScore(), threshold))
			{
				#pragma omp critical
				if(!checkStop())
				{
					found_matches[0] = score;
					found = true;

					// Send the stop signal to all processes
					sendStopSignal(found_matches[0]);
				}
			}

			// If the stop criterion is ALL, all scores are sent to the master
			else if(stop == STOP_ALL)
				found_matches[i] = score;

			// If the stop criteria is NO, and the score is better than the threshold, add the matching to the list
			else if(stop == STOP_NO && Fingerprint::betterOrEqual(score.getScore(), threshold))
			{
				#pragma omp critical
				found_matches.push_back(score);
			}
			else if(stop == STOP_RANKING)
			{
				#pragma omp critical
				if(Fingerprint::better(score, found_matches.back()))
				{
					int j = found_matches.size()-2;
					while(j >= 0 && Fingerprint::better(score, found_matches[j]))
					{
						found_matches[j+1] = found_matches[j];
						--j;
					}

					if(Fingerprint::better(score, found_matches[j+1]))
						found_matches[j+1] = score;
				}
			}
		}

	if(stop == STOP_RANKING)
	{
		unsigned int pos = found_matches.size()-1;

		while(Score::isnull(found_matches[pos]) && pos > 0)
			--pos;

		if(pos != found_matches.size())
			found_matches.erase(found_matches.begin()+pos+1, found_matches.end());
	}

	sendResultsToMaster();

	v = found_matches;
}

void ParallelSlave::sendResultsToMaster()
{
	MPI::COMM_WORLD.Reduce(&avgmatchingtime, 0, 1, MPI::DOUBLE, MPI::SUM, MASTER_ID);

	if(stop == STOP_NO || stop == STOP_ALL || stop == STOP_RANKING)
		sendMatches();
	else if(stop == STOP_MAX)
		sendSingleMatches();
	else if(stop == STOP_YES && !checkStop() && !found && !found_other)
		sendNotFoundSignal();
}



void ParallelSlave::runFirst(vector<Fingerprint *> &vinput, vector< Score > &v)
{
	Score score;
	vector<float> fscore(getNumFingers());
	int local_ranking = (int)ceil(((float)ranking)/(getProcesses()-1));

	iteration++;

	v.clear();
	freeRequests();

	if(stop == STOP_RANKING)
		v.resize(local_ranking, Score());
	else if(stop == STOP_MAX)
		v.push_back(Score());
	
	#pragma omp parallel for schedule(dynamic) firstprivate(score, fscore) shared(v, vinput)
	for(int i = 0; i < fp_list.rows(); ++i)
	{
		for(int j = 0; j < fp_list.cols(); ++j)
			fscore[j] = fp_list[i][j]->match(*vinput[j]);

		score.setScore(aggregateScores(fscore));

		if(stop == STOP_MAX)
		{
			#pragma omp critical
			if(Fingerprint::better(score, v[0]))
			{
				score.setId( fp_list[i][0]->getId() );
				score.setRelIndex(i);
				v[0] = score;
			}
		}
		else if(stop == STOP_NO && Fingerprint::betterOrEqual(score.getScore(), threshold))
		{
			score.setId( fp_list[i][0]->getId() );
			score.setRelIndex(i);

			#pragma omp critical
			v.push_back(score);
		}
		else if(stop == STOP_RANKING)
		{
			#pragma omp critical
			if(Fingerprint::better(score, v.back()))
			{
				score.setId( fp_list[i][0]->getId() );
				score.setRelIndex(i);

				unsigned int pos = binaryPosSearch(&v[0], v.size(), score, true);
				v.insert(v.begin()+pos, score);
				v.pop_back();
			}
		}
	}

	if(stop == STOP_RANKING)
	{
		unsigned int pos = v.size();

		while(pos > 0 && Score::isnull(v[pos-1]))
			--pos;

		if(pos != v.size())
			v.erase(v.begin()+pos, v.end());
	}
}



void ParallelSlave::runSecond(vector<Fingerprint *> &vinput, const std::vector< Score > &candidates, std::vector< Score > &v)
{
	Score score;
	vector<float> fscore(getNumFingers());

	iteration++;

	freeRequests();

	found_matches.clear();
	found_matches.push_back(Score());

	avgmatchingtime = 0.0;

	#pragma omp parallel for schedule(dynamic) firstprivate(score, fscore)
	for(unsigned int i = 0; i < candidates.size(); ++i)
	{
		score = candidates[i];
		double tinicio = MPI::Wtime();

		for(int j = 0; j < fp_list.cols(); ++j)
			fscore[j] = fp_list[candidates[i].getRelIndex()][j]->match(*vinput[j]);

		score.setScore(aggregateScores(fscore));

		tinicio = MPI::Wtime() - tinicio;
#pragma omp atomic
		avgmatchingtime += tinicio;

		// Check if the found score is better than the best found score of the process
		#pragma omp critical
		if(Fingerprint::better(score, found_matches[0]) && Fingerprint::better(score.getScore(), threshold))
			found_matches[0] = score;
	}

	sendResultsToMaster();

	v = found_matches;
}


void ParallelSlave::sendMatches() const
{
	int self_size = found_matches.size();

	// Gathering the number of matches of each process
	MPI::COMM_WORLD.Gather(&self_size, 1, MPI::INT, 0, 0, MPI::INT, MASTER_ID);

	// Gathering all the matches
	MPI::COMM_WORLD.Gatherv(&found_matches[0], self_size, mpi_pair_t, 0, 0, 0, mpi_pair_t, MASTER_ID);
}


void ParallelSlave::sendSingleMatches() const
{
	if(found_matches.size() != 1)
	{
		stringstream ss;
		ss << found_matches.size();
		exit(string("ERROR: in the MAX mode, the number of results must be 1. Instead, it is ") + ss.str(), -1);
	}

	// Sending the match to the master process
	MPI::COMM_WORLD.Gather(&found_matches[0], 1, mpi_pair_t, 0, 0, mpi_pair_t, MASTER_ID);
}




void ParallelSlave::waitForStop()
{
	// Clear the structures from possible preceding runs
	freeRequests();
	stop_requests.push_back(MPI::COMM_WORLD.Irecv(&found_other, 1, MPI::BOOL, MPI_ANY_SOURCE, TAG_OFFSET*iteration + TAG_FOUND));
}
