/**
 * \file    scoreMatrix.cpp
 * \author  Daniel Peralta <dperalta@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Main file for a parallel application that generates the score matrix for a matcher and a database.
 */

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "ParallelHandler.h"
#include "ParallelSlave.h"
#include "ParallelMaster.h"
#include "IOHandler.h"
#include "Score.h"

using namespace std;

void writeAccuracyMeasures(const Matrix<float> &scoresmatrix, const string &filename, const ParallelHandler *parallel_handler);

int main(int argc, char * argv[])
{
	// Variables
	vector<string> cinput;
	vector<Fingerprint *> newF;
	ifstream is;
	ofstream verif;
	double avgmatchingtime = 0;

	// Structures
	ParallelHandler *parallel_handler;
	vector< Score > matches;
	Matrix<float> scoresmatrix;
	IOHandler iohandler(argc, argv);

	// MPI initialization
	parallel_handler = ParallelHandler::getHandler(iohandler.getThreads(), MPI_THREAD_MULTIPLE);

	parallel_handler->setFingerprintType(iohandler.getFpTypes(), argc, argv);
	parallel_handler->loadDBFile(iohandler.getTemplateFiles(), iohandler.getQuality());
	
	newF.resize(parallel_handler->getNumFingers(), 0);

	if(parallel_handler->isMaster())
		scoresmatrix.resize(iohandler.getNumInputFingerprints(), parallel_handler->getDBSize());

	parallel_handler->setStopMode(ParallelHandler::STOP_ALL);

	iohandler.readFileName(cinput, parallel_handler);
// 		cout << parallel_handler->getProcessID() << " : Starting loop" << endl;

	while(cinput[0] != "exit")
	{
// 		cout << "Will read " << cinput << endl;
		if(!parallel_handler->isMaster())
		{
			for(unsigned int i = 0; i < newF.size(); ++i)
			{
				newF[i] = parallel_handler->newFingerprint(i);
			
				if(cinput.size() > 1 && newF[i]->readFile(cinput[i], iohandler.getQuality()))
					cerr << "Error when reading file " << cinput[i] << endl;
				else if(cinput.size() == 1 && newF[i]->readFile(cinput[0], iohandler.getQuality()))
					cerr << "Error when reading file " << cinput[0] << endl;
				else
					newF[i]->initialize();
			}
		}
// 		cout << "Read file " << cinput << endl;

		parallel_handler->run(newF, matches);

		avgmatchingtime += parallel_handler->getAvgMatchingTime();

		if(parallel_handler->isMaster())
		{
			for(unsigned int i = 0; i < matches.size(); i++)
				scoresmatrix[parallel_handler->getIteration()-1][i] = matches[i].getScore();

// 			cout << "Huella " << parallel_handler->getIteration() << " procesada" << endl;
// 			cout << "Huella " << cinput[0] << " procesada" << endl;
		}

		iohandler.readFileName(cinput, parallel_handler);
	}
	

	if(parallel_handler->isMaster())
	{
		cout << "Average matching time: " << avgmatchingtime/scoresmatrix.rows() << endl;
		writeAccuracyMeasures(scoresmatrix, iohandler.getOutputPath(), parallel_handler);
	}

	MPI::COMM_WORLD.Barrier();

	if(!parallel_handler->isMaster())
		for(unsigned int i = 0; i < newF.size(); ++i)
			delete newF[i];
		
	delete parallel_handler;

	MPI::Finalize();

	return 0;
}


void writeAccuracyMeasures(const Matrix<float> &scoresmatrix, const string &filename, const ParallelHandler *parallel_handler)
{
	unsigned int begin = 0, end = scoresmatrix.rows();
	const unsigned int SCORE_INTERVALS = 100;
	const unsigned int RANKS = end-begin;
	unsigned int takes = scoresmatrix.cols() / scoresmatrix.rows();
	unsigned int rows = end-begin, cols = (end-begin)*takes;
	unsigned int begcol = begin*takes;
	unsigned int endcol = begcol+cols;
// 	const unsigned int NUMCLIENTS = rows*takes;
// 	const unsigned int NUMIMPOSTORS = rows*(cols-takes);
	unsigned int validclients = rows;
	unsigned int validimpostors = cols;
	float interval_size;
	unsigned int pos, aux;
	unsigned int dist_clients  [SCORE_INTERVALS+1];
	unsigned int dist_impostors[SCORE_INTERVALS+1];
	unsigned int dist_frr      [SCORE_INTERVALS+1];
	unsigned int dist_far      [SCORE_INTERVALS+1];
	double ddist_clients  [SCORE_INTERVALS+1];
	double ddist_impostors[SCORE_INTERVALS+1];
	double ddist_frr      [SCORE_INTERVALS+1];
	double ddist_far      [SCORE_INTERVALS+1];
	vector<unsigned int> cmc    (RANKS, 0);
	vector<bool> nancli         (rows, true);
	vector<bool> nanimp         (cols, true);
	float scoremax = -1, scoremin = 10000;

	float eer = -1, far100 = -1, far1000 = -1, far10000 = -1, zerofar = -1, zerofrr = -1;
	unsigned int r100 = RANKS, r1000 = RANKS, r10000 = RANKS, r1 = RANKS;

	ofstream distcliimp((filename + "_distcliimp.dat").c_str());
	ofstream frrfar((filename + "_frrfar.dat").c_str());
	ofstream cmcfile((filename + "_cmc.dat").c_str());
	ofstream nanclifile((filename + "_nancli.dat").c_str());
	ofstream nanimpfile((filename + "_nanimp.dat").c_str());
	ofstream matrixfile((filename + "_matrix.dat").c_str());
	ofstream verif((filename + "_verif.dat").c_str());

	if(!distcliimp || !frrfar || !cmcfile || !nanclifile || !nanimpfile || !matrixfile || !verif)
		ParallelHandler::exit(string("Error when opening files ") + filename + "*");

	for(unsigned int i = 0; i <= SCORE_INTERVALS; ++i)
	{
		dist_clients[i] = 0;
		dist_impostors[i] = 0;
		dist_frr[i] = 0;
		dist_far[i] = 0;
	}

	for(unsigned int i = begin; i < end; ++i)
	{
		for(unsigned int j = begcol; j < endcol; ++j)
		{
			if(scoresmatrix[i][j] > scoremax)
				scoremax = scoresmatrix[i][j];
			if(scoresmatrix[i][j] >= 0 && scoresmatrix[i][j] < scoremin)
				scoremin = scoresmatrix[i][j];

			matrixfile << scoresmatrix[i][j] << " ";
		}
		matrixfile << endl;
	}

	matrixfile.close();

	interval_size = (scoremax - scoremin)/SCORE_INTERVALS;

	// Detect the fingerprints that are not valid for the matcher
	for(unsigned int i = begin; i < end; ++i)
	{
		for(unsigned int j = begcol; j < endcol; ++j)
			if(scoresmatrix[i][j] >= 0)
			{
				nancli[i-begin] = false;
				nanimp[j-begcol] = false;
			}

		if(nancli[i-begin])
		{
			nanclifile << i-begin+1 << endl;
			validclients--;
		}
	}

	nanclifile.close();

	for(unsigned int i = 0; i < cols; ++i)
		if(nanimp[i])
		{
			validimpostors--;
			nanimpfile << (i/takes)+1 << " " << (i%takes)+2 << endl;
		}

	nanimpfile.close();

	// Calculate the client score and FRR histograms
	for(unsigned int i = begin; i < end; i++)
	{
		if(!nancli[i-begin])
		{
			for(unsigned int j = begcol; j < takes*i; j++)
			{
				if(scoresmatrix[i][j] >= 0)
				{
					pos = floor((scoresmatrix[i][j]-scoremin)/interval_size);
					dist_impostors[pos]++;

					for(unsigned int k = 0; k <= pos; k++)
						dist_far[k]++;
				}
			}
			for(unsigned int j = takes*i; j < takes*(i+1); j++)
			{
				if(scoresmatrix[i][j] >= 0)
				{
					pos = floor((scoresmatrix[i][j]-scoremin)/interval_size);
					dist_clients[pos+1]++;

					for(unsigned int k = pos+1; k <= SCORE_INTERVALS; k++)
						dist_frr[k]++;
				}
			}
			for(unsigned int j = takes*(i+1); j < endcol; j++)
			{
				if(scoresmatrix[i][j] >= 0)
				{
					pos = floor((scoresmatrix[i][j]-scoremin)/interval_size);
					dist_impostors[pos]++;

					for(unsigned int k = 0; k <= pos; k++)
						dist_far[k]++;
				}
			}
		}
	}

	// Get the relative frequences and the verification accuracy measures
	for(unsigned int j = 0; j <= SCORE_INTERVALS; j++)
	{
		ddist_clients[j]   = dist_clients[j]   * 100.0/(validclients*takes);
		ddist_impostors[j] = dist_impostors[j] * 100.0/(validimpostors*(validclients-1));
		ddist_frr[j]       = dist_frr[j]       * 100.0/(validclients*takes);
		ddist_far[j]       = dist_far[j]       * 100.0/(validimpostors*(validclients-1));

		distcliimp << j*interval_size + scoremin << " " << ddist_clients[j] << " " << ddist_impostors[j] << endl;
		frrfar << j*interval_size + scoremin << " " << ddist_frr[j] << " " << ddist_far[j] << endl;

		// Calculate the EER
		if(eer == -1 && ddist_frr[j] >= ddist_far[j])
		{
			if(j == 0)
				aux = 0;
			else
				aux = j-1;

			if(ddist_frr[j] - ddist_far[j] < ddist_far[aux] - ddist_frr[aux])
				eer = (ddist_frr[j] + ddist_far[j])/2.0;
			else
				eer = (ddist_frr[aux] + ddist_far[aux])/2.0;
		}

		// Calculate the FAR100
		if(far100 == -1 && ddist_far[j] < 1.0)
			far100 = ddist_frr[j];

		// Calculate the FAR1000
		if(far1000 == -1 && ddist_far[j] < 0.1)
			far1000 = ddist_frr[j];

		// Calculate the FAR10000
		if(far10000 == -1 && ddist_far[j] < 0.01)
			far10000 = ddist_frr[j];

		// Calculate the ZeroFAR
		if(zerofar == -1 && ddist_far[j] == 0.0)
			zerofar = ddist_frr[j];

		// Calculate the ZeroFRR
		if(zerofrr == -1 && ddist_frr[j] > 0.0)
		{
			if(j == 0)
				zerofrr = ddist_far[0];
			else
				zerofrr = ddist_far[j-1];
		}
	}

	if(far100 == -1) far100 = 100;
	if(far1000 == -1) far1000 = 100;
	if(far10000 == -1) far10000 = 100;
	if(zerofar == -1) zerofar = 100;
	if(zerofrr == -1) zerofrr = 100;

	distcliimp.close();
	frrfar.close();

	// Calculate the CMC
	// This loop is executed once per each input fingerprint
	for(unsigned int i = begin*takes; i < begin*takes+cols; i++)
	{
		if(!nanimp[i-begin*takes])
		{
			// Calculate the position of the corresponding template fingerprint
			pos = i/takes;
			
			if(!nancli[pos-begin])
			{
				unsigned int bigger = 0;

				// Count all template fingerprints whose score with the input fingerprint is higher or equal to the genuine.
				for(unsigned int j = begin; j < end; j++)
					if(scoresmatrix[j][i] >= scoresmatrix[pos][i])
						bigger++;

				// Increase the number of successful identifications, when the rank is higher or equal to "bigger"
				for(unsigned int j = bigger-1; j < RANKS; j++)
					cmc[j]++;
			}
		}
	}
// cout << "validimpostors = " << validimpostors << endl;
	for(unsigned int j = 0; j < RANKS; j++)
	{
		// Calculate val as the fraction of successful identifications ("cols" is the total)
		float val = (float)cmc[j] / (validimpostors - (rows-validclients)*takes);

		if(val >= 0.99 && r100 == RANKS)
			r100 = j+1;
		if(val >= 0.999 && r1000 == RANKS)
			r1000 = j+1;
		if(val >= 0.9999 && r10000 == RANKS)
			r10000 = j+1;
		if(val == 1.0 && r1 == RANKS)
			r1 = j+1;

		cmcfile << j+1 << " " << val << endl;
	}

	cmcfile.close();
	
	verif << "EER\tFAR100\tFAR1000\tFAR10000\tZero_FAR\tZero_FRR\tR100\tR1000\tR10000\tZeroR" << endl;

	verif << eer << "\t" << far100 << "\t" << far1000 << "\t" << far10000 << "\t" << zerofar << "\t" << zerofrr <<  "\t" << r100 << "\t" << r1000 << "\t" << r10000 << "\t" << r1 << endl;
	
	verif.close();
}

