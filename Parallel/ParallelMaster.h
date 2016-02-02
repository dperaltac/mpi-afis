/**
 * \file    ParallelMaster.h
 * \author  Daniel Peralta <dperalta@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file for the ParallelMaster class.
 */

#ifndef PARALLEL_MASTER_H
#define PARALLEL_MASTER_H

#include <string>
#include <vector>
#include "omp.h"
#include "mpi.h"

#include "ParallelHandler.h"

/**
 * \class ParallelMaster
 *
 * The ParallelMaster class handles the parallel structure of the fingerprint recognizer.
 * It encapsulates the communication and the creation of the threads and processes.
 */
class ParallelMaster : public ParallelHandler
{

public:

	/** Default constructor */
	ParallelMaster();

	/** Default destructor */
	virtual ~ParallelMaster();


	/**
	 * Divides the interval [0,size) in one part for each process, in such a way that the last process may get a smaller part.
	 * The master thread processes the operation, then sends its interval to each process.
	 * \param size Size of the interval that must be split
	 */
	void setBoundaries(unsigned int size) const;


	/**
	 * The master process sends a message to each of the processes that have sent a stop signal.
	 * Then, it receives the values of these matches, and returns them.
	 * \return Vector with the matches found in all processes
	 */
	void retrieveMatchesStop(std::vector< Score > &v);


	/**
	 * Receives and returns all the matches calculated by the slave processes
	 * \return Vector with the matches found in all slave processes
	 */
	void retrieveAllMatches(std::vector< Score > &v) const;


	/**
	 * Receives and returns the best match among all the slave processes
	 * \return Vector with the single best match
	 */
	void retrieveBestMatch(std::vector< Score > &v) const;

	/**
	 * Reads the database part to be processed by the process
	 * \param filelists Vector with the names of the database files for each fingerprint
	 * \param quality Minium quality for the read minutiae. It is an integer in the range [1,100]
	 */
	virtual void loadDBFile(const std::vector<std::string> &filelists, unsigned int quality);

	/**
	 * Not supported.
	 */
	virtual void loadDBFile(const std::string &filename, const std::string &classfile, unsigned int quality);

	/**
	 * Reads and initializes the database part to be processed by the process
	 * \param filelists Vector with the names of the database files for each fingerprint
	 * \param quality Minium quality for the read minutiae. It is an integer in the range [1,100]
	 */
	virtual void loadInitializeDBFile(const std::vector<std::string> &filelists, unsigned int quality);

	/**
	 * Not supported.
	 */
	virtual void loadInitializeDBFile(const std::string &filename, const std::string &classfile, unsigned int quality);


	/**
	 * Main loop of the algorithm. Returns the matches found by the slave processes
	 * \param v Returned vector with the matches found by the slave processes
	 * \param numfps Number of fingerprints in the search database. Used to calculate statistics.
	 */
	virtual void run(std::vector< Score > &v, unsigned int numfps);


	/**
	 * Main loop of the algorithm. Returns the matches found by the slave processes
	 * \param newF Fingerprint to be found in the database (not used)
	 * \param v Returned vector with the matches found by the slave processes
	 */
	virtual void run(Fingerprint *newF, std::vector< Score > &v);


	/**
	 * Main loop of the algorithm. Returns the matches found by the slave processes
	 * \param vinput Vector with the input fingerprints for a user (not used)
	 * \param v Returned vector with the matches found by the slave processes
	 */
	virtual void run(std::vector<Fingerprint *> &vinput, std::vector< Score > &v);


	/**
	 * Same as \a run
	 * \param newF Fingerprint to be found in the database (not used)
	 * \param candidates Vector of candidates (not used)
	 * \param v Vector of matches
	 */
	virtual void runSecond(Fingerprint *newF, const std::vector< Score > &candidates, std::vector< Score > &v);


	/**
	 * Same as \a run
	 * \param vinput Vector with the input fingerprints for a user (not used)
	 * \param candidates Vector of candidates (not used)
	 * \param v Vector of matches
	 */
	virtual void runSecond(std::vector<Fingerprint *> &vinput, const std::vector< Score > &candidates, std::vector< Score > &v);
	

	void waitForStop();

protected:

	bool *found_others;
	bool found;

private:

	/**
	 * Copy constructor
	 * \param other Object to copy from
	 */
	ParallelMaster(const ParallelMaster& other) {};

	/**
	 * Assignment operator
	 * \param other Object to assign from
	 * \return A reference to this object
	 */
	virtual ParallelMaster& operator=(const ParallelMaster& other) {return *this;};

};

inline void ParallelMaster::runSecond(Fingerprint *newF, const std::vector< Score > &candidates, std::vector< Score > &v) { run(v, candidates.size()); }

inline void ParallelMaster::runSecond(std::vector<Fingerprint *> &vinput, const std::vector< Score > &candidates, std::vector< Score > &v) { run(v, candidates.size()); }

inline void ParallelMaster::run(Fingerprint *newF, std::vector< Score > &v) { run(v, totalsize); }

inline void ParallelMaster::run(std::vector<Fingerprint *> &vinput, std::vector< Score > &v) { run(v, totalsize); }

inline void ParallelMaster::loadDBFile(const std::string &filename, const std::string &classfile, unsigned int quality)
{
	exit("ERROR: ParallelMaster does not support classification", -1);
}

inline void ParallelMaster::loadInitializeDBFile(const std::string &filename, const std::string &classfile, unsigned int quality)
{
	exit("ERROR: ParallelMaster does not support classification", -1);
}

inline void ParallelMaster::loadInitializeDBFile(const std::vector<std::string> &filenames, unsigned int quality)
{
	loadDBFile(filenames, quality);
}

#endif // PARALLEL_MASTER_H
