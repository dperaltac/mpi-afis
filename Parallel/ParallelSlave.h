/**
 * \file    ParallelSlave.h
 * \author  Daniel Peralta <dperalta@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file for the ParallelSlave class.
 */

#ifndef PARALLEL_SLAVE_H
#define PARALLEL_SLAVE_H

#include <string>
#include <vector>
#include "omp.h"
#include "mpi.h"

#include "ParallelHandler.h"
#include "Fingerprint.h"

/**
 * \class ParallelSlave
 *
 * The ParallelSlave class handles the parallel structure of the fingerprint recognizer.
 * It encapsulates the search through a partition of the database.
 */
class ParallelSlave : public ParallelHandler {

public:

	/** Default constructor */
	ParallelSlave();

	/** Default destructor */
	virtual ~ParallelSlave();

	/**
	 * The process sends its results to the master, according with the selected stop mode
	*/
	void sendResultsToMaster();


	/**
	 * All processes send the found matches to the \a ParallelMaster
	 */
	void sendMatches() const;

	/**
	 * Each process sends the found match to the \a ParallelMaster
	 */
	void sendSingleMatches() const;


	/**
	 * Runs the matching algorithm for each fingerprint of the assigned database part, and sends the results to the \a ParallelMaster
	 * \param newF Fingerprint to be found in the database
	 * \param v Vector with the matches found by the current process
	 */
	virtual void run(Fingerprint *newF, std::vector< Score > &v);


	/**
	 * Runs the matching algorithm for each fingerprint of the assigned database part, and sends the results to the \a ParallelMaster
	 * \param vinput Vector with the input fingerprints for a user
	 * \param v Vector with the matches found by the current process
	 */
	virtual void run(std::vector<Fingerprint *> &vinput, std::vector< Score > &v);


	/**
	 * Runs the matching algorithm for each fingerprint of the assigned database part
	 * \param newF Fingerprint to be found in the database
	 * \param v Vector with the matches found by the current process
	 */
	virtual void runFirst(Fingerprint *newF, std::vector< Score > &v);


	/**
	 * Runs the matching algorithm for each fingerprint of the assigned database part
	 * \param vinput Vector with the input fingerprints for a user
	 * \param v Vector with the matches found by the current process
	 */
	virtual void runFirst(std::vector<Fingerprint *> &vinput, std::vector< Score > &v);


	/**
	 * Explores the database looking for the fingerprint, using only the fingerprints signaled in \p candidates
	 * \param newF Fingerprint to be found in the database
	 * \param candidates Vector of candidates
	 * \param v Vector with the matches found by the current process
	 */
	virtual void runSecond(Fingerprint *newF, const std::vector< Score > &candidates, std::vector< Score > &v);


	/**
	 * Explores the database looking for the fingerprint, using only the fingerprints signaled in \p candidates
	 * \param vinput Vector with the input fingerprints for a user
	 * \param candidates Vector of candidates
	 * \param v Vector with the matches found by the current process
	 */
	virtual void runSecond(std::vector<Fingerprint *> &vinput, const std::vector< Score > &candidates, std::vector< Score > &v);


	/**
	 * Reads the database part to be processed by the process.
	 * \param filename Name of the database file
	 * \param quality Minimum quality for the read minutiae. It is an integer in the range [0,100]
	 */
	void loadDBFile(const std::string &filename, unsigned int quality);



	/**
	 * Reads the database part to be processed by the process
	 * \param filenames Vector with the names of the database files for each fingerprint
	 * \param quality Minimum quality for the read minutiae. It is an integer in the range [0,100]
	 */
	void loadDBFile(const std::vector<std::string> &filenames, unsigned int quality);
	
	
	/**
	 * Not supported.
	 */
	virtual void loadDBFile(const std::string &filename, const std::string &classfile, unsigned int quality);


	/**
	 * Reads and initializes the database part to be processed by the process.
	 * \param filename Name of the database file
	 * \param quality Minimum quality for the read minutiae. It is an integer in the range [0,100]
	 */
	void loadInitializeDBFile(const std::string &filename, unsigned int quality);



	/**
	 * Reads and initializes the database part to be processed by the process
	 * \param filenames Vector with the names of the database files for each fingerprint
	 * \param quality Minimum quality for the read minutiae. It is an integer in the range [0,100]
	 */
	void loadInitializeDBFile(const std::vector<std::string> &filenames, unsigned int quality);
	
	
	/**
	 * Not supported.
	 */
	virtual void loadInitializeDBFile(const std::string &filename, const std::string &classfile, unsigned int quality);


	/**
	 * Creates a structure to detect if one of the other processes sends a stop signal.
	 * The signal can be checked with the method \a checkStop()
	 */
	virtual void waitForStop();

	/**
	 * Initializes all the fingerprints in the database and synchronizes some values
	 */
	void initializeDatabase();

	unsigned int getBegin() const;
	unsigned int getEnd() const;

protected:

	std::vector< Score > found_matches; //!< Set of found matches
	unsigned int begin; //!< First fingerprint in the database for this slave
	unsigned int end; //!< Last fingerprint in the database for this slave
	Matrix<Fingerprint *> fp_list; //!< Template fingerprint database
	bool found_other; //!< Whether another slave has found a match
	bool found; //!< Whether this slave has found a match


	/**
	* Sets the variables to read the database parts.
	* \param num_files Number of fingerprints
	* \param num_fingers Number of fingers
	 */
	void setListLimits(unsigned int num_files, unsigned int num_fingers);

	/**
	 * Sends a stop signal to every process
	 */
	void sendStopSignal(const Score &match) const;

	/**
	 * Sends a not found signal to every process
	 */
	void sendNotFoundSignal() const;

	/**
	 * Checks if one or several processes have sent a stop signal since the last call to \a waitForStop()
	 * The index of these processes is stored within the object.
	 * \return True if a stop signal has been received, false otherwise
	 */
	bool checkStop();

private:

	/**
	 * Copy constructor
	 * \param other Object to copy from
	 */
	ParallelSlave(const ParallelSlave& other) {};

	/**
	 * Assignment operator
	 * \param other Object to assign from
	 * \return A reference to this object
	 */
	virtual ParallelSlave& operator=(const ParallelSlave& other) {return *this;};

};

inline unsigned int ParallelSlave::getBegin() const {return begin;}
inline unsigned int ParallelSlave::getEnd() const {return end;}

inline void ParallelSlave::loadDBFile(const std::string &filename, unsigned int quality)
{
	loadDBFile(std::vector<std::string>(1, filename), quality);
}


inline void ParallelSlave::loadDBFile(const std::string &filename, const std::string &classfile, unsigned int quality)
{
	exit("ERROR: ParallelMaster does not support classification", -1);
}

inline void ParallelSlave::loadInitializeDBFile(const std::string &filename, unsigned int quality)
{
	loadInitializeDBFile(std::vector<std::string>(1, filename), quality);
}


inline void ParallelSlave::loadInitializeDBFile(const std::string &filename, const std::string &classfile, unsigned int quality)
{
	loadDBFile(filename, classfile, quality);
}

inline void ParallelSlave::run(Fingerprint *newF, std::vector< Score > &v)
{
	std::vector<Fingerprint *> aux(1,newF);
	run(aux, v);
}

inline void ParallelSlave::runFirst(Fingerprint *newF, std::vector< Score > &v)
{
	std::vector<Fingerprint *> aux(1,newF);
	runFirst(aux, v);
}

inline void ParallelSlave::runSecond(Fingerprint *newF, const std::vector< Score > &candidates, std::vector< Score > &v)
{
	std::vector<Fingerprint *> aux(1,newF);
	runSecond(aux, candidates, v);
}

inline bool ParallelSlave::checkStop()
{
	return (stop_requests.front() == MPI::REQUEST_NULL || stop_requests.front().Test());
}

#endif // PARALLEL_SLAVE_H
