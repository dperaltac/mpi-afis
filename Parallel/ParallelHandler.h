/**
 * \file    ParallelHandler.h
 * \author  Daniel Peralta <dperalta@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file for the ParallelHandler class.
 */

#ifndef PARALLEL_HANDLER_H
#define PARALLEL_HANDLER_H

#include "Fingerprint.h"
#include <string>
#include <vector>
#include "omp.h"
#include "mpi.h"
#include "Score.h"

/**
 * \class ParallelHandler
 *
 * The ParallelHandler class handles the parallel structure of the fingerprint recognizer.
 * It encapsulates some of the common processing for both master and slave processes processes.
 */
class ParallelHandler
{

public:

	typedef enum {STOP_YES = 1, STOP_NO, STOP_MAX, STOP_ALL, STOP_RANKING} stop_t;
	typedef enum {FP_UNKNOWN, FP_JIANG, FP_MCC} fingerprint_t;
	typedef enum {FUSION_SUM,FUSION_MAX,FUSION_MIN,FUSION_PROD} fusion_t;
	typedef enum {DPDDFF_SS,DPDDFF_SD,DPDDFF_DS,DPDDFF_DD} variant_t;

	static const unsigned int MASTER_ID = 0;

	/** Default constructor */
	ParallelHandler();

	/** Default destructor */
	virtual ~ParallelHandler();

	/**
	 * Changes the threshold for the search in the database
	 * \param t New threshold
	 */
	void setThreshold(float t);

	/**
	 * Changes the ranking for the search in the database
	 * \param r New ranking
	 */
	void setRanking(unsigned int r);

	/**
	 * Changes the used matching algorithm
	 * \param fp New fingerprint type
	 * \param argc Number of arguments for the fingeprint matcher
	 * \param argv Arguments for the fingeprint matcher
	 * \param finger Index of the finger
	 */
	void setFingerprintType(fingerprint_t fp, int argc, char *argv[], int finger = 0);

	/**
	 * Changes the used matching algorithm
	 * \param fp Vector of new fingerprint type (one type per finger)
	 * \param argc Number of arguments for the fingeprint matcher
	 * \param argv Arguments for the fingeprint matcher
	 */
	void setFingerprintType(const std::vector<fingerprint_t> &fp, int argc, char *argv[]);

	/**
	 * Returns the used matching algorithm
	 * \param finger Index of the finger
	 * \return Fingerprint type
	 */
	fingerprint_t getFingerprintType(int finger = 0) const;

	/**
	 * Returns the used matching algorithms
	 * \return Vector of fingerprint types
	 */
	std::vector<fingerprint_t> getFingerprintTypes() const;

	/**
	 * Changes the stop mode for the search in the database
	 * \param s New stop mode
	 */
	void setStopMode(stop_t s);

	/**
	 * Get the number of processes
	 * \return Number of MPI processes
	 */
	static int getProcesses();

	/**
	 * Get the number of threads
	 * \return Number of OpenMP threads
	 */
	int getThreads() const;

	/**
	 * Get the number of performed matches in the search, i.e. the number of fingeprints in the database that have been explored
	 * \return Number of explored fingerprints
	 */
	int getExplored() const;

	/**
	 * Get the size of the fingerprint database
	 * \return Size of the fingerprint database
	 */
	int getDBSize() const;

	/**
	 * Get the penetration rate, as the number of explored fingeprints (\a getExplored()) divided by the size of the database (\a getDBSize())
	 * \return Penetration rate
	 */
	virtual float getPenetrationRate() const;

	/**
	 * Get the number of the current process
	 * \return MPI ID of the current process. This ID can go from 0 to \a getProcesses().
	 */
	static int getProcessID();

	/**
	 * Get the iteration number (i.e. the number of fingerprints that have been searched in the database since the program started)
	 * \return Iteration number
	 */
	int getIteration() const;

	/**
	 * Get the average matching time per fingerprint
	 * \return Average matching time in seconds
	 */
	double getAvgMatchingTime() const;

	/**
	 * Resets the iteration number to zero (i.e. the number of fingerprints that have been searched in the database since the program started)
	 */
	void resetIteration();

	/**
	 * Determines if the current process is the master process
	 * \return true if the current process is the master processMPI ID of the current process, false otherwise
	 */
	bool isMaster() const;

	/**
	 * Initialize MPI and OpenMP
	 * \param threads Number of OpenMP threads
	 * \param classification Indicates whether classification is performed
	 * \param level Required parallelism level for MPI
	 * \return A pointer to the subclass of ParallelHandler
	 * \post This pointer is allocated dynamically with new, so it must be deallocated with delete
	 */
	static ParallelHandler *getHandler(int threads, bool classification = false, int level = MPI_THREAD_FUNNELED);

	/**
	 * Exits the program, showing an error message and calling MPI::Finalize()
	 * \param message Error message to be shown in the standard error output.
	 * \param status Number that is returned to the SO
	 */
	static void exit(const std::string & message, int status = 0);


	/**
	 * Main loop of the process, that explores the database looking for the fingerprint.
	 * The behavior differs depending on the process type (master or slave)
	 * \param newF Fingerprint to be found in the database
	 * \param v Vector of matches
	 * \returns The vector of found matches in \p v
	 */
	virtual void run(Fingerprint *newF, std::vector< Score > &v) = 0;


	/**
	 * Main loop of the process, that explores the database looking for the fingerprint.
	 * The behavior differs depending on the process type (master or slave)
	 * \param vinput Vector with the input fingerprints for a user
	 * \param v Vector of matches
	 * \returns The vector of found matches in \p v
	 */
	virtual void run(std::vector<Fingerprint *> &vinput, std::vector< Score > &v) = 0;


	/**
	 * Explores the database looking for the fingerprint, using only the fingerprints signaled in \p candidates
	 * The behavior differs depending on the process type (master or slave)
	 * \param newF Fingerprint to be found in the database
	 * \param candidates Vector of candidates
	 * \param v Vector of matches
	 * \returns The vector of found matches in \p v
	 */
	virtual void runSecond(Fingerprint *newF, const std::vector< Score > &candidates, std::vector< Score > &v) = 0;


	/**
	 * Explores the database looking for the fingerprint, using only the fingerprints signaled in \p candidates
	 * The behavior differs depending on the process type (master or slave)
	 * \param vinput Vector with the input fingerprints for a user
	 * \param candidates Vector of candidates
	 * \param v Vector of matches
	 * \returns The vector of found matches in \p v
	 */
	virtual void runSecond(std::vector<Fingerprint *> &vinput, const std::vector< Score > &candidates, std::vector< Score > &v) = 0;


	/**
	 * Creates a fingerprint of the appropriate type, according to the information set in \a setFingerprintType
	 * \param finger Index of the finger
	 * \return Pointer to a subclass of \a Fingerprint
	 * \post The returned pointer is allocated with new, so it must be deallocated with delete
	 */
	virtual Fingerprint *newFingerprint(int finger = 0) const;

	/**
	 * Aggregates the scores according to the type of fusion in \a fusion
	 * \param fscore Vector of scores
	 * \return Fused score
	 */
	float aggregateScores(const std::vector<float> &fscore) const;


	/**
	 * Every slave process sends its database size to the master, that sums them up.
	 */
	void calculateDBSize();


	/**
	 * Reads the database part to be processed by the process
	 * \param filename Name of the database file
	 * \param quality Minium quality for the read minutiae. It is an integer in the range [1,100]
	 */
	virtual void loadDBFile(const std::string &filename, unsigned int quality);

	/**
	 * Reads the database part to be processed by the process
	 * \param filelists Vector with the names of the database files for each fingerprint
	 * \param quality Minium quality for the read minutiae. It is an integer in the range [1,100]
	 */
	virtual void loadDBFile(const std::vector<std::string> &filelists, unsigned int quality) = 0;

	/**
	 * Reads the database part to be processed by the process.
	 * \param filename Name of the database file
	 * \param classfile Name of the file that contains the class for each fingerprint in \p filename. If empty, no classification is performed.
	 * \param quality Minimum quality for the read minutiae. It is an integer in the range [0,100]
	 */
	virtual void loadDBFile(const std::string &filename, const std::string &classfile, unsigned int quality) = 0;

	/**
	 * Reads the database part to be processed by the process.
	 * \param filelists Vector with the names of the database files for each fingerprint
	 * \param classfiles Name of the file that contains the class for each fingerprint in \p filename. If empty, no classification is performed.
	 * \param quality Minimum quality for the read minutiae. It is an integer in the range [0,100]
	 */
	virtual void loadDBFile(const std::vector<std::string> &filelists, const std::vector<std::string> &classfiles, unsigned int quality);


	/**
	 * Reads and initializes the database part to be processed by the process
	 * \param filename Name of the database file
	 * \param quality Minium quality for the read minutiae. It is an integer in the range [1,100]
	 */
	virtual void loadInitializeDBFile(const std::string &filename, unsigned int quality);

	/**
	 * Reads and initializes the database part to be processed by the process
	 * \param filelists Vector with the names of the database files for each fingerprint
	 * \param quality Minium quality for the read minutiae. It is an integer in the range [1,100]
	 */
	virtual void loadInitializeDBFile(const std::vector<std::string> &filelists, unsigned int quality) = 0;

	/**
	 * Reads and initializes the database part to be processed by the process.
	 * \param filename Name of the database file
	 * \param classfile Name of the file that contains the class for each fingerprint in \p filename. If empty, no classification is performed.
	 * \param quality Minimum quality for the read minutiae. It is an integer in the range [0,100]
	 */
	virtual void loadInitializeDBFile(const std::string &filename, const std::string &classfile, unsigned int quality) = 0;

	/**
	 * Reads and initializes the database part to be processed by the process.
	 * \param filelists Vector with the names of the database files for each fingerprint
	 * \param classfiles Name of the file that contains the class for each fingerprint in \p filename. If empty, no classification is performed.
	 * \param quality Minimum quality for the read minutiae. It is an integer in the range [0,100]
	 */
	virtual void loadInitializeDBFile(const std::vector<std::string> &filelists, const std::vector<std::string> &classfiles, unsigned int quality);

	/**
	* Set the fusion type for the fingerprint database
	* \param fus New fusion type
	*/
	void setFusion(fusion_t fus);

	/**
	* Returns the number of fingers or the number of matchers. This is the dimension of the
	* input vector for \a run.
	* \return Number of fingers or matchers
	*/
	virtual unsigned int getNumFingers() const;

	/**
	 * Every slave process sends the number of explored fingerprints to the master, which sums them up.
	 */
	void calculateExplored();

protected:

	static const int TAG_FOUND = 0;
	static const int TAG_BOUNDARIES = 1;
	static const int TAG_MATCH = 2;
	static const int TAG_OFFSET = 3;

	std::vector<MPI::Request> stop_requests; //!< Set of stop requests from each slave. Used with when the seach is stopped when a score threshold is reached.
	stop_t stop; //!< Stopping criterion for the search.
	float threshold; //!< Threshold to stop the search or select candidates.
	unsigned int ranking; //!< Ranking to stop the search or select candidates.
	int iteration; //!< Number of input fingerprints searched.
	int explored; //!< Number of matchings performed in the last seach.
	int totalsize; //!< Total number of fingerprints (in this node or in all nodes)
	double avgmatchingtime; //!< Average matching time, in seconds
	unsigned int numfingers; //!< Number of fingers per identity
	fusion_t fusion; //!< Fusion type

	std::vector<fingerprint_t> template_fp; //!< Matcher used for each finger
	MPI::Datatype mpi_pair_t; //!< MPI Datatype for the scores

	/**
	 * Frees all the requests in the \a stop_requests vector.
	 */
	virtual void freeRequests();


private:

	/**
	 * Copy constructor
	 * \param other Object to copy from
	 */
	ParallelHandler(const ParallelHandler& other) {};

	/**
	 * Assignment operator
	 * \param other Object to assign from
	 * \return A reference to this object
	 */
	virtual ParallelHandler& operator=(const ParallelHandler& other) {return *this;};

};


inline unsigned int ParallelHandler::getNumFingers() const {return numfingers;}
inline void ParallelHandler::setThreshold(float t) {threshold = t;}
inline void ParallelHandler::setRanking(unsigned int r) {ranking = r;}
inline std::vector<ParallelHandler::fingerprint_t> ParallelHandler::getFingerprintTypes() const {return template_fp;}
inline ParallelHandler::fingerprint_t ParallelHandler::getFingerprintType(int finger) const {return template_fp[finger];}
inline void ParallelHandler::setStopMode(stop_t s) {stop = s;}
inline int ParallelHandler::getProcesses() {return MPI::COMM_WORLD.Get_size();}
inline int ParallelHandler::getThreads() const {return omp_get_num_threads();}
inline int ParallelHandler::getProcessID() {return MPI::COMM_WORLD.Get_rank();}
inline int ParallelHandler::getIteration() const {return iteration;}
inline void ParallelHandler::resetIteration() {iteration = 0;}
inline int ParallelHandler::getExplored() const {return explored;}
inline int ParallelHandler::getDBSize() const {return totalsize;}
inline double ParallelHandler::getAvgMatchingTime() const {return avgmatchingtime;}
inline bool ParallelHandler::isMaster() const {return (getProcessID() == MASTER_ID);}

inline void ParallelHandler::loadDBFile(const std::string &filename, unsigned int quality)
{
	loadDBFile(std::vector<std::string>(1, filename), quality);
}

inline void ParallelHandler::loadDBFile(const std::string &filename, const std::string &classfile, unsigned int quality)
{
	loadDBFile(std::vector<std::string>(1, filename), quality);
}

inline void ParallelHandler::loadInitializeDBFile(const std::string &filename, unsigned int quality)
{
	loadInitializeDBFile(std::vector<std::string>(1, filename), quality);
}

inline void ParallelHandler::loadInitializeDBFile(const std::string &filename, const std::string &classfile, unsigned int quality)
{
	loadInitializeDBFile(std::vector<std::string>(1, filename), quality);
}

inline float ParallelHandler::getPenetrationRate() const
{
	return ((float)explored)/totalsize;
}

inline void ParallelHandler::setFusion(ParallelHandler::fusion_t fus) {fusion = fus;}

#endif // PARALLEL_HANDLER_H
