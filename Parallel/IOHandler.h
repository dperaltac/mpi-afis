/**
 * \file    IOHandler.h
 * \author  Daniel Peralta <dperalta@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * This file defines class \c IOHandler
 */

#ifndef IOHANDLER_H
#define IOHANDLER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <iterator>

#include "ParallelHandler.h"

class IOHandler
{

public:
	// Data types
	typedef enum {OUT_SHORT, OUT_LONG, OUT_BYFILE} output_t;

	/**
	 * Parameter constructor
	 * @param argc Number of command line parameters
	 * @param argv Command line parameters
	 */
	IOHandler(int argc, char *argv[]);

	/**
	 * Destructor
	 */
	~IOHandler();

	/**
	* The master process reads a string from the file \a finger and broadcasts it.
	* The slave processes receive the string from the master.
	* If it reaches the end of file, it returns "exit".
	* \param finger Number of the input file
	* \param cinput Read string
	* \param parallel_handler Handler for the parallel environment
	*/
	void readFileName(int finger, std::string &cinput, const ParallelHandler *parallel_handler);

	/**
	* The master process reads a string from the first input file and broadcasts it.
	* The slave processes receive the string from the master.
	* If it reaches the end of file, it returns "exit".
	* \param cinput Read string
	* \param parallel_handler Handler for the parallel environment
	*/
	void readFileName(std::string &cinput, const ParallelHandler *parallel_handler);

	/**
	* The master process reads a string from each input file and broadcasts it.
	* The slave processes receive all the strings from the master.
	* If it reaches the end of file, it returns a vector of "exit".
	* \param cinput Read strings
	* \param parallel_handler Handler for the parallel environment
	*/
	void readFileName(std::vector<std::string> &cinput, const ParallelHandler *parallel_handler);

	/**
	* Reads one or more matrices with the features of the input fingerprints.
	* \returns Vector of matrices with the features of the input fingerprints. Each row corresponds to a fingerprint, in the same order as in the input file. Each column corresponds to a feature. Each matrix corresponds to a set of features.
	*/
	std::vector< Matrix<double> > readFeatureMatrices();

	/**
	* Read a string from the first input file.
	* \param cinput Read string
	*/
	void readFileNameLocal(std::string &cinput);

	/**
	* Read a string from the input file \a finger.
	* \param finger Number of the input file
	* \param cinput Read string
	*/
	void readFileNameLocal(int finger, std::string &cinput);



	/**
	* Returns the type of fingerprint coded in \a cinput
	* \param cinput String that contains a fingerprint type
	* \returns Type of fingerprint
	*/
	static ParallelHandler::fingerprint_t getFpType(const std::string &cinput);

	/**
	* Returns the type of fusion coded in \a cinput
	* \param cinput String that contains a fusion type
	* \returns Type of fusion
	*/
	static ParallelHandler::fusion_t getFusionType(const std::string &cinput);

	/**
	* Returns the type of stop criterion for the search coded in \a cinput,
	* and the associated variables
	* \param cinput String that contains a stop criterion
	* \returns Type of stop criterion and the associated variables.
	*/
	static ParallelHandler::stop_t getStopType(const std::string &cinput);

	/**
	* Returns the type of output coded in \a cinput
	* \param cinput String that contains an output type
	* \returns Output type
	*/
	static output_t getOutputType(const std::string &cinput);

	/**
	* Returns the DPD-DFF variant coded in \a cinput
	* \param cinput String that contains a variant
	* \returns DPD-DFF variant
	*/
	static ParallelHandler::variant_t getVariant(const std::string &cinput);

	/**
	* Returns the type of fingerprint
	* \param finger Number of the finger
	* \returns Type of fingerprint
	*/
	ParallelHandler::fingerprint_t getFpType(int finger = 0) const;

	/**
	* Returns the type of fusion
	* \returns Type of fusion
	*/
	ParallelHandler::fusion_t getFusionType() const;

	/**
	* Returns the types of the fingerprints
	* \returns Types of fingerprints
	*/
	std::vector<ParallelHandler::fingerprint_t> getFpTypes() const;

	/**
	* Returns the type of stop criterion for the search
	* \returns Type of stop criterion
	*/
	ParallelHandler::stop_t getStopType() const;

	/**
	* Returns the type of output
	* \returns Output type
	*/
	output_t getOutputType() const;

	/**
	* Returns the variant for DPD-DFF
	* \returns Variant (one of SS, SD, DS, DD)
	*/
	ParallelHandler::variant_t getVariant() const;

	/**
	 * Prints a message in stderr explaining the use of the program
	 */
	static void printSyntax();
	
	/**
	 * Prints the value of all the parameters.
	 */
	void printParameters() const;

	unsigned int getThreads() const;
	unsigned int getRanking() const;
	double getThreshold(int finger=0) const;
	unsigned int getQuality() const;
	bool getClassification() const;

	std::vector<std::string> getTemplateFiles() const;
	std::string getTemplateFile(int finger = 0) const;
	std::vector<std::string> getClassFiles() const;
	std::string getClassFile(int finger = 0) const;
	std::vector<std::string> getInputFiles() const;

	std::string getOutputPath() const;
	std::vector<std::string> getClassifierPaths() const;
	std::string getClassifierPath(int finger = 0) const;

	void printInitialOutput();
	void printIdentificationOutput(const std::string &cinput, const std::vector<Score> &matches, double avgmatchingtime, const std::vector<double> &times, int iteration, float num_candidates, float penetrationrate);
	void printFinalOutput(const ParallelHandler *parallel_handler);

	int getNumInputFingerprints(int finger = 0) const;

	template<class T>
	std::string printTabSeparated(const std::vector<T> &v, const char *sep="\t") const;

	void printScores(const std::vector<Score> &matches) const;
	
	void clearCache() const;

protected:

// 	vector<string> argv;

	std::vector<ParallelHandler::fingerprint_t> fp_type;
	ParallelHandler::fusion_t fusion;
	unsigned int threads;
	ParallelHandler::stop_t stop;
	unsigned int ranking;
	std::vector<double> threshold;
	unsigned int quality;
	output_t output;
	ParallelHandler::variant_t variant;

	std::vector<std::string> templatefiles;
	std::vector<std::string> inputfiles;
	std::vector<std::string> classfiles;
	std::vector<std::string> featurefiles;

	std::vector<std::ifstream *> inputfd;

	std::ofstream scorefile;
	std::ofstream timefile;
	std::string outpath;
	double meantime;
	
	bool inputclasses;
	std::vector<std::string> classifierpaths;
};


inline void IOHandler::readFileName(std::string &cinput, const ParallelHandler *parallel_handler)
{
	readFileName(0, cinput, parallel_handler);
}

inline void IOHandler::readFileNameLocal(std::string &cinput) { readFileNameLocal(0, cinput); }
inline ParallelHandler::fingerprint_t IOHandler::getFpType(int finger) const {return fp_type[finger];}
inline ParallelHandler::fusion_t IOHandler::getFusionType() const {return fusion;}
inline std::vector<ParallelHandler::fingerprint_t> IOHandler::getFpTypes() const {return fp_type;}
inline ParallelHandler::stop_t IOHandler::getStopType() const {return stop;}
inline IOHandler::output_t IOHandler::getOutputType() const {return output;}
inline unsigned int IOHandler::getThreads() const {return threads;}
inline unsigned int IOHandler::getRanking() const {return ranking;}
inline bool IOHandler::getClassification() const {return !classfiles.empty();}
inline double IOHandler::getThreshold(int finger) const {return threshold[finger];}
inline unsigned int IOHandler::getQuality() const {return quality;}
inline ParallelHandler::variant_t IOHandler::getVariant() const {return variant;}
inline std::vector<std::string> IOHandler::getTemplateFiles() const {return templatefiles;}
inline std::vector<std::string> IOHandler::getClassFiles() const {return classfiles;}
inline std::string IOHandler::getTemplateFile(int finger) const {return templatefiles[finger];}
inline std::string IOHandler::getClassFile(int finger) const {return classfiles[finger];}
inline std::vector<std::string> IOHandler::getInputFiles() const {return inputfiles;}
inline std::string IOHandler::getOutputPath() const {return outpath;}
inline std::vector<std::string> IOHandler::getClassifierPaths() const {return classifierpaths;}
inline std::string IOHandler::getClassifierPath(int finger) const {return classifierpaths[finger];}

template<class T>
std::string IOHandler::printTabSeparated(const std::vector<T> &v, const char *sep) const
{
	std::stringstream ss;

	typename std::vector<T>::const_iterator i = v.begin();
	ss << *i;

	for(i = i+1; i != v.end(); ++i)
		ss << sep << *i;

	return ss.str();
}

#endif
