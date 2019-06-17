/**
 * \file    Score.h
 * \author  Daniel Peralta <dperalta@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file for the Score class.
 */

#ifndef SCORE_H
#define SCORE_H

#ifdef MPI_VERSION
#include "mpi.h"
#endif

#include <string>
#include <cstring>

class Score
{
	
protected:

	static const int IDLENGTH = 200;
	float score;
	int relindex;
	int absindex;
	int process;
	char id[IDLENGTH];
	
public:
	
	Score();
	
	Score(const Score &s);
	
	Score(float s, int reli, int absi, int proc, const std::string &id);
	
// 	virtual ~Score() {}
	
	float getScore() const;
	
	int getRelIndex() const;
	
	int getAbsIndex() const;
	
	int getProcess() const;
	
	std::string getId() const;
	
	void setScore(float s);
	
	void setRelIndex(int s);
	
	void setAbsIndex(int s);
	
	void setProcess(int s);
	
	void setId(const std::string &s);

	Score & operator=(const Score &s);

	#ifdef MPI_VERSION
	static MPI::Datatype getDatatype();
	#endif

	static bool isnull(const Score &s);

	bool operator<(const Score &s) const;
	bool operator>(const Score &s) const;
	static bool better(const Score &s1, const Score &s2);
	
};

inline float Score::getScore() const {return score;}
inline int Score::getRelIndex() const {return relindex;}
inline int Score::getAbsIndex() const {return absindex;}
inline int Score::getProcess() const {return process;}
inline std::string Score::getId() const {return std::string(id);}
inline void Score::setScore(float s) {score = s;}
inline void Score::setRelIndex(int s) {relindex = s;}
inline void Score::setAbsIndex(int s) {absindex = s;}
inline void Score::setProcess(int s) {process = s;}
inline void Score::setId(const std::string &s) {strcpy(id, s.c_str());}
inline bool Score::isnull(const Score &s) {return s.score == -1;}
inline bool Score::operator<(const Score &s) const {return score < s.score;}
inline bool Score::operator>(const Score &s) const {return score > s.score;}
inline bool Score::better(const Score &s1, const Score &s2) {return s1.score > s2.score;}


#endif
