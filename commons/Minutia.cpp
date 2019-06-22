/**
 * \file    Minutia.cpp
 * \author  Joaquin Derrac <jderrac@decsai.ugr.es>
 * \author  Salvador Garcia <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * The implementation of the minutia class.
 */

#include <iostream>
#include "Minutia.h"

Minutia::Minutia() : index(-1), X(0), Y(0), T(0), quality(-1), type(OTH)
{
}

Minutia::~Minutia()
{
}

Minutia::Minutia(int newIndex, int valX, int valY, float valT, unsigned int valQ, typeMin valM) :
	index(newIndex),
	X(valX),
	Y(valY),
	T(valT),
	quality(valQ),
	type(valM)
{
}

Minutia::Minutia(const Minutia& other) :
	index(other.index),
	X(other.X),
	Y(other.Y),
	T(other.T),
	quality(other.quality),
	type(other.type)
{
}

Minutia& Minutia::operator=(const Minutia& rhs)
{
	if (this == &rhs) return *this;

	X=rhs.X;
	Y=rhs.Y;
	T=rhs.T;
	quality=rhs.quality;
	index=rhs.index;
	type=rhs.type;

	return *this;
}

/**
 * Prints a minutia
 * \param output Output stream
 * \param M Minutia to print
 */
std::ostream& operator<<(std::ostream& output, const Minutia& M)
{
	output << M.index << ": (" <<  M.X << ", " << M.Y << ", " << M.T << ", " << M.quality << ", " <<
	(M.type==RIG?"RIG)":"BIF)");
	return output;
}
