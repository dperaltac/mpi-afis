/**
 * \file    Functions.h
 * \author  Joaquin Derrac <jderrac@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file defining some useful functions
 */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <utility>
#include <vector>
#include <string>

//********************************************
//    Rounding
//********************************************


struct sort_pred {
    bool operator()(const std::pair<float,std::pair<int,int> > &left, const std::pair<float,std::pair<int,int> > &right) {
        return left.first > right.first;
    }
};

/**
 * Rounds a float to its nearest integer. Supports both negative and positive values
 * \param r Value to round
 * \return r value rounded to the nearest integer
 */
int roundInt(float r);

/**
 * Computes the (nearest) inner angle between two vectors, in radians
 * \param x11 X coordinate of the first point of the first vector
 * \param y11 Y coordinate of the first point of the first vector
 * \param x12 X coordinate of the second point of the first vector
 * \param y12 Y coordinate of the second point of the first vector
 * \param x21 X coordinate of the first point of the second vector
 * \param y21 Y coordinate of the first point of the second vector
 * \param x22 X coordinate of the second point of the second vector
 * \param y22 Y coordinate of the second point of the second vector
 * \return Inner angle between the vectors, in radians [0,PI]
 */
float computeInnerAngle(int x11,int y11,int x12,int y12,int x21,int y21,int x22,int y22);

/**
 * Computes the angle between a vector and the X axis, clockwise, starting at 9:00, in degrees.
 * \param x1 X coordinate of the first point
 * \param x2 X coordinate of the second point
 * \param y1 Y coordinate of the first point
 * \param y2 Y coordinate of the second point
 * \return Angle between the vector and the X axis [0,360)
 */
float computeAngle(int x1,int x2,int y1,int y2);

/**
 * Computes the angle between a vector and the X axis, clockwise, starting at 9:00, in radians.
 * \param x1 X coordinate of the first point
 * \param x2 X coordinate of the second point
 * \param y1 Y coordinate of the first point
 * \param y2 Y coordinate of the second point
 * \return Angle between the vector and the X axis [0,2*PI)
 */
float computeAngleRad(int x1,int x2,int y1,int y2);

/**
 * Computes the distance from a point to a line
 * \param cx X coordinate of the point
 * \param cy XY coordinate of the point
 * \param ax X coordinate of the first point of the segment
 * \param ay Y coordinate of the first point of the segment
 * \param bx X coordinate of the second point of the segment
 * \param by Y coordinate of the second point of the segment
 * \return distance from the point to the infinite line and segment
 */

void DistanceFromLine(double cx, double cy, double ax, double ay ,
                                          double bx, double by, double &distanceSegment,
                                          double &distanceLine);

/**
 * Computes the distance from a point to a segment
 * \param cx X coordinate of the point
 * \param cy XY coordinate of the point
 * \param ax X coordinate of the first point of the segment
 * \param ay Y coordinate of the first point of the segment
 * \param bx X coordinate of the second point of the segment
 * \param by Y coordinate of the second point of the segment
 * \return distance from the point to the segment
 */

double DistanceFromSegment(double cx, double cy, double ax, double ay ,
                                          double bx, double by);

/**
 * Computes the left probability under the gaussian curve (0, 1)
 * \param z value

 * \return The probability
 */
inline float doLeft(double z)
{
	return (1.0+erf(z*0.707106781))*0.5;
}


/**
 * Computes the square power of the given number
 * \param z Number to be multiplied by itself

 * \return \p z to the square power
 */
template<class T>
inline T square(T z) {return z*z;}


/**
 * Turns a big-endian 2-byte number into a little-endian one, or viceversa
 * \param s Number to be converted
 *
 * \return Opposite endian version of \p s
 */
short shortIntegerSwap( short s );


/**
 * Turns a big-endian 2-byte number into a little-endian one, or viceversa
 * \param s Number to be converted
 *
 * \return Opposite endian version of \p s
 */
int longIntegerSwap (int i);


/**
* Function used in MatcherDeng
*
**/
int dtris2 ( int point_num, int base, double point_xy[], int *tri_num, int tri_vert[], int tri_nabe[] );


/**
 * Gets a 2-bytes integer number from its binary representation
 * \param a Most significant byte
 * \param b Less significant byte
 *
 * \return 2-bytes integer equivalent to 256*a + b
 */
inline unsigned short int getShortInt (unsigned char a, unsigned char b)
{
	return 256*static_cast<unsigned char>(a) + static_cast<unsigned char>(b);
}


template<typename T, typename P>
T removeIf(T beg, T end, P pred)
{
	T dest = beg;
	for (T itr = beg;itr != end; ++itr)
		if (!pred(*itr))
			*(dest++) = *itr;
	return dest;
}



/** Computes the difference between two angles in radians in the interval [-PI, PI)
	* \param t1 es the angle 1
	* \param t2 is the angle 2
	* \return the value of the computed angle
	*/
float dFi(float t1, float t2);


template<class T>
unsigned int binaryPosSearch(const T *v, unsigned int size, const T &value, bool descending=true)
{
	unsigned int pos = size/2;
	
	if(size == 0)
		return 0;
	else if(size == 1)
	{
		if(value < v[0])
			return 1;
		else
			return 0;
	}
	else if(value < v[pos])
	{
		// Insert right
		return pos+1 + binaryPosSearch(v+pos+1, size-pos-1, value, descending);
	}
	else if(v[pos] < value)
	{
		// Insert left
		return binaryPosSearch(v, pos, value, descending);
	}
	else
		return pos;
}


std::vector<std::string> stringSplit(const std::string &source, const char *delimiter = ",", bool keepEmpty = false);
std::vector<double> stringSplitDouble(const std::string &source, const char *delimiter = ",", bool keepEmpty = false);





/**
	* Returns one string for each line contained in file \p filename.
	* \param filename Name of the file
	* \return Vector of lines contained in file \p filename
	*/
std::vector<std::string> readFileByLines(const std::string &filename);


/**
	* Reads the whole file as a single character string, removing any end of lines.
	* \param filename Name of the file
	* \param result Output value
	* \return Chars contained in file \p filename
	*/
void readFileByChars(const std::string &filename, char * &result);

#endif
