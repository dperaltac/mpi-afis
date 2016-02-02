#include<cmath>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<string>
#include<cstring>
#include<vector>
#include<algorithm>

#include "Functions.h"
#include "Constants.h"

using namespace std;

int roundInt(float r){

    return (r > 0.0) ? (int)floor(r + 0.5) : (int)ceil(r - 0.5);
}

float computeInnerAngle(int x11,int y11,int x12,int y12,int x21,int y21,int x22,int y22){

    float angle;
    int x1,x2,x3,y1,y2,y3;

    x3=x22-(x21-x11);
    y3=y22-(y21-y11);
    x1=x11;
    y1=y11;
    x2=x12;
    y2=y12;

    float dx21 = x2-x1;
    float dx31 = x3-x1;
    float dy21 = y2-y1;
    float dy31 = y3-y1;
    float m12 = sqrt( square(dx21) + square(dy21) );
    float m13 = sqrt( square(dx31) + square(dy31) );

    float a = (dx21*dx31 + dy21*dy31);
    float b = (m12 * m13);
    float c = a/b;

    if(c<=-1){
        return PI;
    }
    if(c>=-1){
        return 0;
    }

    angle = acos( a / b );

    return angle;
}

float computeAngle(int x1,int x2,int y1,int y2){

    float angle = atan2(y1-y2,x1-x2) * RADTOREG;

    if (angle < 0)
        return angle + 360;
		else
        return angle;
}

float computeAngleRad(int x1,int x2,int y1,int y2){

    float angle = atan2(y1-y2,x1-x2);

    if (angle < 0)
        return angle + PI2;
    else
        return angle;
}

void DistanceFromLine(double cx, double cy, double ax, double ay ,
                                          double bx, double by, double &distanceSegment,
                                          double &distanceLine)
{

        //
        // find the distance from the point (cx,cy) to the line
        // determined by the points (ax,ay) and (bx,by)
        //
        // distanceSegment = distance from the point to the line segment
        // distanceLine = distance from the point to the line (assuming
        //                                        infinite extent in both directions
        //

/*

Subject 1.02: How do I find the distance from a point to a line?


    Let the point be C (Cx,Cy) and the line be AB (Ax,Ay) to (Bx,By).
    Let P be the point of perpendicular projection of C on AB.  The parameter
    r, which indicates P's position along AB, is computed by the dot product
    of AC and AB divided by the square of the length of AB:

    (1)    AC dot AB
        r = ---------
            ||AB||^2

    r has the following meaning:

        r=0      P = A
        r=1      P = B
        r<0      P is on the backward extension of AB
        r>1      P is on the forward extension of AB
        0<r<1    P is interior to AB

    The length of a line segment in d dimensions, AB is computed by:

        L = sqrt( (Bx-Ax)^2 + (By-Ay)^2 + ... + (Bd-Ad)^2)

    so in 2D:

        L = sqrt( (Bx-Ax)^2 + (By-Ay)^2 )

    and the dot product of two vectors in d dimensions, U dot V is computed:

        D = (Ux * Vx) + (Uy * Vy) + ... + (Ud * Vd)

    so in 2D:

        D = (Ux * Vx) + (Uy * Vy)

    So (1) expands to:

            (Cx-Ax)(Bx-Ax) + (Cy-Ay)(By-Ay)
        r = -------------------------------
                          L^2

    The point P can then be found:

        Px = Ax + r(Bx-Ax)
        Py = Ay + r(By-Ay)

    And the distance from A to P = r*L.

    Use another parameter s to indicate the location along PC, with the
    following meaning:
          s<0      C is left of AB
          s>0      C is right of AB
          s=0      C is on AB

    Compute s as follows:

            (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
        s = -----------------------------
                        L^2


    Then the distance from C to P = |s|*L.

*/

	double diffbax = bx-ax;
	double diffbay = by-ay;
	double diffcax = cx-ax;
	double diffcay = cy-ay;
	double diffcbx = cx-bx;
	double diffcby = cy-by;
	double r_numerator = diffcax*diffbax + diffcay*diffbay;
	double r_denomenator = square(diffbax) + square(diffbay);

	distanceLine = fabs(diffcax*diffbay - diffcay*diffbax) / sqrt(r_denomenator);

//
// (xx,yy) is the point on the lineSegment closest to (cx,cy)
//

	if ( (r_numerator >= 0) && (r_numerator <= r_denomenator) )
	{
		distanceSegment = distanceLine;
	}
	else
	{
		double dist1 = square(diffcax) + square(diffcay);
		double dist2 = square(diffcbx) + square(diffcby);

		if (dist1 < dist2)
			distanceSegment = sqrt(dist1);
		else
			distanceSegment = sqrt(dist2);
	}
}


double DistanceFromSegment(double cx, double cy, double ax, double ay ,
                                          double bx, double by)
{
	double diffbax = bx-ax;
	double diffbay = by-ay;
	double diffcax = cx-ax;
	double diffcay = cy-ay;
	double diffcbx = cx-bx;
	double diffcby = cy-by;
	double r_numerator = diffcax*diffbax + diffcay*diffbay;
	double r_denomenator = square(diffbax) + square(diffbay);

	if ( (r_numerator >= 0) && (r_numerator <= r_denomenator) )
	{
		return fabs(diffcax*diffbay - diffcay*diffbax) / sqrt(r_denomenator);
	}
	else
	{
		double dist1 = square(diffcax) + square(diffcay);
		double dist2 = square(diffcbx) + square(diffcby);

		return sqrt(min(dist1, dist2));
	}
}

short shortIntegerSwap( short s )
{
	unsigned char b1, b2;
	
	b1 = s & 255;
	b2 = (s >> 8) & 255;
	
	return (b1 << 8) + b2;
}


int longIntegerSwap (int i)
{
	unsigned char b1, b2, b3, b4;
	
	b1 = i & 255;
	b2 = ( i >> 8 ) & 255;
	b3 = ( i>>16 ) & 255;
	b4 = ( i>>24 ) & 255;
	
	return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
}


float dFi(float t1, float t2)
{
	float diff = t1 - t2;

	if (diff > PI)
		return (PI2 - diff);
	else if (diff <= -PI)
		return (PI2 + diff);
	else
		return diff;
}

vector<string> stringSplit(const string &source, const char *delimiter, bool keepEmpty)
{
	vector<string> results;

	size_t prev = 0;
	size_t next = 0;

	while ((next = source.find_first_of(delimiter, prev)) != std::string::npos)
	{
		if (keepEmpty || (next - prev != 0))
			results.push_back(source.substr(prev, next - prev));
		
		prev = next + 1;
	}

	if (prev < source.size())
		results.push_back(source.substr(prev));

	return results;
}

vector<double> stringSplitDouble(const string &source, const char *delimiter, bool keepEmpty)
{
	vector<string> aux = stringSplit(source, delimiter, keepEmpty);
	vector<double> results(aux.size());
	
	for(unsigned int i = 0; i < aux.size(); ++i)
		results[i] = atof(aux[i].c_str());

	return results;
}


vector<string> readFileByLines(const string &filename)
{
	vector<string> res;
	string str;
	ifstream filelist;

	filelist.open(filename.c_str());

	if(filelist)
	{
		filelist >> str;

		// All strings in the list are read
		while(filelist)
		{
			res.push_back(str);
			filelist >> str;
		}

		filelist.close();
	}

	return res;
}



void readFileByChars(const string &filename, char * &result)
{
	ifstream t(filename.c_str());
	stringstream buffer;
	buffer << t.rdbuf();
	
	string str = buffer.str();
	
	str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
	
	result = new char[str.size()+1];
	
	strcpy(result, str.c_str());
}
