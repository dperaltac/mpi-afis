/**
 * \file    Cylinder.cpp
 * \author  Salvador Garcia <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * The implementation of the cylinder class.
 */

#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>

#include "Cylinder.h"
#include "Functions.h"

unsigned int Cylinder::NUMCELLS;
float Cylinder::DELTAS;
float Cylinder::RELSHIFT;
unsigned int Cylinder::MINCELLS;
float Cylinder::MUPSI;
unsigned int Cylinder::Ns;
bool Cylinder::bit;

unsigned int Cylinder::MAXNVCELLS;

Cylinder::Cylinder() : Minutia(), valid(false), cmVector(), cmBit1(), cmBit2() {

	if(bit)
	{
		cmBit1.resize(NUMCELLS,false);
		cmBit2.resize(NUMCELLS,false);
	}
	else
		cmVector.resize(NUMCELLS,0);
}

Cylinder::~Cylinder(){
//     cmVector.clear();
//     cmBit1.clear();
//     cmBit2.clear();
}

Cylinder::Cylinder(int newIndex, int valX, int valY, float valT) : 
	Minutia(newIndex, valX, valY, valT, 0, RIG),
	valid(false), cmVector(), cmBit1(), cmBit2()
{
	if(bit)
	{
		cmBit1.resize(NUMCELLS,false);
		cmBit2.resize(NUMCELLS,false);
	}
	else
		cmVector.resize(NUMCELLS,0);
}

Cylinder::Cylinder(const Cylinder& other) : 
	Minutia(other),
	valid(other.valid),
	cmVector(other.cmVector),
	cmBit1(other.cmBit1),
	cmBit2(other.cmBit2)
{
}

Cylinder::Cylinder(const Minutia& other) : Minutia(other), valid(false), cmVector(), cmBit1(), cmBit2() {

	if(bit)
		{
			cmBit1.resize(NUMCELLS,false);
			cmBit2.resize(NUMCELLS,false);
		}
		else
			cmVector.resize(NUMCELLS,0);
}

int Cylinder::configureAlgorithm(unsigned int ns, bool pbit, bool nhs)
{
	Ns = ns;
	bit = pbit;
	
	DELTAS = (2.0*R)/Ns;
	RELSHIFT = 0.5*(Ns+1);
	
	if(nhs)
	{
		MUPSI = 0.005;
		NUMCELLS = 0;
		
		for (unsigned int i=1; i<=Ns; i++)
			for (unsigned int j=1; j<=Ns; j++)
				if(square(i-RELSHIFT) + square(j-RELSHIFT) <= square(Ns/2))
					NUMCELLS++;
			
		NUMCELLS *= Nd;
	}
	else
	{
		MUPSI = 0.01;
		NUMCELLS = Nd*Ns*Ns;
	}
	
	MINCELLS = floor(MINME*NUMCELLS);
	
	if(nhs)
		MAXNVCELLS = floor((1.0-MINVC) * NUMCELLS);
	else
		MAXNVCELLS = floor((1.0-MINVC*(2-4/PI)) * NUMCELLS);
	
	return 0;
}


int Cylinder::configureAlgorithm(int argc, char *argv[])
{
	int res;
	unsigned int ns;
	bool nhs;
	bool pbit;
	
	if(argc < 4)
	{
		cerr << "ERROR: expected 4 parameters for MCC fingerprints: <Ns> {LSS|LSA|LSSR|LSAR|NHS} {y|n} {b|r}. Using default." << endl;
		ns = 8;
		res = 1;
		pbit = false;
		nhs = false;
	}
	else
	{
		ns = atoi(argv[0]);
		nhs = (strcmp(argv[1], "NHS") == 0);
		
		if(argv[3][0] == 'b')
			pbit = true;
		else if(argv[3][0] == 'r')
			pbit = false;
		else
		{
			cerr << "ERROR: the bit/real parameter should be {b|r}. Using default (real encoding)." << endl;
			pbit = false;
		}
		
		res = 0;
	}
	
	configureAlgorithm(ns, pbit, nhs);
	
	return res;
}


void Cylinder::setMinutia(const Minutia& other){

    X=other.getX();
    Y=other.getY();
    T=other.getT();
    index=other.getIndex();
    valid = false;

		if(bit)
		{
			cmBit1.resize(NUMCELLS,false);
			cmBit2.resize(NUMCELLS,false);
		}
		else
			cmVector.resize(NUMCELLS,0);
}

Cylinder& Cylinder::operator=(const Cylinder& rhs){

    if (this == &rhs) return *this;

    X=rhs.X;
    Y=rhs.Y;
    T=rhs.T;
    index=rhs.index;
    cmVector=rhs.cmVector;
    valid=rhs.valid;
    cmBit1=rhs.cmBit1;
    cmBit2=rhs.cmBit2;
		bit = rhs.bit;
		Ns = rhs.Ns;

    return *this;
}

/**
 * Prints a cylinder
 * \param output Output stream
 * \param C Cylinder to print
 */
std::ostream& operator<<(std::ostream& output, const Cylinder& C) {

//     unsigned int i, j, k;

    if (C.bit == false) {
//         output << C.index << ": (" <<  C.X << ", " << C.Y << ", " << C.T << (C.valid?") Valid":") Invalid") << endl << "cmVector real: " << endl;
//         for (k=1; k<=C.Nd; k++) {
//             output << "Layer " << k << endl;
//             for (i=1; i<=C.Ns; i++) {
//                 for (j=1; j<=C.Ns; j++) {
//                     output << C.getCM(Cylinder::lin(i,j,k)) << " ";
//                 }
//                 output << endl;
//             }
//         }
			
			output << C.index << ": (" <<  C.X << ", " << C.Y << ", " << C.T << (C.valid?") Valid":") Invalid") << endl << "cmVector real: " << endl;
			
			for(unsigned int i=0; i < C.NUMCELLS; ++i)
				output << C.getCM(i) << " ";
			
			output << endl;
    } else {
//         output << C.index << ": (" <<  C.X << ", " << C.Y << ", " << C.T << (C.valid?") Valid":") Invalid") << endl << "cmVector bit: " << endl;
//         for (k=1; k<=C.Nd; k++) {
//             output << "Layer " << k << endl;
//             for (i=1; i<=C.Ns; i++) {
//                 for (j=1; j<=C.Ns; j++) {
//                     output << C.getB1(Cylinder::lin(i,j,k)) << "(" << C.getB2(Cylinder::lin(i,j,k)) << ") ";
//                 }
//                 output << endl;
//             }
//         }
    }
    return output;
}

void Cylinder::pij(int i, int j, float & i_out, float & j_out)
{
	float relative_i = i - RELSHIFT;
	float relative_j = j - RELSHIFT;
	float cosT = cos(getrT());
	float sinT = sin(getrT());

	i_out = DELTAS * (cosT*relative_i + sinT*relative_j) + this->X;
	j_out = DELTAS * (cosT*relative_j - sinT*relative_i) + this->Y;
}

bool Cylinder::xiM(float i, float j, vector <point2d> const & convex)
{
	double distanceSegment;

	if (dss(i,j) > R*R)
	{
		return false;
	}
	else if (pnpoly(convex,i,j))
	{
		return true;
	}
	else
	{
		unsigned int last = convex.size() - 1;
		
		distanceSegment = DistanceFromSegment(i,j,convex[last].x,convex[last].y,convex[0].x,convex[0].y);
		
		if (distanceSegment <= OMEGA)
			return true;
			
		for (unsigned int k=0; k<last; ++k)
		{
			distanceSegment = DistanceFromSegment(i,j,convex[k].x,convex[k].y,convex[k+1].x,convex[k+1].y);
			
			if (distanceSegment <= OMEGA)
				return true;
		}

		return false;
	}
}

bool Cylinder::xiM(float i, float j)
{
	return (dss(i,j) <= R*R);
}

bool Cylinder::pnpoly(vector <point2d> const & convex, float testx, float testy)
{
	int pivote;
	unsigned int last = convex.size()-1;
	bool right;

	pivote = (testy-convex[last].y)*(convex[0].x-convex[last].x) -
						(testx-convex[last].x)*(convex[0].y-convex[last].y);

	if(pivote == 0)
		return true;

	right = (pivote < 0);

	for (unsigned int i=0; i<last; i++) {
		pivote = (testy-convex[i].y)*(convex[i+1].x-convex[i].x) - (testx-convex[i].x)*(convex[i+1].y-convex[i].y);

		if ((pivote < 0 && !right) || (pivote > 0 && right))
			return false;
		else if(pivote == 0)
			return true;
	}

	return true;
}

float Cylinder::dFi(float ang1, float ang2)
{
    float diff = ang1-ang2;

    if (diff >= -PI && diff < PI) {
        return diff;
    } else if (diff < -PI) {
        return PI2 + diff;
    } else {
        return -PI2 + diff;
    }
}

float Cylinder::contributionD(float angle, float k_angle) const
{
	float alpha = dFi(dFi(angle, getrT()), k_angle) * ISIGMAD;

	/*Area of the curve computation*/
	return doLeft(alpha + 0.5) - doLeft(alpha - 0.5);
}



void Cylinder::computeCMVector(vector <Cylinder> const & cylinders, vector <point2d> const & convex)
{
	unsigned int count;
	float cmVal;
	
	float i_abs, j_abs;

	count = 0;
	
	if (bit==false)
	{
		for (unsigned int i=1, pos = 0; i<=Ns; i++)
			for (unsigned int j=1; j<=Ns; j++)
			{
				pij(i,j,i_abs,j_abs);
				
				if (xiM(i_abs,j_abs,convex))
				{
					for (unsigned int k=1; k<=Nd; k++, pos++)
						cmVector[pos] = cm(i_abs,j_abs,k,cylinders,convex);
				}
				else
				{
					for (unsigned int k=1; k<=Nd; k++, pos++)
						cmVector[pos] = -1.0;
					
					count += Nd;
				}
			}
	}
	else
	{
		for (unsigned int i=1, pos = 0; i<=Ns; i++)
			for (unsigned int j=1; j<=Ns; j++)
			{
				pij(i,j,i_abs,j_abs);
				
				if (xiM(i_abs,j_abs,convex))
				{
					for (unsigned int k=1; k<=Nd; k++, pos++)
					{
						cmVal = cm(i_abs,j_abs,k,cylinders,convex);
						
						cmBit1[pos] = (cmVal == 1.0);
						cmBit2[pos] = (cmVal != -1.0);
					}
				}
				else
				{
					for (unsigned int k=1; k<=Nd; k++, pos++)
					{
						cmBit1[pos] = false;
						cmBit2[pos] = false;
					}
					
					count += Nd;
				}
			}
	}

	if (count >= MAXNVCELLS)
	{
		valid = false;
	}
	else
	{
		count = 0;
		
		for (unsigned int jj=0; count < MINM && (count+cylinders.size()-jj) >= MINM && jj<cylinders.size(); jj++)
			if (cylinders[jj].dss(X, Y) <= RNEIGHBORHOODRADIUS*RNEIGHBORHOODRADIUS)
				count++;

		valid = (count >= MINM);
	}
}


void Cylinder::computeCMVector(vector <Cylinder> const & cylinders)
{
	unsigned int count;
	float cmVal;
	
	float i_abs, j_abs;

	count = 0;
	
	if (bit==false)
	{
		for (unsigned int i=1, pos = 0; i<=Ns; i++)
			for (unsigned int j=1; j<=Ns; j++)
			{
				pij(i,j,i_abs,j_abs);
				
				if (xiM(i_abs,j_abs))
				{
					for (unsigned int k=1; k<=Nd; k++, pos++)
						cmVector[pos] = cm(i_abs,j_abs,k,cylinders);
				}
				else
				{
					for (unsigned int k=1; k<=Nd; k++, pos++)
						cmVector[pos] = -1.0;
					
					count += Nd;
				}
			}
	}
	else
	{
		for (unsigned int i=1, pos = 0; i<=Ns; i++)
			for (unsigned int j=1; j<=Ns; j++)
			{
				pij(i,j,i_abs,j_abs);
				
				if (xiM(i_abs,j_abs))
				{
					for (unsigned int k=1; k<=Nd; k++, pos++)
					{
						cmVal = cm(i_abs,j_abs,k,cylinders);
						
						cmBit1[pos] = (cmVal == 1.0);
						cmBit2[pos] = (cmVal != -1.0);
					}
				}
				else
				{
					for (unsigned int k=1; k<=Nd; k++, pos++)
					{
						cmBit1[pos] = false;
						cmBit2[pos] = false;
					}
					
					count += Nd;
				}
			}
	}

	if (count >= MAXNVCELLS)
	{
		valid = false;
	}
	else
	{
		count = 0;
		
		for (unsigned int jj=0; count < MINM && (count+cylinders.size()-jj) >= MINM && jj<cylinders.size(); jj++)
			if (cylinders[jj].dss(X, Y) <= RNEIGHBORHOODRADIUS*RNEIGHBORHOODRADIUS)
				count++;

		valid = (count >= MINM);
	}
}

void Cylinder::computeNHSBitVectors(vector <Cylinder> const & cylinders)
{
	unsigned int count, pos;
	float cmVal;
	float i_abs, j_abs;

	count = 0;
	pos = 0;
	
	for (unsigned int i=1; i<=Ns; i++)
		for (unsigned int j=1; j<=Ns; j++)
			for (unsigned int k=1; k<=Nd; k++)
			{
				pij(i,j,i_abs,j_abs);
				
				if(xiM(i_abs,j_abs))
				{
					cmVal = cm(i,j,k,cylinders);
				
					setB1(pos, cmVal == 1.0);
			
					if (cmVal == -1.0)
						count++;
					
					pos++;
				}
			}

	if (count >= MAXNVCELLS)
	{
		setValidity(false);
	}
	else
	{
		count = 0;
		
		for (unsigned int jj=0; count < MINM && (count+cylinders.size()-jj) >= MINM && jj<cylinders.size(); jj++)
			if (cylinders[jj].dss(X, Y) <= RNEIGHBORHOODRADIUS*RNEIGHBORHOODRADIUS)
				count++;

		setValidity(count >= MINM);
	}
}


float Cylinder::cm(float i_abs, float j_abs, unsigned int k, vector <Cylinder> const & cylinders, vector <point2d> const & convex)
{
	float sum = 0.0;
	float angle = this->getrT();
	float dfikk = dFik(k);
		
	for (vector<Cylinder>::const_iterator ii=cylinders.begin(); ii != cylinders.end(); ++ii)
		if (this->index != ii->index && ii->dss(i_abs,j_abs) <= NEIGHBORHOODRADIUS*NEIGHBORHOODRADIUS)
			sum += ii->contributionS(i_abs,j_abs) * ii->contributionD(angle,dfikk);
		
	if(bit)
		return Cylinder::psiBit(sum);
	else
		return Cylinder::psi(sum,MUPSI,TAUPSI);
}

float Cylinder::cm(float i_abs, float j_abs, unsigned int k, vector <Cylinder> const & cylinders)
{
	float sum = 0.0;
	float angle = this->getrT();
	float dfikk = dFik(k);
	
	for (vector<Cylinder>::const_iterator ii=cylinders.begin(); ii != cylinders.end(); ++ii)
		if (this->index != ii->index && ii->dss(i_abs,j_abs) <= NEIGHBORHOODRADIUS*NEIGHBORHOODRADIUS)
			sum += ii->contributionS(i_abs,j_abs) * ii->contributionD(angle,dfikk);
	
	if(bit)
		return Cylinder::psiBit(sum);
	else
		return Cylinder::psi(sum,MUPSI,TAUPSI);
}


float Cylinder::nhs(const Cylinder &c) const
{
	int hamming = 0;
	const float P = 30;

	if (abs(dFi(getrT(),c.getrT())) > 0.785398163397 || dss(c.getX(), c.getY()) > DELTAXY*DELTAXY)
		return 0;
	
	for(unsigned int i = 0; i < NUMCELLS; ++i)
		if (getB1(i) != c.getB1(i))
			hamming++;
		
	return pow(1.0-((float)hamming)/NUMCELLS, P);
}


float Cylinder::similarity(const Cylinder & c) const
{
	unsigned int count = 0;
	float norma_b = 0, normb_a = 0, norm_diff = 0;
	float ca_b, cb_a;

	if (abs(dFi(getrT(),c.getrT())) > DELTAZETA)
		return 0;
	
	if (bit == false)
	{
		for (unsigned int i=0; i<NUMCELLS; i++)
		{
			ca_b = cmVector[i];
			cb_a = c.getCM(i);

			if (ca_b>=0 && cb_a>=0)
			{
				count++;

				norma_b += ca_b*ca_b;
				normb_a += cb_a*cb_a;
				norm_diff += ca_b*cb_a;
			}
		}

		//Check if two cylinders are matchable
		if(count >= MINCELLS)
		{
			norm_diff = sqrt(norma_b + normb_a - 2.0*norm_diff);
			return 1.0 - (norm_diff/(sqrt(norma_b)+sqrt(normb_a)));
		}
		else
			return 0;
	}
	else
	{
		int counta_b = 0, countb_a = 0, count_diff = 0;

		for (unsigned int i=0; i<cmBit1.size(); ++i)
			if (getB2(i) && c.getB2(i))
			{
				count++;
				
				if (getB1(i))
				{
					counta_b++;
					
					if(!c.getB1(i))
						count_diff++;
				}
				
				if (c.getB1(i))
				{
					countb_a++;
					
					if(!getB1(i))
						count_diff++;
				}
			}
		
		//Check if two cylinders are matchable
		if (count >= MINCELLS)
			return (1 - (sqrt(count_diff)/(sqrt(counta_b)+sqrt(countb_a))));
		else
			return 0;
	}
}


