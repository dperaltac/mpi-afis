/**
 * \file    FingerprintJiang.cpp
 * \author  Salvador García <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * The implementation of the FingerprintJiang class.
 */

#include <algorithm>
#include <cmath>
#include "Parameters.h"
#include "Constants.h"
#include "Functions.h"
#include "Minutia.h"
#include "Matrix.h"
#include "FingerprintJiang.h"
#include <iomanip>
#include <iostream>
using namespace std;


float FingerprintJiang::match(const Fingerprint &f) const
{
	float accumulated, maxAccumulated = 0.0;
	unsigned int maxJ = 0, maxK = 0;
	float maxML;
	int posK;
	int coordx, coordy;
	float angle;

	FingerprintJiang &p = (FingerprintJiang &)f;
	unsigned int num_minutiae = minutiae.size(), pnum_minutiae = p.minutiae.size();
	
	float bestS[BEST];
	unsigned int bestI[BEST], bestJ[BEST];

	if(num_minutiae <= static_cast<unsigned int>(NEWFVSIZE) || pnum_minutiae <= static_cast<unsigned int>(NEWFVSIZE))
		return -1;

	Matrix<float> Fg1(num_minutiae, 3), Fg2(pnum_minutiae, 3);
	Matrix<float> sl(num_minutiae, pnum_minutiae), ml(num_minutiae, pnum_minutiae);

	for (int i=0; i<BEST; ++i) {
			bestS[i] = 0.0;
			bestJ[i] = 0;
			bestI[i] = 0;
	}

	// Similarity computation
	for (unsigned int i=0; i < num_minutiae; ++i) {
		for (unsigned int j=0; j < pnum_minutiae; ++j) {
			sl(i,j) = similarity(i,p,j);
			posK = BEST-1;

			while(posK>=0 && bestS[posK] < sl(i,j))
				--posK;

			// Keeping the BEST similar minutiae
			if (posK < BEST-1) {
				++posK;

				for (int l=BEST-1; l>posK; --l) {
						bestS[l]=bestS[l-1];
						bestI[l]=bestI[l-1];
						bestJ[l]=bestJ[l-1];
				}
				bestS[posK] = sl(i,j);
				bestI[posK] = i;
				bestJ[posK] = j;
			}
		}
	}

	for (int i=0; i<BEST; ++i) {
// 		cout << "Best match: " << bestI[i] << "," << bestJ[i] << endl;
// 	cout << "Coordinates: (" << minutiae[bestI[i]].getX() << "," << minutiae[bestI[i]].getY() << ") (" << p.minutiae[bestJ[i]].getX() << "," <<p.minutiae[bestJ[i]].getY() << ")" << endl;
// 	cout << "Spatial X diff: " << minutiae[bestI[i]].getX()-p.minutiae[bestJ[i]].getX() << endl;
// 	cout << "Spatial Y diff: " << minutiae[bestI[i]].getY()-p.minutiae[bestJ[i]].getY() << endl;
// 	cout << "Angular diff: " << minutiae[bestI[i]].getrT()-p.minutiae[bestJ[i]].getrT() << endl;
		coordx = minutiae[bestI[i]].getX();
		coordy = minutiae[bestI[i]].getY();
		angle  = minutiae[bestI[i]].getcrnT();

		// Convert to polar coordinates
		for (unsigned int j=0; j < num_minutiae; ++j) {
			Fg1(j,0) = distanceMatrix[j][bestI[i]];
			Fg1(j,1) = dFi(atan2(minutiae[j].getY() - coordy, minutiae[j].getX() - coordx), minutiae[bestI[i]].getcrnT());
			Fg1(j,2) = dFi(minutiae[j].getcrnT(), angle);
		}

		coordx = p.minutiae[bestJ[i]].getX();
		coordy = p.minutiae[bestJ[i]].getY();
		angle  = p.minutiae[bestJ[i]].getcrnT();

		// Convert to polar coordinates
		for (unsigned int j=0; j < pnum_minutiae; ++j) {
			Fg2(j,0) = p.distanceMatrix[j][bestJ[i]];
			Fg2(j,1) = dFi(atan2(p.minutiae[j].getY() - coordy, p.minutiae[j].getX() - coordx), p.minutiae[bestJ[i]].getcrnT());
			Fg2(j,2) = dFi(p.minutiae[j].getcrnT(), angle);
		}

		// Compute the ml matrix
		for (unsigned int j=0; j < num_minutiae; ++j)
			for (unsigned int k=0; k < pnum_minutiae; ++k)
				ml(j,k) = ML(Fg1[j],Fg2[k],sl[j][k]);

// 		// Avoid minutiae being doubly used for matching
		accumulated = 0.0;

		do
		{
			maxML = 0;

// 			Find the maximum value of the matrix that has not been explored yet
			for (unsigned int j=0; j < num_minutiae; ++j)
				for (unsigned int k=0; k < pnum_minutiae; ++k)
				{
					if(ml(j,k) > maxML)
					{
						maxJ = j;
						maxK = k;
						maxML = ml(maxJ, maxK);
					}
				}

			// If a value was found, set its whole row and column to 0
			if(maxML > 0)
			{
				accumulated += maxML;

				for (unsigned int l=0; l < num_minutiae; ++l)
					ml(l,maxK)=0.0;

				for (unsigned int l=0; l < pnum_minutiae; ++l)
					ml(maxJ,l)=0.0;
			}

		} while (maxML > 0);

		if (accumulated > maxAccumulated)
			maxAccumulated = accumulated;
	}

	return (maxAccumulated / std::max(num_minutiae, pnum_minutiae));
}

void FingerprintJiang::computeFeatureVectors()
{
	float angle;
	int coordx, coordy;
	int neighbor;

	if(minutiae.size() < NN+1)
		return;

	Flk.resize(minutiae.size(), NEWFVSIZE);

	for (unsigned int i=0; i < minutiae.size(); ++i) {

		angle = minutiae[i].getcrnT();
		coordx = minutiae[i].getX();
		coordy = minutiae[i].getY();

		for (int k = 0, kNN = NN; k<NN; ++k, ++kNN)
		{
			neighbor = neighbourhood[i][k];

			// Distances computation
// 			Flk(i,k) = distanceMatrix(i, neighbor);

			//radial angle computation
			Flk(i,k) = dFi(atan2(coordy - minutiae[neighbor].getY(),coordx - minutiae[neighbor].getX()), angle);

			//minutia direction computation
			Flk(i,kNN) = dFi(angle, minutiae[neighbor].getcrnT());

			//ridge count
// 			Flk(i,kNN3) = (float)ridgeCount(i,neighbor);

			// type of neighbour minutiae
// 			if (minutiae[neighbor].getType() == BIF)
// 				Flk(i,kNN4) = 0.0;
// 			else
// 				Flk(i,kNN4) = 1.0;
		}

// 		if (minutiae[i].getType() == BIF)
// 			Flk(i, NN4) = 0.0;
// 		else
// 			Flk(i, NN4) = 1.0;
	}
}

float FingerprintJiang::similarity (int i, const FingerprintJiang &p, int j) const
{
	float sum = 0.0;
	int k;

	for (k = 0; k<NN && sum < BL; ++k)
	{
		sum += abs(distanceMatrix(i,neighbourhood[i][k])-p.distanceMatrix(j,p.neighbourhood[j][k])) * W1
		     + abs(dFi(Flk(i,k),p.Flk(j,k))) * W2;
	}

	for (; k<NN2 && sum < BL; ++k)
		sum += abs(dFi(Flk(i,k),p.Flk(j,k))) * W3;

// 	for (; i<NN4 && sum < BL; ++i)
// 		sum += abs(Fli[i]-Flj[i]) * W4;
//
// 	for (; i<FVSIZE && sum < BL; ++i)
// 		sum += abs(Fli[i]-Flj[i]) * W5;

	if (sum < BL)
		return 1.0 - (sum / BL);
	else
		return 0.0;
}

float FingerprintJiang::ML (const float *Fgi, const float *Fgj, float slij) const
{
	if (abs(Fgi[0]-Fgj[0])<BG1 && abs(Fgi[1]-Fgj[1])<BG2 && abs(Fgi[2]-Fgj[2])<BG3)
		return 0.5*(1.0 + slij);
	else
		return 0.0;
}

void FingerprintJiang::initialize(){
	computeDistances();
	computeNeighbourhood();
// 	computeRidgeCount();
	computeFeatureVectors();
}
