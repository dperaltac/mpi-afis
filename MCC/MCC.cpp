/**
 * \file    MCC.cpp
 * \author  Salvador García <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * The implementation of the MCC class.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <unistd.h>

#include "Constants.h"
#include "Minutia.h"
#include "MCC.h"
#include "Functions.h"
#include "Matrix.h"
#include "Munkres.h"


MCC::typeConsolidation MCC::consolidation;
bool MCC::convexhull;
unsigned int MCC::Ns;
bool MCC::bit;
const string MCC::valid_args = "-N:HC:B";

using namespace std;

MCC::MCC() : Fingerprint()
{
}

MCC::MCC(const MCC &fp) : Fingerprint(fp)
{
	cylinders = fp.cylinders;
}


int MCC::configureAlgorithm(unsigned int ns, int cons, bool ch, int fus, bool pbit)
{
// 	Fingerprint::configureAlgorithm();

	Ns = ns;
	consolidation = (typeConsolidation)cons;
	convexhull = ch;
	bit = pbit;

	Cylinder::configureAlgorithm(Ns, bit, cons == NHS);

	return 0;
}


int MCC::configureAlgorithm(int argc, char *argv[])
{
	Fingerprint::configureAlgorithm(argc, argv);
	int c;

	// getopt variables initialization
	optind = 1;
	opterr = 0;

	// Default values
	Ns = 8;
	consolidation = LSSR;
	convexhull = true;
	bit = false;

	while((c = getopt (argc, argv, valid_args.c_str())) != -1)
	{
		switch (c)
		{
			case 'N':
				Ns = atoi(optarg);
				break;
			case 'C':
				consolidation = getConsolidationType(optarg);
				break;
			case 'B':
				bit = true;
				break;
			case 'H':
				convexhull = false;
				break;
			case '?':
				if (valid_args.find(optopt) != string::npos)
					cerr << "Syntax error: option " << (char)optopt << " requires an argument. Using default." << endl;
		}
	}

	Cylinder::configureAlgorithm(Ns, bit, consolidation == NHS);

	return 0;
}


float MCC::match(const Fingerprint &f) const
{
	MCC &p = (MCC &)f;

	unsigned int num_cylinders = cylinders.size();
	unsigned int pnum_cylinders = p.cylinders.size();
	unsigned int n_P = computeNP(num_cylinders, pnum_cylinders);

	if(n_P > num_cylinders || n_P > pnum_cylinders)
		return -1;

	if(consolidation == NHS)
	{
		float sum = 0.0, elem, maxelem;

		for (unsigned int i=0; i<num_cylinders;i++)
		{
			maxelem = 0.0;

			for (unsigned int j=0; j<pnum_cylinders;j++)
			{
				elem = cylinders[i].nhs(p.cylinders[j]);

				if(maxelem < elem)
					maxelem = elem;
			}

			sum += maxelem;
		}

		return sum / num_cylinders;
	}

	Matrix <float> gamma (num_cylinders, pnum_cylinders);
	vector <pair <int,int> > pairs(n_P);

	for (unsigned int i=0; i<num_cylinders;i++)
		for (unsigned int j=0; j<pnum_cylinders;j++){
			gamma(i,j) = cylinders[i].similarity(p.cylinders[j]);
// 			cout << cylinders[i].getIndex() << "\t" << p.cylinders[j].getIndex() << "\t" << gamma(i,j) << endl;
		}

// 	cout << "GAMMA MATRIX:" << gamma << endl;

	switch (consolidation)
	{
		case LSS:
			consolidationLSS (pairs, gamma, n_P, num_cylinders, pnum_cylinders);
			break;
		case LSA:
			consolidationLSA (pairs, gamma, n_P, num_cylinders, pnum_cylinders);
			break;
		case LSSR:
			consolidationLSSR(pairs, gamma, n_P, cylinders,     p.cylinders   );
			break;
		case LSAR:
			consolidationLSAR(pairs, gamma, n_P, cylinders,     p.cylinders   );
			break;
		case LGS:
			consolidationLGS (pairs, gamma, n_P, num_cylinders, pnum_cylinders);
			break;
		case LGSR:
			consolidationLGSR(pairs, gamma, n_P, cylinders,     p.cylinders   );
			break;
		default:
			cerr << "MCC::match(): " << consolidation << " is not a valid consolidation" << endl;
			return -1;
	}


// 	for(unsigned int i = 0; i < pairs.size(); ++i)
// 		cout << "(" << cylinders[pairs[i].first].getIndex() << ",\t" << p.cylinders[pairs[i].second].getIndex() << ",\t" << gamma(pairs[i].first, pairs[i].second) << ") " << endl;
// 	cout << endl;

	return computeGlobalScore(pairs, gamma);
}


MCC::typeConsolidation MCC::getConsolidationType(const string &cons)
{
	if(cons == "LSS")
		return LSS;
	else if(cons == "LSSR")
		return LSSR;
	else if(cons == "LSA")
		return LSA;
	else if(cons == "LSAR")
		return LSAR;
	else if(cons == "LGS")
		return LGS;
	else if(cons == "LGSR")
		return LGSR;
	else if(cons == "NHS")
		return NHS;
	else
	{
		cerr << "WARNING: consolidation expected to be one of {LSS|LSA|LSSR|LSAR|LGS|LGSR|NHS} and found " << cons << " instead. Using LSSR as default." << endl;
		return LSSR;
	}
}


void MCC::computeConvexHull(vector <point2d> &convex)
{
	vector <point2d> pointSet(minutiae.size());

	convex.clear();

	for (unsigned int i=0; i<minutiae.size(); i++) {
			pointSet[i].x = minutiae[i].getX();
			pointSet[i].y = minutiae[i].getY();
	}

	GrahamScanConvexHull()(pointSet, convex);
}

void MCC::computeCylindersNHS()
{
	unsigned int num_minutiae = minutiae.size();

	cylinders.resize(num_minutiae);

	for (unsigned int ii=0; ii<num_minutiae;ii++)
		cylinders[ii].setMinutia(minutiae[ii]);

	for (unsigned int ii=0; ii<num_minutiae; ii++)
		cylinders[ii].computeNHSBitVectors(cylinders);
}


void MCC::computeCylinders()
{
	unsigned int num_minutiae = minutiae.size();

	cylinders.resize(num_minutiae);

	for (unsigned int ii=0; ii<num_minutiae;ii++)
		cylinders[ii].setMinutia(minutiae[ii]);

	if(convexhull)
	{
		vector<point2d> convex;

		computeConvexHull(convex);

		for (unsigned int ii=0; ii<num_minutiae; ii++)
			cylinders[ii].computeCMVector(cylinders, convex);
	}
	else
	{
		for (unsigned int ii=0; ii<num_minutiae; ii++)
			cylinders[ii].computeCMVector(cylinders);
	}
}

void MCC::compactCylinders()
{
	vector<Cylinder>::iterator i=cylinders.begin();

	while(i != cylinders.end())
	{
		if (i->getValidity())
			++i;
		else
			i = cylinders.erase(i);
	}
}

float MCC::computeGlobalScore(const vector <pair <int,int> > & elements, Matrix <float> & gamma) const
{
	float sum = 0.0;

	for (vector<pair <int,int> >::const_iterator i=elements.begin(); i != elements.end(); ++i)
			sum += gamma(i->first, i->second);

	return sum / elements.size();
}

float MCC::rho(Cylinder const & t_a, Cylinder const & t_b, Cylinder const & k_a, Cylinder const & k_b)
{
	float d1, d2, d3;

	d1 = fabs(t_a.ds(k_a.getX(),k_a.getY()) - t_b.ds(k_b.getX(),k_b.getY()));
	d2 = fabs(t_a.dFi(t_a.dFi(t_a.getrT(),k_a.getrT()),t_b.dFi(t_b.getrT(),k_b.getrT())));
	d3 = fabs(t_a.dFi(dR(t_a,k_a),dR(t_b,k_b)));

	return Cylinder::psi(d1,MUP1,TAUP1)*Cylinder::psi(d2,MUP2,TAUP2)*Cylinder::psi(d3,MUP3,TAUP3);
}

float MCC::rho_tun(Cylinder const & t_a, Cylinder const & t_b, Cylinder const & k_a, Cylinder const & k_b)
{
	float d1, d2, d3;

	d1 = fabs(t_a.ds(k_a.getX(),k_a.getY()) - t_b.ds(k_b.getX(),k_b.getY()))/(t_a.ds(k_a.getX(),k_a.getY()) + t_b.ds(k_b.getX(),k_b.getY()));
	d2 = fabs(t_a.dFi(t_a.dFi(t_a.getrT(),k_a.getrT()),t_b.dFi(t_b.getrT(),k_b.getrT())));
	d3 = fabs(t_a.dFi(dR(t_a,k_a),dR(t_b,k_b)));

	return Cylinder::psi(d1,MUP1,TAUP1)*Cylinder::psi(d2,MUP2,TAUP2)*Cylinder::psi(d3,MUP3,TAUP3);
}

void MCC::consolidationLSS(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, unsigned int nA, unsigned int nB) const
{
	unsigned int k;
	float scores[nP];

	best.resize(nP);

	for(unsigned int i = 0; i < nP; i++)
		scores[i] = -1.0;

	for (unsigned int i=0; i<nA; i++)
	{
		for (unsigned int j=0; j<nB; j++)
		{
			for (k=0; k<nP && scores[k]>=gamma(i,j); k++);

			if (k < nP)
			{
				for (unsigned int l=nP-1; l>k; l--)
				{
					scores[l] = scores[l-1];
					best[l] = best[l-1];
				}
				scores[k] = gamma(i,j);
				best[k].first = i;
				best[k].second = j;
			}
		}
	}
// 	cout << "Valid cylinders: " << nA << "," << nB << "," << nP << endl;
// 	for(unsigned int i = 0; i < best.size(); ++i)
// 		cout << i << "\t(" << best[i].first << ",\t" << best[i].second << ",\t" << scores[i] << ") " << endl;
// 	cout << endl;
}

void MCC::consolidationLSA(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, unsigned int nA, unsigned int nB) const
{
    unsigned int i,j,k,l;
    pair <int,int> tmp;
    float *scores = new float[nP];

    best.reserve(nP);
    best.assign(nP, tmp);
    for (i=0; i<nP; i++) {
        scores[i] = -1;
    }

    Matrix<float> matrix(nA, nB);

    for (i=0; i<nA; i++) {
        for (j=0; j<nB; j++) {
            if (gamma(i,j) > 0)
                matrix(i,j) = 1.0 / gamma(i,j);
            else
                matrix(i,j) = 666;
        }
    }

    // Apply Munkres algorithm to matrix.
    Munkres m;
    m.solve(matrix);

    for (i=0; i<nA; i++) {
        for (j=0; j<nB; j++) {
            if (matrix(i,j) == 0) {
                for (k=0; k<nP && scores[k]>=gamma(i,j); k++);
                if (k < nP) {
                    for (l=nP-1; l>k; l--) {
                        scores[l] = scores[l-1];
                        best[l].first = best[l-1].first;
                        best[l].second = best[l-1].second;
                    }
                    scores[k] = gamma(i,j);
                    best[k].first = i;
                    best[k].second = j;
                }
            }
        }
    }

    delete [] scores;
}

void MCC::consolidationLGS(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, unsigned int nA, unsigned int nB) const
{
		unsigned int k;
    float scores[nP];
	bool rowsSelected [nA];
	bool colsSelected [nB];

    best.resize(nP);

		for(unsigned int i = 0; i < nP; i++)
			scores[i] = -1.0;

	for (unsigned int i=0; i<nA; i++)
		rowsSelected[i] = true;
	for (unsigned int i=0; i<nB; i++)
		colsSelected[i] = true;


    for (unsigned int i=0; i<nA; i++) {
        for (unsigned int j=0; j<nB; j++) {
            for (k=0; k<nP && scores[k]>=gamma(i,j); k++);

            if (k < nP && rowsSelected[i] && colsSelected[j]) {
                for (unsigned int l=nP-1; l>k; l--) {
                    scores[l] = scores[l-1];
                    best[l] = best[l-1];
                }
                scores[k] = gamma(i,j);
                best[k].first = i;
                best[k].second = j;
				rowsSelected[i] = false;
				colsSelected[j] = false;
            }
        }
    }
}

void MCC::consolidationLSSR(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, vector <Cylinder> const & cA, vector <Cylinder> const & cB) const
{
	unsigned int nR = min(cA.size(), cB.size());
	float efficiency;
	float sum;
	const float LAMBDAWEIGHT = (1.0-WR)/(nR-1.0);
	vector <pair <int,int> > inter(nR);
	vector<float> scoresB(nP, -1);
	vector<float> lambdaT(nR);
	vector<float> lambdaT1(nR);

	Matrix<float> rhotab(nR,nR);

	best.resize(nP);

	consolidationLSS(inter, gamma, nR, cA.size(), cB.size());

// 	for(unsigned int i = 0; i < inter.size(); ++i)
// 		cout << i << "\t(" << cA[inter[i].first].getIndex() << ",\t" << cB[inter[i].second].getIndex() << ",\t" << gamma[inter[i].first][inter[i].second] << ") " << endl;
// 	cout << endl;

	for (unsigned int i=0; i<nR; i++)
	{
		lambdaT[i] = gamma(inter[i].first,inter[i].second);

		for (unsigned int k=0; k<nR; k++)
			if (k!=i)
				rhotab(i, k) = rho(cA[inter[i].first],cB[inter[i].second],cA[inter[k].first],cB[inter[k].second]);
			else
				rhotab(i, k) = 0;
	}

// 	cout << "lambdaT: " << lambdaT.size() << endl;
// cout << LAMBDAWEIGHT << endl;

	for (unsigned int i=0; i<NREL; i++)
	{
		lambdaT.swap(lambdaT1);

		for (unsigned int j=0; j<nR; j++)
		{
			sum = 0.0;
			for (unsigned int k=0; k<nR; k++)
				sum += rhotab(j,k) * lambdaT1[k];

			lambdaT[j] = WR*lambdaT1[j] + LAMBDAWEIGHT*sum;
		}
	}

	unsigned int j;
	for (unsigned int i=0; i<nR; i++)
	{
		efficiency = lambdaT[i] / gamma(inter[i].first,inter[i].second);

		for (j=0; j<nP && scoresB[j]>=efficiency; j++);

		if (j < nP)
		{
			for (unsigned int k=nP-1; k>j;k--)
			{
				scoresB[k] = scoresB[k-1];
				best[k] = best[k-1];
			}

			scoresB[j] = efficiency;
			best[j] = inter[i];
			gamma(inter[i].first,inter[i].second) = lambdaT[i];
		}
	}

}

void MCC::consolidationLSAR(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, vector <Cylinder> const & cA, vector <Cylinder> const & cB) const
{
    vector <pair <int,int> > inter;
    unsigned int nR = min(cA.size(),cB.size());
    float *scoresI = new float[nR];
    float *scoresB = new float[nP];
    float *lambdaT = new float[nR];
    float *efficiency = new float[nR];
    unsigned int i, j, k;
    pair <int,int> tmp;
    float sum;

    inter.reserve(nR);
    inter.assign(nR, tmp);
    best.reserve(nP);
    best.assign(nP, tmp);
    for (i=0; i<nP; i++) {
        scoresB[i] = -1;
    }

    consolidationLSA(inter, gamma,nR,(int)cA.size(),(int)cB.size());
    for (i=0; i<nR; i++) {
        scoresI[i] = gamma(inter[i].first,inter[i].second);
        lambdaT[i] = scoresI[i];
    }

    for (i=0; i<NREL; i++) {
        for (j=0; j<nR; j++) {
            sum = 0;
            for (k=0; k<nR; k++) {
                if (k!=j) {
                    sum += rho(cA[inter[j].first],cB[inter[j].second],cA[inter[k].first],cB[inter[k].second]) * lambdaT[k];
                }
            }
            lambdaT[j] = WR*lambdaT[j] + (1.0-WR)*sum/((float)nR-1.0);
        }
    }

    for (i=0; i<nR; i++) {
        efficiency[i] = lambdaT[i] / scoresI[i];
    }

    for (i=0; i<nR; i++) {
        for (j=0; j<nP && scoresB[j]>=efficiency[i]; j++);
        if (j < nP) {
            for (k=nP-1; k>j;k--) {
                scoresB[k] = scoresB[k-1];
                best[k].first = best[k-1].first;
                best[k].second = best[k-1].second;
            }
            scoresB[j] = efficiency[i];
            best[j].first = inter[i].first;
            best[j].second = inter[i].second;
        }
    }

    for (i=0; i<nR; i++) {
        gamma(inter[i].first,inter[i].second) = lambdaT[i];
    }

    delete [] scoresI;
    delete [] scoresB;
    delete [] lambdaT;
    delete [] efficiency;
}

void MCC::consolidationLGSR(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, vector <Cylinder> const & cA, vector <Cylinder> const & cB) const
{
    unsigned int nR = min(cA.size(), cB.size());
    float efficiency;
    unsigned int j;
    float sum;
		const float LAMBDAWEIGHT = (1.0-WR)/(nR-1.0);
    vector <pair <int,int> > inter(nR);
    vector<float> scoresB(nP, -1);
    vector<float> lambdaT(nR);
    vector<float> lambdaT1(nR);

		best.resize(nP);

    consolidationLGS(inter, gamma, nR, cA.size(), cB.size());

    for (unsigned int i=0; i<nR; i++)
        lambdaT[i] = gamma(inter[i].first,inter[i].second);


    for (unsigned int i=0; i<NREL; i++) {
			lambdaT.swap(lambdaT1);

        for (j=0; j<nR; j++) {
            sum = 0.0;
            for (unsigned int k=0; k<nR; k++) {
                if (k!=j) {
                    sum += rho_tun(cA[inter[j].first],cB[inter[j].second],cA[inter[k].first],cB[inter[k].second]) * lambdaT1[k];
                }
            }
            lambdaT[j] = WR*lambdaT1[j] + LAMBDAWEIGHT*sum;
        }
    }

    for (unsigned int i=0; i<nR; i++) {
        efficiency = lambdaT[i] / gamma(inter[i].first,inter[i].second);
        for (j=0; j<nP && scoresB[j]>=efficiency; j++);
        if (j < nP) {
            for (unsigned int k=nP-1; k>j;k--) {
                scoresB[k] = scoresB[k-1];
                best[k] = best[k-1];
            }
            scoresB[j] = efficiency;
            best[j] = inter[i];
            gamma(inter[i].first,inter[i].second) = lambdaT[i];
        }
    }
}

void MCC::addStructureBD(){

}

void MCC::loadCylinder(string name){
    ifstream reader;
    int nCylinders;
    float * bufferF  =  new float [Ns*Ns*Nd];
    bool * bufferB1 = new bool[Ns*Ns*Nd];
    bool * bufferB2 = new bool[Ns*Ns*Nd];
    int auxI;
    float auxF;


	cerr << "Leyendo cilindro" << endl;

    reader.open(name.c_str(), fstream::binary);
    if (!reader.is_open())
    {
       cerr << "Error opening file" << endl;
       return;
    }



        reader.read ((char *) (&nCylinders), sizeof (int));


	//cout << "Tenemos cilindros = " << nCylinders << endl;

	Cylinder cmVectorAUX;
	for(int i=0; i<nCylinders; ++i){

            //read index
            reader.read ((char *) (&auxI), sizeof (int));
            cmVectorAUX.setIndex(auxI);

            //read X Y T variables
            reader.read ((char *) (&auxI), sizeof (int));
            cmVectorAUX.setX(auxI);

            reader.read ((char *) (&auxI), sizeof (int));
            cmVectorAUX.setY(auxI);

            reader.read ((char *) (&auxF), sizeof (float));
            cmVectorAUX.setT(auxF);



	    if(bit){
	        reader.read ((char *) bufferB1, sizeof (bool)*Ns*Ns*Nd);
	        reader.read ((char *) bufferB2, sizeof (bool)*Ns*Ns*Nd);

		for(unsigned int j=0; j< Ns*Ns*Nd;++j){
			cmVectorAUX.setB1(j, bufferB1[j]);
			cmVectorAUX.setB2(j, bufferB2[j]);
		}

	    }else{

	        //bufferF
	        reader.read ((char *) bufferF, sizeof (float)*Ns*Ns*Nd);

		//cmVectorAUX.reserve(Ns*Ns*Nd);

		for(unsigned int j=0; j< Ns*Ns*Nd;++j){
		//	cout << bufferF[j] << "  ";
			cmVectorAUX.setCM(j, bufferF[j]);
			cmVectorAUX.setValidity(true);
		}


		//	delete [] bufferF;

	    }

		//cout << endl;
		cylinders.push_back(cmVectorAUX);

	}  // endfor

		//cout << cylinders[22];
	cout << "Fin lectura = " << (int)cylinders.size() << endl;

    reader.close();
}

void MCC::loadCylinder(ifstream &reader){

    int nCylinders;
    float * bufferF  =  new float [Ns*Ns*Nd];
    int auxI;
    float auxF;
    bool * bufferB1 = new bool[Ns*Ns*Nd];
    bool * bufferB2 = new bool[Ns*Ns*Nd];

//cylinders.clear();
        reader.read ((char *) (&nCylinders), sizeof (int));


	//cout << "Tenemos cilindros = " << nCylinders << endl;

		Cylinder cmVectorAUX;
	for(int i=0; i<nCylinders; ++i){

            //read index
            reader.read ((char *) (&auxI), sizeof (int));
            cmVectorAUX.setIndex(auxI);

            //read X Y T variables
            reader.read ((char *) (&auxI), sizeof (int));
            cmVectorAUX.setX(auxI);

            reader.read ((char *) (&auxI), sizeof (int));
            cmVectorAUX.setY(auxI);

            reader.read ((char *) (&auxF), sizeof (float));
            cmVectorAUX.setT(auxF);




	    if(bit){
	        reader.read ((char *) bufferB1, sizeof (bool)*Ns*Ns*Nd);
	        reader.read ((char *) bufferB2, sizeof (bool)*Ns*Ns*Nd);

		for(unsigned int j=0; j< Ns*Ns*Nd;++j){
			cmVectorAUX.setB1(j, bufferB1[j]);
			cmVectorAUX.setB2(j, bufferB2[j]);
		}

	    }else{

	        //bufferF
	        reader.read ((char *) bufferF, sizeof (float)*Ns*Ns*Nd);

		//cmVectorAUX.reserve(Ns*Ns*Nd);

		for(unsigned int j=0; j< Ns*Ns*Nd;++j){
		//	cout << bufferF[j] << "  ";
			cmVectorAUX.setCM(j, bufferF[j]);
			cmVectorAUX.setValidity(true);
		}


		//	delete [] bufferF;

	    }

		//cout << endl;
		cylinders.push_back(cmVectorAUX);
	}

		//cout << cylinders[22];
	//cout << "Fin lectura = " << (int)cylinders.size() << endl;


}


void MCC::writeCylinder(string name){
    ofstream writer;
    string id;
    Minutia min;
    int indexMin;
    int XMin;
    int YMin;
    float TMin;
    float * floatBuffer;
    bool * boolBuffer, * boolBuffer2;

    writer.open(name.c_str(), fstream::binary | fstream::app);  // at the end of the file. | fstream::app

    if (!writer.is_open())
    {
       cerr << "Error opening file" << endl;
       return;
    }




// save the cylinders

		int nCylinders = (int)cylinders.size();

		writer.write ((const char *) (&nCylinders), sizeof (nCylinders));

		//cout<< "Numer of cilinders = " << nCylinders <<  endl;
		// For each cylinder.
		for(int i=0; i<nCylinders; ++i){

			    //write index
		    indexMin=cylinders[i].getIndex();
		    writer.write ((const char *) &indexMin, sizeof (indexMin));

		    //write X Y T variables
		    XMin=cylinders[i].getX();
		    writer.write ((const char *) &XMin, sizeof (XMin));

		    YMin=cylinders[i].getY();
		    writer.write ((const char *) &YMin, sizeof (YMin));

		    TMin=cylinders[i].getT();
		    writer.write ((const char *) &TMin, sizeof (TMin));


			if(bit){
				//cout << "soy de bit" << endl;
				int sizeCylinder= cylinders[i].getB1Vector().size();

				boolBuffer= new bool [sizeCylinder];
				boolBuffer2= new bool [sizeCylinder];

				for(int j=0; j< sizeCylinder;++j){
					boolBuffer[j]=cylinders[i].getB1(j);
					boolBuffer2[j]=cylinders[i].getB2(j);
				}

				writer.write ((const char *) boolBuffer, sizeof(bool)*sizeCylinder);
				writer.write ((const char *) boolBuffer2, sizeof(bool)*sizeCylinder);


			}else{
				//cout << "soy de reales" << endl;

				int sizeCylinder= cylinders[i].getCMVector().size();
				floatBuffer= new float [sizeCylinder];

				for(int j=0; j< sizeCylinder;++j){
					floatBuffer[j]=cylinders[i].getCM(j);
				}

				writer.write ((const char *) floatBuffer, sizeof(float)*sizeCylinder);

			}
		}

		//cout << cylinders[22];




    writer.close();

}


void MCC::initialize()
{
	if(consolidation == NHS)
		computeCylindersNHS();
	else
		computeCylinders();

	compactCylinders();
}


/**
 * Prints a fingerprint
 * \param output Output stream
 * \param F Fingerprint
 */
ostream& operator<<(ostream& output, const MCC& F)
{
	output << (Fingerprint &)F << endl;
	
	for(unsigned int i = 0; i < F.cylinders.size(); ++i)
		output << "Cylinder " << i << ": " << F.cylinders[i] << endl;
	output << endl;
	
	return output;
}
