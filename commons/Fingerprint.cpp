/**
 * \file    Fingerprint.cpp
 * \author  Joaquin Derrac <jderrac@decsai.ugr.es>
 * \author  Salvador Garcia <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * The implementation of the fingerprint class.
 */

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <unistd.h>

#include "Fingerprint.h"
#include "Constants.h"
#include "Functions.h"
#include "File19794.h"

using namespace std;

const string Fingerprint::valid_args = "";

Fingerprint::Fingerprint() :
	id(""),
	minutiae(),
	distanceMatrix(),
	neighbourhood(),
	ridgeCount(),
	fpclass('U'),
	w(0),
	h(0),
	binaryImage(),
	featurevector(),
	selected_classifier(-1)
{
}

Fingerprint::~Fingerprint()
{
}

Fingerprint::Fingerprint(const Fingerprint& other) :
	id(other.id),
	minutiae(other.minutiae),
	distanceMatrix(other.distanceMatrix),
	neighbourhood(other.neighbourhood),
	ridgeCount(other.ridgeCount),
	fpclass(other.fpclass),
	w(other.w),
	h(other.h),
	binaryImage(other.binaryImage),
	featurevector(other.featurevector),
	selected_classifier(other.selected_classifier)
{
}

Fingerprint& Fingerprint::operator=(const Fingerprint& other)
{
	if (this == &other) return *this;

	id=other.id;
	minutiae=other.minutiae;
	distanceMatrix = other.distanceMatrix;
	neighbourhood = other.neighbourhood;
	ridgeCount= other.ridgeCount;
	fpclass = other.fpclass;
	w=other.w;
	h=other.h;
	binaryImage=other.binaryImage;
	featurevector = other.featurevector;
	selected_classifier = other.selected_classifier;

	return *this;
}


Fingerprint::Fingerprint(const Matrix<int> &xyt, const string &identifier) :
id(identifier),
minutiae(),
distanceMatrix(),
neighbourhood(),
ridgeCount(),
fpclass('U'),
w(0),
h(0),
binaryImage(),
featurevector(),
selected_classifier(-1)
{
	int num_minutiae = xyt.rows();
	
	if(xyt.cols() != 3)
	{
		cerr << "ERROR in Fingerprint constructor: xyt must have 3 columns (for x, y and t)" << endl;
		return;
	}
	
	minutiae.resize(num_minutiae);
	
	for(int i = 0; i < num_minutiae; ++i)
	{
		minutiae[i] = Minutia(i, xyt(i, 0), xyt(i, 1), xyt(i, 2)*2, 100, typeMin::OTH);
	}
}


void Fingerprint::addMinutia(Minutia add)
{
	minutiae.push_back(add);
	distanceMatrix.clear();
	neighbourhood.clear();
	ridgeCount.clear();
}

void Fingerprint::dropMinutia(int index)
{
	if(index<(int)minutiae.size())
		minutiae.erase(minutiae.begin()+index);

	distanceMatrix.clear();
	neighbourhood.clear();
	ridgeCount.clear();
}

void Fingerprint::setMinutia(Minutia newM, int index)
{
	if(index<(int)minutiae.size())
		minutiae[index]=newM;

	distanceMatrix.clear();
	neighbourhood.clear();
	ridgeCount.clear();
}

void Fingerprint::setMinutiae(const vector<Minutia> & val)
{
	minutiae = val;

	distanceMatrix.clear();
	neighbourhood.clear();
	ridgeCount.clear();
}

int Fingerprint::readMinutiaeFile(const string & name, unsigned int quality)
{
	if(name.find(".bir") == name.size()-4 || name.find(".ist") == name.size()-4)
		return read19794file(name, quality);
	else if(name.find(".xyt") == name.size()-4)
		return readXYTFile(name, quality);
	else
		return readXYTFile(name+".xyt", quality);
}

int Fingerprint::readMinutiaeBinFile(const string & name, unsigned int quality)
{
	if(name.find(".bir") == name.size()-4 || name.find(".ist") == name.size()-4)
	{
		int result = read19794file(name, quality);
		string binname = name.substr(0, name.size()-3) + "brw";

		if(result)
			return result;
		else
			return readBinaryFile(binname);
	}
	else
		return readNIGOSfile(name, quality);
}

int Fingerprint::readBinaryFile(const string & name)
{
	ifstream bFile;
	char * memblock;
	int sizeF;

	bFile.open(name.c_str(), fstream::binary);

	if(bFile.is_open()){

		//read the binary fingerprint image
		binaryImage.resize(h,w);
		bFile.seekg (0, ios::end);
		sizeF = bFile.tellg();

		memblock = new char [sizeF];
		bFile.seekg (0, ios::beg);
		bFile.read (memblock, sizeF);
		bFile.close();

		for (int i=0; i<h; i++)
			for (int j=0; j<w; j++)
				binaryImage(i,j)=(memblock[(i*w)+j]==0);

		delete[] memblock;

		return 0;
	}
	else
	{
		cerr << "Could not open binary file " << name << endl;
		bFile.close();
		return -1;
	}
}


int Fingerprint::readNIGOSfile(const string & name, unsigned int quality)
{
	string basename;

	clear();

	if(name.find(".xyt") == name.size()-4)
		basename = name.substr(0, name.size()-4);
	else
		basename = name;

	int result = readMinutaeNIGOSfile(basename, quality);

	if(result)
		return result;
	else
		return readBinaryFile(basename + ".brw");
}

void Fingerprint::clear(){
	minutiae.clear();
	distanceMatrix.clear();
	neighbourhood.clear();
	ridgeCount.clear();
	binaryImage.clear();
}


int Fingerprint::readXYTFile(const string &name, unsigned int quality){

	ifstream xFile;

	string buffer;

	int X,Y;
	int Tx;
	unsigned int Qxyt;
	Minutia newM;

	clear();
	xFile.open(name.c_str());

	if(xFile.is_open()){

		id = name;

// 	cout << "\tReading file " << id << endl;

		// Set standard resolution
		if(w == 0 || h == 0)
		{
			w = 288;
			h = 384;
		}

		getline(xFile,buffer);

		//read the minutiae
		for(int i=0; xFile; i++){

				// Read the angle from XYT file, which has a better resolution (X and Y are the same as in .min)
				// Convert the angle to standard coordinates (as in ISO 19794-2)
				sscanf (buffer.c_str(), "%d %d %d %d", &X, &Y, &Tx, &Qxyt);

				if(Qxyt >= quality || Qxyt == 0)
				{
					newM.setIndex(i);
					newM.setX(X);
					newM.setY(Y);
					newM.setT(Tx*2);
					newM.setQuality(Qxyt);
					newM.setType(OTH);

					addMinutia(newM);
				}

				getline(xFile,buffer);
		}

		xFile.close();

		return 0;
	}
	else
	{
		xFile.close();

		return -1;
	}
}

int Fingerprint::readMinutaeNIGOSfile(const string &name, unsigned int quality){

	ifstream mFile, xFile;

	string buffer;

	int number;
	int X,Y;
	float Q;
	int Tx, Tm;
	unsigned int Qxyt;
	typeMin type;
	char ctype[4];
	Minutia newM;

// 	cout << "\tClearing" << endl;
	clear();

// 	cout << "\tAll clear" << endl;

	mFile.open((name+".min").c_str());
	xFile.open((name+".xyt").c_str());
	if(xFile.is_open() && mFile.is_open()){

		id = name;

 //	cout << "\tReading file " << id << endl;

		//First line contains the resolution of the image
		getline (mFile,buffer);
		sscanf (buffer.c_str(), "Image (w,h) %d %d", &w, &h);
// cout << "w = " << w << ", h = " << h << endl;
		//Second line can be ignored
		getline (mFile,buffer);

		//Third line has the number of minutia
		getline (mFile,buffer);
		sscanf (buffer.c_str(), "%d", &number);

		//Next lines are ommited
		getline (mFile,buffer);
		getline(xFile,buffer);

		//read the minutiae
		for(int i=0; i < number; i++){

				// Read the angle from XYT file, which has a better resolution (X and Y are the same as in .min)
				// Convert the angle to standard coordinates (as in ISO 19794-2)
				sscanf (buffer.c_str(), "%d %d %d %d", &X, &Y, &Tx, &Qxyt);

				//read index, quality and type from mFile
				getline(mFile,buffer);

				if(Qxyt >= quality || Qxyt == 0)
				{

// 	cout << "\tReading minutia " << i << " of " << number << endl;
					// Remove the spaces
					buffer.erase(removeIf(buffer.begin(), buffer.end(), ::isspace), buffer.end());

					buffer = buffer.substr(buffer.find(":")+1);

					// Read the minutiae components
					sscanf (buffer.c_str(), "%d,%d:%d:%f:%3s", &X, &Y, &Tm, &Q, ctype);

					if (!strcmp(ctype, ("RIG")))
							type = RIG;
					else
							type = BIF;

					newM.setIndex(i);
					newM.setX(X);
					newM.setY(Y);
					newM.setT(Tx*2);
					newM.setQuality(Qxyt);
					newM.setType(type);

					addMinutia(newM);
				}

				getline(xFile,buffer);
		}

		mFile.close();
		xFile.close();

		return 0;
	}
	else
	{
		mFile.close();
		xFile.close();

		return -1;
	}
}


int Fingerprint::read19794file(const string & path, unsigned int quality)
{
	File19794 f(path, quality);

	clear();
	id = path;
	w = f.getWidth();
	h = f.getHeight();
	minutiae = f.getMinutiae();
	fpclass = f.getClass();

	return 0;
}


void Fingerprint::computeDistances(){

    int length = minutiae.size();
    int coordx, coordy;

// 		cout << "Length = " << length << endl;

    distanceMatrix.resize(length,length);

    for (int i=0; i<length; ++i) {
	coordx = minutiae[i].getX();
	coordy = minutiae[i].getY();
        distanceMatrix(i,i) = 0.0;
        for (int j=i+1; j<length; ++j) {
            distanceMatrix(i,j)=sqrt((float)square(coordx - minutiae[j].getX())+square(coordy -minutiae[j].getY()));
            distanceMatrix(j,i)=distanceMatrix(i,j);
        }
    }

}

float Fingerprint::getDistance(unsigned int i, unsigned int j) const{

    if (i < minutiae.size() && j < minutiae.size() && !distanceMatrix.empty()) {
        return distanceMatrix.Get(i,j);
    } else {
        return -1;
    }
}

/**
 * Establish an order between pair, using the second component
 * \param i First pair
 * \param j Second pair
 * \return bool order between pairs
 */
bool orderPair (pair <int,float> i,pair <int,float> j) { return (i.second<j.second); }

void Fingerprint::computeNeighbourhood (){

	if (distanceMatrix.size()>0) {

		unsigned int num_minutiae = minutiae.size();
		vector<pair <int,float> > line(num_minutiae-1);

		neighbourhood.clear();
		neighbourhood.resize(num_minutiae,num_minutiae-1);

		for (unsigned int i=0; i<num_minutiae; ++i) {

			for (unsigned int j=0; j<num_minutiae; ++j) {
				if (i>j) {
					line[j].first = j;
					line[j].second = distanceMatrix(i,j);
				}
				else if (i<j) {
					line[j-1].first = j;
					line[j-1].second = distanceMatrix(i,j);
				}
			}
			sort (line.begin(),line.end(),orderPair);

			for (unsigned int j=0; j<line.size(); ++j) {
				neighbourhood[i][j] = line[j].first;
			}
		}
	}
}

int Fingerprint::getNeighbour (unsigned int min, unsigned int i) const
{
	if (min < minutiae.size() && i < minutiae.size()-1)
		return neighbourhood[min][i];
	else
		return -1;
}

void Fingerprint::computeRidgeCount(){

    int i, j, k;
    int x1,y1,x2,y2;
    int aux;
    int width, height;
    float gradient;
    float acc;
    int count;
    bool valley;
    int offset;
    int length = minutiae.size();
    bool horizontal;

		if(length == 0 || binaryImage.empty())
		{
			binaryImage.clear();
			return;
		}

    ridgeCount.resize(length,length);

    for (i=0; i<length; i++) {
        ridgeCount(i,i)=0;
        for (j=i+1; j<length; j++) {
            x1=minutiae[i].getX();
            y1=minutiae[i].getY();
            x2=minutiae[j].getX();
            y2=minutiae[j].getY();

            if (abs(x2-x1) > abs(y2-y1)) {
                horizontal = true;
                if (x1 > x2) {
                    aux = x1;
                    x1 = x2;
                    x2 = aux;
                    aux = y1;
                    y1 = y2;
                    y2 = aux;
                }
                width = x2-x1;
                height = y2-y1;
                gradient = (float)height / (float)width;
            } else {
                horizontal = false;
                if (y1 > y2) {
                    aux = y1;
                    y1 = y2;
                    y2 = aux;
                    aux = x1;
                    x1 = x2;
                    x2 = aux;
                }
                width = x2-x1;
                height = y2-y1;
                gradient = (float)width / (float)height;
            }
            acc = 0.0;
            count = 0;
            valley = false;
            if (horizontal) {
                for (k=0; k<width; k++) {
                    offset = roundInt(acc);
                    if (binaryImage(y1+offset+1,x1+k-1)) { //currently is a ridge
                        if (valley) { //a new ridge is found
                            count++;
                            valley = false;
                        }
                    } else { //currently is a valley
                        if (!valley) { //a new valley is found
                            valley = true;
                        }
                    }
                    acc += gradient;
                }
            } else {
                for (k=0; k<height; k++) {
                    offset = roundInt(acc);

                    if (binaryImage(y1+k+1,x1+offset-1)) { //currently is a ridge
                        if (valley) { //a new ridge is found
                            count++;
                            valley = false;
                        }
                    } else { //currently is a valley
                        if (!valley) { //a new valley is found
                            valley = true;
                        }
                    }
                    acc += gradient;
                }
            }

            ridgeCount(i,j) = max(count-1, 0);
            ridgeCount(j,i) = ridgeCount(i,j);
        }
    }

    // The binary image is not necessary any more
    //binaryImage.clear();
}

int Fingerprint::getRidgeCount(unsigned int i, unsigned int j) const{

    if (i < minutiae.size() && j < minutiae.size() && !ridgeCount.empty()) {
        return ridgeCount.Get(i,j);
    } else {
        return -1;
    }
}


/**
 * Prints a fingerprint
 * \param output Output stream
 * \param F Fingerprint
 */
ostream& operator<<(ostream& output, const Fingerprint& F){

    int i, j;

    output << "Fingerprint #" << F.id << endl;
    output << "minutias: " << endl;

    for(int i=0;i<(int)F.size();i++){
        output << F.getMinutia(i) << endl;
    }

    output << endl;

    if (F.distanceMatrix.size()>0) {
        output << "Distances: " << endl;
        for (i=0; i<(int)F.minutiae.size(); i++) {
            for (j=0; j<(int)F.minutiae.size(); j++) {
                output << setprecision(6) << setw(7) << F.distanceMatrix.Get(i,j) << " ";
            }
            output << endl;
        }
    }

    if (!F.neighbourhood.empty()) {
        output << "Neighbourhood: " << endl;
        for (i=0; i<(int)F.minutiae.size(); i++) {
            for (j=0; j<(int)F.minutiae.size()-1; j++) {
                output << setw(3) << F.neighbourhood[i][j] << " ";
            }
            output << endl;
        }
    }

    if (F.ridgeCount.size()>0) {
        output << "Ridge Count: " << endl;
        for (i=0; i<(int)F.minutiae.size(); i++) {
            for (j=0; j<(int)F.minutiae.size(); j++) {
                output << setw(3) << F.ridgeCount.Get(i,j) << " ";
            }
            output << endl;
        }
    }

    return output;
}


void Fingerprint::setFeatureVector(const std::vector<double> & val, int selected_classifier)
{
	featurevector = val;
	this->selected_classifier = selected_classifier;

}

void Fingerprint::setFeatureVector(const double *val, unsigned int size, int selected_classifier)
{
	featurevector.assign(val, val+size);
	this->selected_classifier = selected_classifier;
}
