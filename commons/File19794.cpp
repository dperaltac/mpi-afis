/**
* \file    File19794.cpp
* \author  Daniel Peralta <dperalta@decsai.ugr.es>
* \version 1.0
*
* \section DESCRIPTION
*
* The implementation of the File19794 class
*/

#include <fstream>
#include <string>
#include <cstring>

#include "File19794.h"
#include "Functions.h"

using namespace std;

File19794::File19794()
{
	id = "";
	fpclass = 'U';
}

File19794::File19794(const string &file, unsigned int quality)
{
	fpclass = 'U';
	readFile(file, quality);
}


File19794::~File19794(){}


File19794::File19794(const File19794& other)
{
	id = other.id;
	minutiae = other.minutiae;
	ridgeCount = other.ridgeCount;
	file_size = other.file_size;
	capture_equipment[0] = other.capture_equipment[0];
	capture_equipment[1] = other.capture_equipment[1];
	width = other.width;
	height = other.height;
	x_resolution = other.x_resolution;
	y_resolution = other.y_resolution;
	num_finger_views = other.num_finger_views;
	finger_position = other.finger_position;
	view_impression = other.view_impression;
	finger_quality = other.finger_quality;
	num_minutiae = other.num_minutiae;
	ext_data_length = other.ext_data_length;
	ext_data_info = other.ext_data_info;
	fpclass = other.fpclass;
}


File19794& File19794::operator=(const File19794& other)
{
	if(&other != this)
	{
		id = other.id;
		minutiae = other.minutiae;
		ridgeCount = other.ridgeCount;
		file_size = other.file_size;
		capture_equipment[0] = other.capture_equipment[0];
		capture_equipment[1] = other.capture_equipment[1];
		width = other.width;
		height = other.height;
		x_resolution = other.x_resolution;
		y_resolution = other.y_resolution;
		num_finger_views = other.num_finger_views;
		finger_position = other.finger_position;
		view_impression = other.view_impression;
		finger_quality = other.finger_quality;
		num_minutiae = other.num_minutiae;
		ext_data_length = other.ext_data_length;
		ext_data_info = other.ext_data_info;
		fpclass = other.fpclass;
	}

	return *this;
}

int File19794::readFile(const string & path, unsigned int quality)
{
	char block[8];
	ifstream bFile(path.c_str(), fstream::binary);

	id = path;

	/********** READ HEADER ***********/

	// Read the block
	bFile.read (block, 8);
	block[3] = '\0';
	block[7] = '\0';

	if(!bFile)
	{
		cerr << "ERROR: could not open file " << path << endl;
		return -1;
	}

	// Check standard version
	if (strcmp(block, "FMR") || (strcmp(block+4, "030") && strcmp(block+4, " 20")))
	{
		cerr << "Error while reading file " << path << " in 19794-2 format: expected \"FMR 030\" or \"FMR  20\", read \"" << block << " " << block+4 << "\"" << endl;
		return -2;
	}
	else if (strcmp(block+4, " 20") == 0)
		return readFile2005(bFile, quality);
	else if (strcmp(block+4, "030") == 0)
		return readFile2011(bFile, quality);

	return -1;
}


int File19794::readFile2005(ifstream &bFile, unsigned int quality)
{
	char *memblock, block[12];
	unsigned int sizeF;
	typeMin type_min;
	unsigned int edl_sum;
	ExtDataAreaInfo edai;

	// Read the actual file size
	bFile.seekg (0, ios::end);
	sizeF = bFile.tellg();
	
	bFile.seekg (8, ios::beg);
	bFile.read((char*)&file_size, 4);

	// Convert big-endian to little-endian
	file_size = longIntegerSwap(file_size);

	if(file_size != sizeF)
	{
		cerr << "Error while reading file " << id << " in 19794-2 format: file states " << file_size << " bytes, but contains " << sizeF << " instead" << endl;
		return -3;
	}

	// Read the rest of the record block
// 	bFile.seekg (12, ios::beg);
	bFile.read (block, 12);

	capture_equipment[0] = block[0];
	capture_equipment[1] = block[1];

	width        = getShortInt(block[2], block[3]);
	height       = getShortInt(block[4], block[5]);
	x_resolution = getShortInt(block[6], block[7]);
	y_resolution = getShortInt(block[8], block[9]);
	num_finger_views = static_cast<unsigned char>(block[10]);

	if(block[11] != '\0')
		cerr << "WARNING while reading file " << id << " in 19794-2 format: the block reserved byte is not set to zero." << endl;



	/********** READ FINGER RECORD ***********/
// 	bFile.seekg (24, ios::beg);
	bFile.read(block, 4);
	finger_position = block[0];
	view_impression = block[1];
	finger_quality  = static_cast<unsigned char>(block[2]);
	num_minutiae    = static_cast<unsigned char>(block[3]);


	/********** READ MINUTIAE ***********/

	memblock = new char[6*num_minutiae + 2];

// 	bFile.seekg (28, ios::beg);
	bFile.read (memblock, 6*num_minutiae + 2);	// Read two bytes more to get the size of the extended area

// 	minutiae.resize(num_minutiae);
	minutiae.clear();

// 	cout << "MINUTIAE: " << num_minutiae << endl;

	for (unsigned int i = 0, desp = 0; i < num_minutiae; ++i, desp+=6)
	{
		unsigned int coordX, coordY;
		unsigned char min_type = (memblock[desp] & 0xC0) >> 6;
		float angle = 0;
		unsigned int Q;

		Minutia newmin;

		Q = (unsigned char)memblock[5 + desp];

		if(Q >= quality || Q == 0)
		{
			if(min_type == 0x1)
				type_min = RIG;
			else if(min_type == 0x2)
				type_min = BIF;
			else
			{
// 				fprintf(stderr, "Fingerprint file %s: skipping minutia number %d of type OTHER (memblock: %x   min_type: %x)\n", id.c_str(), i, memblock[desp], min_type);
// 				continue;
				type_min = OTH;
			}

			coordX = ((unsigned int)(memblock[    desp] & 0x3F) << 8) + (unsigned char)memblock[1 + desp];
			coordY = ((unsigned int)(memblock[2 + desp] & 0x3F) << 8) + (unsigned char)memblock[3 + desp];

			angle = (unsigned char)memblock[4 + desp] * 1.40625;

			newmin.setIndex  (i);
			newmin.setX      (coordX);
			newmin.setY      (coordY);
			newmin.setType   (type_min);
			newmin.setT      (angle);
			newmin.setQuality(Q);

			minutiae.push_back(newmin);

// 			cout << coordX << " " << coordY << " " << roundInt(angle) << " " << Q << endl;
		}
	}

	ext_data_length = getShortInt(memblock[6*num_minutiae], memblock[6*num_minutiae+1]);

	delete [] memblock;

	if(ext_data_length > 0)
	{
		memblock = new char[ext_data_length];

		/********** READ EXTENDED DATA INFO ***********/
	// 	bFile.seekg (28 + 6*num_minutiae + 2, ios::beg);
		bFile.read (memblock, ext_data_length);

		// This variable counts the number of read bytes in the extended area
		edl_sum = 0;
		ext_data_info.clear();

		while(edl_sum < ext_data_length)
		{
			edai.type =   getExtendedRecordType(memblock[ext_data_info.size()*4  ], memblock[ext_data_info.size()*4+1]);
			edai.length = getShortInt          (memblock[ext_data_info.size()*4+2], memblock[ext_data_info.size()*4+3]);
			ext_data_info.push_back(edai);

			edl_sum += edai.length + 4;
		}

		edl_sum = ext_data_info.size()*4;

		for(unsigned int i = 0; i < ext_data_info.size(); i++)
		{
			switch(ext_data_info[i].type)
			{
				case RESERVED:
					cerr << "WARNING: the fingerprint file " << id << " has an extended area with a reserved code" << endl;
				break;
				case RIDGECOUNT:
				break;
				case COREDELTA:
				break;
				case ZONALQUALITY:
				break;
				case CLASS:
					fpclass = readEAClass(memblock+edl_sum, ext_data_info[i].length);
				break;
				case VENDORDEFINED:
					cerr << "WARNING: the fingerprint file " << id << " has an extended area with a vendor-defined code" << endl;
				break;
			}

			edl_sum += ext_data_info[i].length;
		}

		delete [] memblock;
	}

	bFile.close();

	return 0;
}


int File19794::readFile2011(ifstream &bFile, unsigned int quality)
{
	char *memblock;
	unsigned int sizeF;
	typeMin type_min;
	unsigned int edl_sum;
	ExtDataAreaInfo edai;
	short int finger_representations;
	unsigned char certification_block, quality_length, minutiae_field_length, certification_length;
	int representation_length, shift;

	// Read the actual file size
	bFile.seekg (0, ios::end);
	sizeF = bFile.tellg();

	bFile.seekg (8, ios::beg);
	bFile.read((char*)&file_size, 4);

	// Convert big-endian to little-endian
	file_size = longIntegerSwap(file_size);

	if(file_size != sizeF)
	{
		cerr << "Error while reading file " << id << " in 19794-2 format: file states " << file_size << " bytes, but contains " << sizeF << " instead" << endl;
		return -3;
	}

	// Read the number of finger representations
// 	bFile.seekg (12, ios::beg);
	bFile.read((char*)&finger_representations, 2);
	finger_representations = shortIntegerSwap(finger_representations);

	if(finger_representations != 1)
	{
		cerr << "Error while reading file " << id << " in 19794-2 format: the file contains" << finger_representations << " finger representations, but the maximum allowed by this program is 1" << endl;
		return -4;
	}

// 	bFile.seekg (14, ios::beg);
	bFile.read ((char*)&certification_block, 1);

// 	bFile.seekg (15, ios::beg);
	bFile.read((char*)&representation_length, 4);
	representation_length = longIntegerSwap(representation_length);


	memblock = new char[representation_length];

	/********** READ FINGERPRINT RECORD ***********/
// 	bFile.seekg (19, ios::beg);
	bFile.read (memblock, representation_length);

	quality_length = static_cast<unsigned char>(memblock[14]);

	if(certification_block == 0x1)
	{
		certification_length = static_cast<unsigned char>(memblock[15+5*quality_length]);
		shift = 15+5*quality_length+1+3*certification_length;
	}
	else
		shift = 15+5*quality_length;

	finger_position = static_cast<unsigned char>(memblock[shift]);

	// Skip representation number, which is 0

	x_resolution = getShortInt(memblock[shift+2], memblock[shift+3]);
	y_resolution = getShortInt(memblock[shift+4], memblock[shift+5]);
	view_impression = static_cast<unsigned char>(memblock[shift+6]);

	width = getShortInt(memblock[shift+7], memblock[shift+8]);
	height = getShortInt(memblock[shift+9], memblock[shift+10]);

	minutiae_field_length = static_cast<unsigned char>(memblock[shift+11] >> 4);

	/********** READ FINGER RECORD ***********/
	num_minutiae    = static_cast<unsigned char>(memblock[shift+12]);
	

	/********** READ MINUTIAE ***********/

// 	minutiae.resize(num_minutiae);
	minutiae.clear();

	cout << "MINUTIAE: " << num_minutiae << endl;

	for (unsigned int i = 0, desp = shift+13; i < num_minutiae; ++i, desp+=minutiae_field_length)
	{
		unsigned int coordX, coordY;
		unsigned char min_type = (memblock[desp] & 0xC0) >> 6;
		float angle = 0;
		unsigned int Q;

		Minutia newmin;

		if(minutiae_field_length == 6)
			Q = (unsigned char)memblock[5 + desp];
		else
			Q = 0;

		if(Q >= quality || Q == 0)
		{
			if(min_type == 0x1)
				type_min = RIG;
			else if(min_type == 0x2)
				type_min = BIF;
			else
			{
				fprintf(stderr, "Fingerprint file %s: skipping minutia number %d of type OTHER (memblock: %x   min_type: %x)\n", id.c_str(), i, memblock[desp], min_type);
				continue;
			}

			coordX = ((unsigned int)(memblock[    desp] & 0x3F) << 8) + (unsigned char)memblock[1 + desp];
			coordY = ((unsigned int)(memblock[2 + desp] & 0x3F) << 8) + (unsigned char)memblock[3 + desp];

			angle = (unsigned char)memblock[4 + desp] * 1.40625;

			newmin.setIndex  (i);
			newmin.setX      (coordX);
			newmin.setY      (coordY);
			newmin.setType   (type_min);
			newmin.setT      (angle);
			newmin.setQuality(Q);

			minutiae.push_back(newmin);

// 			cout << coordX << " " << coordY << " " << roundInt(angle) << " " << Q << endl;
		}
	}

	ext_data_length = getShortInt(memblock[shift+13+minutiae_field_length*num_minutiae], memblock[shift+13+minutiae_field_length*num_minutiae+1]);

	delete [] memblock;

	if(ext_data_length > 0)
	{
		memblock = new char[ext_data_length];

		/********** READ EXTENDED DATA INFO ***********/
	// 	bFile.seekg (28 + 6*num_minutiae + 2, ios::beg);
		bFile.read (memblock, ext_data_length);

		// This variable counts the number of read bytes in the extended area
		edl_sum = 0;
		ext_data_info.clear();

		while(edl_sum < ext_data_length)
		{
			edai.type =   getExtendedRecordType(memblock[ext_data_info.size()*4  ], memblock[ext_data_info.size()*4+1]);
			edai.length = getShortInt          (memblock[ext_data_info.size()*4+2], memblock[ext_data_info.size()*4+3]);
			ext_data_info.push_back(edai);

			edl_sum += edai.length + 4;
		}

		edl_sum = ext_data_info.size()*4;

		for(unsigned int i = 0; i < ext_data_info.size(); i++)
		{
			switch(ext_data_info[i].type)
			{
				case RESERVED:
					cerr << "WARNING: the fingerprint file " << id << " has an extended area with a reserved code" << endl;
				break;
				case RIDGECOUNT:
				break;
				case COREDELTA:
				break;
				case ZONALQUALITY:
				break;
				case CLASS:
					fpclass = readEAClass(memblock+edl_sum, ext_data_info[i].length);
				break;
				case VENDORDEFINED:
					cerr << "WARNING: the fingerprint file " << id << " has an extended area with a vendor-defined code" << endl;
				break;
			}

			edl_sum += ext_data_info[i].length;
		}

		delete [] memblock;
	}

	bFile.close();

	return 0;
}


int File19794::writeXYTFiles(const string &oroot) const
{
	ofstream xytfile((oroot + ".xyt").c_str());

	if(!xytfile)
	{
		cerr << "Error when writing file " << oroot << ".xyt" << endl;
		return -1;
	}

	for(vector<Minutia>::const_iterator i = minutiae.begin(); i != minutiae.end(); ++i)
	{
		xytfile << i->getX() << " " << i->getY() << " " << roundInt(i->getT()) << " " << i->getQuality() << endl;
	}

	return 0;
}


int File19794::addClass(char newfpclass)
{
	if(id == "")
	{
		cerr << "ERROR when adding a class: no 19794-2 file was previously opened" << endl;
		return -1;
	}

	if(fpclass != 'U')
	{
		cerr << "ERROR when adding a class to " << id << ": this file already has the class \'" << fpclass << "\'" << endl;
		return -1;
	}

	fstream bFile(id.c_str(), fstream::binary | ios::in | ios::out);
	unsigned int new_size = file_size + 1 + 4;
	unsigned int new_size_be = longIntegerSwap(new_size);
	unsigned short int new_ea_length_be = shortIntegerSwap(ext_data_length + 1 + 4);
	unsigned int old_ea_size = file_size - (28 + 6*num_minutiae + 2 + ext_data_info.size()*4);
	char *old_ea = new char[old_ea_size];
	char class_ea_type[2] = {0x01, 0x01};
	unsigned short int class_ea_size_be = shortIntegerSwap(1);

	if(!bFile.is_open())
	{
		cerr << "ERROR: could not open file " << id << endl;
		return -1;
	}

	fpclass = newfpclass;

	// Increment the file size
	bFile.seekg(8, ios::beg);
	bFile.write((char *)&new_size_be, 4);

	// Increment the extended area length
	bFile.seekg(28 + 6*num_minutiae, ios::beg);
	bFile.write((char *)&new_ea_length_be, 2);

	// Store the old extended areas
	bFile.seekg(file_size - old_ea_size, ios::beg);
	bFile.read (old_ea, old_ea_size);

	// Add extended area info field
	bFile.seekg(file_size - old_ea_size, ios::beg);
	bFile.write(class_ea_type, 2);
	bFile.write((char *)&class_ea_size_be, 2);

	// Restore the previous extended area contents
	bFile.write(old_ea, old_ea_size);
	
	// Append the class extended area block
	bFile.write(&fpclass, 1);

	file_size = new_size;
	
	return 0;
}


File19794::extended_record_t File19794::getExtendedRecordType(unsigned char byte1, unsigned char byte2) const
{
	if(byte1 == 0x01 && byte2 == 0x01)
		return CLASS;
	
	if(byte1 == 0x00)
	{
		if(byte2 == 0x00)
			return RESERVED;
		else if(byte2 == 0x01)
			return RIDGECOUNT;
		else if(byte2 == 0x02)
			return COREDELTA;
		else if(byte2 == 0x03)
			return ZONALQUALITY;
		else
			return RESERVED;
	}
	else
	{
		if(byte2 == 0x00)
			return RESERVED;
		else
			return VENDORDEFINED;
	}
}

