/**
 * \file    File19794.h
 * \author  Daniel Peralta <dperalta@decsai.ugr.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file for the File19794 class
 */

#ifndef FILE19794_H
#define FILE19794_H

#include <vector>
#include <fstream>
#include "Minutia.h"
#include "Matrix.h"

/**
 * @class File19794
 *
 * The File19794 class handles the fingerprint minutiae files with the ISO 19794-2 standard format. This class handles both 2005 and 2011 versions of che standard.
 */
class File19794{
	
public:
	
	/** Default constructor */
	File19794();
	
	/**
	 * Primitive constructor
	 * \param file Full path of the minutiae file
	 * \param quality Minimum quality for the read minutiae
	 */
	File19794(const std::string &file, unsigned int quality = 0);
	
	/** Default destructor */
	virtual ~File19794();
	
	/**
	 * Copy constructor
	 * \param other Object to copy from
	 */
	File19794(const File19794& other);
	
	/**
	 * Assignment operator
	 * \param other Object to assign from
	 * \return A reference to this
	 */
	File19794& operator=(const File19794& other);
	
	/**
	 * Get the id
	 * \return The identification of the fingerprint
	 */
	std::string getId() const;
	
	/**
	 * Get the image width
	 * \return The width of the fingerprint image
	 */
	int getWidth() const;

	/**
	 * Get the image height
	 * \return The height of the fingerprint image
	 */
	int getHeight() const;

	/**
	 * Get the number of minutiae
	 * \return The number of minutiae
	 */
	int getNumMinutiae() const;
	
	/**
	 * Access minutiae
	 * \return A copy of the vector of minutiae
	 */
	std::vector<Minutia> getMinutiae() const;
	
	/**
	 * Access the ridge matrix
	 * \return A copy of the ridge matrix
	 */
	Matrix<int> getRidgeCountMatrix() const;
	
	/**
	 * Access to a minutia
	 * \param index Index of the minutia
	 * \return A copy of the minutia selected
	 */
	Minutia getMinutia(int index) const;
	
	/**
	 * Size of the file
	 * \return The number of bytes in the file
	 */
	unsigned int size() const;
	
	/**
	 * Set id
	 * \param val Identification of the file
	 */
	void setId(std::string val);

	/**
	 * Read a data file in 19794-2 format. It requires a single binary (.bir, .ist) file. Both 2005 and 2011 standard versions are supported.
	 * \param name String with the file name (with extension)
	 * \param quality Minimum quality for the read minutiae
	 * \return A negative integer if an error has ocurred. 0 otherwise
	 */
	int readFile(const std::string & name, unsigned int quality = 0);

	/**
	 * Read a data file in 19794-2(2005) format. It requires a single binary (.ist) file.
	 * \param bFile Input binary stream for the file
	 * \param quality Minimum quality for the read minutiae
	 * \return A negative integer if an error has ocurred. 0 otherwise
	 */
	int readFile2005(std::ifstream &bFile, unsigned int quality = 0);

	/**
	 * Read a data file in 19794-2(2011) format. It requires a single binary (.ist) file.
	 * \param bFile Input binary stream for the file
	 * \param quality Minimum quality for the read minutiae
	 * \return A negative integer if an error has ocurred. 0 otherwise
	 */
	int readFile2011(std::ifstream &bFile, unsigned int quality = 0);


	/**
	 * Writes the read 19794 file contents to files in XYT format
	 * \param oroot Base name for the output files. The extensions are added to this name.
	 * \return A negative integer if an error has ocurred. 0 otherwise
	 */
	int writeXYTFiles(const std::string &oroot) const;

	/**
	 * Appends an extended record to the 19794-2 file, that contains the fingerprint class.
	 * \param fpclass Class of the fingerprint
	 * \return A negative integer if an error has ocurred. 0 otherwise
	 */
	int addClass(char fpclass);


	/**
	 * Get the fingerprint class that was read in the file
	 * \return The fingerprint class
	 */
	char getClass() const;

	
protected:


	typedef enum {RESERVED, RIDGECOUNT, COREDELTA, ZONALQUALITY, CLASS, VENDORDEFINED} extended_record_t;

	typedef struct
	{
		extended_record_t type;
		unsigned short int length;
	} ExtDataAreaInfo;

	
	std::string id;                            //!< Path of the original file
	std::vector<Minutia> minutiae;            //!< Set of Minutiae
	Matrix<int> ridgeCount;               //!< Ridge count between minutiae

	// Header fields
	unsigned int file_size;               //!< Size of the file
	unsigned char capture_equipment[2];
	int width;                            //!< Width of the image
	int height;                           //!< Heigth of the image
	unsigned int x_resolution;
	unsigned int y_resolution;
	unsigned int num_finger_views;

	// Finger record fields
	unsigned char finger_position;
	unsigned char view_impression;
	unsigned int  finger_quality;
	unsigned int  num_minutiae;

	// Extended records
	unsigned int ext_data_length;
	std::vector<ExtDataAreaInfo> ext_data_info;
	char fpclass;


	/**
	 * Reads the class Extended Area from the file.
	 * \param memblock Pointer to the beginning of the extended area
	 * \param length Length of the extended area
	 * \return The read class
	 * \see getClassCharacter
	 */
	char readEAClass(char *memblock, unsigned int length) const;

	/**
	 * Returns the type of an extended area according with the type code in the file.
	 * \param byte1 First byte of the code
	 * \param byte2 Second byte of the code
	 * \return The type corresponding to the inserted code
	 */
	extended_record_t getExtendedRecordType(unsigned char byte1, unsigned char byte2) const;
};

inline std::string File19794::getId() const{ return id; }

inline int File19794::getWidth() const { return width; }

inline int File19794::getHeight() const { return height; }

inline int File19794::getNumMinutiae() const { return minutiae.size(); }

inline std::vector<Minutia> File19794::getMinutiae() const { return minutiae; }

inline Minutia File19794::getMinutia(int index) const { return minutiae[index];}

inline Matrix<int> File19794::getRidgeCountMatrix() const { return ridgeCount; }

inline unsigned int File19794::size() const{ return file_size;}

inline void File19794::setId(std::string val) { id = val;}

inline char File19794::getClass() const {return fpclass;}

inline char File19794::readEAClass(char *memblock, unsigned int length) const {return memblock[0];}


#endif
