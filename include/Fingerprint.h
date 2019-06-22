/**
 * \file    Fingerprint.h
 * \author  Joaquin Derrac <jderrac@decsai.ugr.es>
 * \author  Salvador Garcia <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file for the fingerprint class
 */

#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include <vector>
#include <cmath>
#include <iostream>

#include "Score.h"
#include "Matrix.h"
#include "Minutia.h"

/**
 * @class Fingerprint
 *
 * The Fingerprint class represents a Fingerprint as a vector of minutaes, three matrixes of neighbors, distances and ridges between minutae,
 * and a bit representation of its binary image.
 */
class Fingerprint
{
	public:

		/** Default constructor */
		Fingerprint();

		/** Default destructor */
		virtual ~Fingerprint();
		
		/**
		 * Copy constructor
		 * \param other Object to copy from
		 */
		Fingerprint(const Fingerprint& other);
		
		/** Constructor to load the minutiae of a fingerprint
		 * \param xyt Matrix containing 1 row per minutia, and 3 columns (X, Y and T). The quality is set to 100, and the type to "other".
		 * \param identifier Optional name for the fingerprint
		 */
		Fingerprint(const Matrix<int> &xyt, const std::string &identifier = "");

		/**
			* Assignment operator
			* \param other Object to assign from
			* \return A reference to this
			*/
		Fingerprint& operator=(const Fingerprint& other);

		/** Parameter configuration for compatibility with the command line
			* \return 0 if there is no error, a different number otherwise
			*/
		static int configureAlgorithm(int argc, char *argv[]);
		
		/**
		* Returns the valid arguments for getopt
		* \returns String of valid arguments for getopt
		*/
		static std::string getValidArgs();

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
			* Access minutiae
			* \return A copy of the vector of minutiae
			*/
		std::vector<Minutia> getMinutiae() const;

		/**
			* Access the feature vector
			* \return A copy of the feature vector
			*/
		std::vector<double> getFeatureVector() const;

		/**
			* Access neighbourhood
			* \return A copy of the vector of neighbours
			*/
		Matrix<int> getNeighbourhood() const;

		/**
			* Access the distance matrix
			* \return A copy of the distance matrix
			*/
		Matrix<float> getDistanceMatrix() const;

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
		Minutia getMinutia(unsigned int index) const;

		/**
			* Access to a single feature
			* \param index Index of the feature in the feature vector
			* \return The value of the feature
			*/
		double getFeature(unsigned int index) const;

		/**
			* Get the distance between two minutiae
			* \param i First minutia
			* \param j Second minutia
			* \return The distance between the two minutiae referred as indexes. -1 if the distance
			* is not computed or the indexes are out of range.
			*/
			float getDistance(unsigned int i, unsigned int j) const;

		/**
			* Get the i-th neighbour of a minutia
			* \param min Minutia index
			* \param i Neighbour index
			* \return The index of the i-th neighbour of the minutia min. -1 if this neighbour is not computed
			* or is out of range.
			*/
			int getNeighbour(unsigned int min, unsigned int i) const;

		/**
			* Get the ridge count between two minutiae
			* \param i First minutia
			* \param j Second minutia
			* \return The ridge count between the two minutiae . -1 if the ridge count is not computed or the
			* indexes are out of range.
			*/
			int getRidgeCount(unsigned int i, unsigned int j) const;

		/**
			* Get the orientation in a certain pixel of the image
			* \param x coordinate of the pixel
			* \param y coordinate of the pixel
			* \return The orientation from 0 (north direction) to 15, in steps of 11.25 degrees.
			*/
			int getOrientation(unsigned int x, unsigned int y);

		/**
			* Size of the fingerprint
			* \return The number of minutiae in the fingerprint
			*/
			unsigned int size() const;

		/**
			* Set id
			* \param val Identification of the fingerprint
			*/
		void setId(std::string val);

		/**
			* Set the fingerprint class
			* \param fpc New class of the fingerprint
			*/
		void setClass(char fpc);

		/**
			* Get the fingerprint class
			* \return The fingerprint class
			*/
		char getClass() const;

		/**
			* Set minutiae. Cancels the distance and neighbourhood information.
			* \param val Vector of minutiae
			*/
		void setMinutiae(const std::vector<Minutia> & val);

		/**
			* Set a minutia. Cancels the distance and neighbourhood information.
			* \param newM Minutia to set
			* \param index Position to insert
			*/
		void setMinutia(Minutia newM, int index);

		/**
			* Set the feature vector.
			* \param val Feature vector
			* \param selected_classifier Index of the classifier that corresponds to the feature vector \p val
			*/
		void setFeatureVector(const std::vector<double> & val, int selected_classifier = 0);

		/**
			* Set the feature vector.
			* \param val Feature vector
			* \param size Size of the feature vector
			* \param selected_classifier Index of the classifier that corresponds to the feature vector \p val
			*/
		void setFeatureVector(const double *val, unsigned int size, int selected_classifier = 0);

		/**
			* Set a feature.
			* \param newM New value of the feature
			* \param index Position of the feature in the feature vector
			*/
		void setFeature(double val, int index);

		/**
			* Get the classifier used with the current feature vector.
			* \return Index of the selected classifier, -1 if none is set.
			*/
		int getSelectedClassifier() const;

		/**
			* Add a minutia to the end of the vector of minutiae. Cancels the distance and neighbourhood information.
			* \param add Minutia to add
			*/
		void addMinutia(Minutia add);

		/**
			* Remove the minutia in the position specified. Cancels the distance and neighbourhood information.
			* \param index Position of the minutia
			*/
		void dropMinutia(int index);

		/**
			* Read a data file in the necessary format according to the Fingerprint subclass.
			* \param name String with the file base name (with no extension, the function looks for XYT+MIN files)
			* \param quality Minimum quality for the read minutiae
			* \return -1 if an error has ocurred. 0 otherwise
			*/
		virtual int readFile(const std::string & name, unsigned int quality = 0);

		/**
			* Read a data file in NIGOS or 19794-2 format plus the binary image (.brw). The used format depends on the extension in @a name.
			* \param name String with the file base name (with no extension, the function looks for XYT+MIN files)
			* \param quality Minimum quality for the read minutiae
			* \return -1 if an error has ocurred. 0 otherwise
			*/
		int readMinutiaeBinFile(const std::string & name, unsigned int quality = 0);
		
		/**
			* Read a binary file, that contains the information about the fingerprint ridges and valleys. It requires the .brw file.
			* \param name String with the file full name
			* \return -1 if an error has ocurred. 0 otherwise
			*/
		int readBinaryFile(const std::string & name);

		/**
			* Read a data file in NIGOS or 19794-2 format. The used format depends on the extension in @a name.
			* \param name String with the file base name (with no extension, the function looks for XYT+MIN files)
			* \param quality Minimum quality for the read minutiae
			* \return -1 if an error has ocurred. 0 otherwise
			*/
		int readMinutiaeFile(const std::string & name, unsigned int quality = 0);

		/**
			* Read a data file in NIGOS format. It requires the .xyt, .brw and .min files.
			* \param name String with the file base name (without extensions)
			* \param quality Minimum quality for the read minutiae
			* \return -1 if an error has ocurred. 0 otherwise
			*/
		int readNIGOSfile(const std::string & name, unsigned int quality = 0);

		/**
			* Read a data file in XYT format. It requires the .xyt file.
			* \param name String with the file name (with the .xyt extention)
			* \param quality Minimum quality for the read minutiae
			* \return -1 if an error has ocurred. 0 otherwise
			*/
		int readXYTFile(const std::string & name, unsigned int quality = 0);

		/**
			* Read a data file in NIGOS format. It requires the .xyt and .min file.
			* \param name String with the file base name (without extensions)
			* \param quality Minimum quality for the read minutiae
			* \return -1 if an error has ocurred. 0 otherwise
			*/
		int readMinutaeNIGOSfile(const std::string & name, unsigned int quality = 0);
		
		/**
			* Read a data file in 19794-2 format. It requires a single binary (.bir, .ist) file. Both 2005 and 2011 standard versions are supported.
			* \param name String with the file name (with extension)
			* \return A negative integer if an error has ocurred. 0 otherwise
			*/
		int read19794file(const std::string & name, unsigned int quality = 0);
				
		/**
			* Compute the distance matrix for all minutiae
			*/
		void computeDistances();

		/**
			* Compute the neighbourhood of all minutiae. Requires the distance matrix to be computed.
			*/
		void computeNeighbourhood();

		/**
			* Compute the ridge count between minutiae. Requires the distance matrix to be computed.
			*/
		void computeRidgeCount();

		/**
			* Removes all the fingerprint contents.
			*/
		void clear();

		/** Calls all the methods that need to be called before the matching
			*/
		virtual void initialize() = 0;


  /**
   * Free the memory used for storing the binary image
   */
//  void freeBinaryImage();


		/**
		* Match the Fingerprint with a second one.
		* \param f A second Fingerprint
		* \return A float value estimating the cost of the match
		*/
		virtual float match(const Fingerprint &f) const 
		{
			std::cerr << "WARNING: Fingerprint::match is not implemented and should not be called" << std::endl;
			return -1.0;
		}

		/**
		* Determines which score is better
		* \param score1 A score obtained from the method \a match
		* \param score2 A score obtained from the method \a match
		* \return true iv \a score1 is better than \a score2, false otherwise
		*/
		static bool better(float score1, float score2);

		/**
		* Determines which score is better
		* \param score1 A score obtained from the method \a match
		* \param score2 A score obtained from the method \a match
		* \return true iv \a score1 is better than \a score2, false otherwise
		*/
		static bool betterOrEqual(float score1, float score2);

		/**
		* Determines which score is better
		* \param score1 A score obtained from the method \a match
		* \param score2 A score obtained from the method \a match
		* \return true iv \a score1 is better than \a score2, false otherwise
		*/
		static bool better(const Score &score1, const Score &score2);

		/**
		* Determines which score is better
		* \param score1 A score obtained from the method \a match
		* \param score2 A score obtained from the method \a match
		* \return true iv \a score1 is better than \a score2, false otherwise
		*/
		static bool betterOrEqual(const Score &score1, const Score &score2);

		//friend operators
		friend std::ostream& operator<<(std::ostream& output, const Fingerprint& F);

		
	protected:
		std::string id; //!< Identidication of the fingerprint
		std::vector <Minutia> minutiae; //!< Set of Minutiae
		Matrix<float> distanceMatrix; //!< Distances between minutiae
		Matrix<int> neighbourhood; //!< Neighbourhood of the minutiae
		Matrix<int> ridgeCount; //!< Ridge count between minutiae
		char fpclass; //!< Fingerprint class

		int w; //!< Width of the image
		int h; //!< Heigth of the image
		Matrix<bool> binaryImage; //!< Binary image of the fingerprint
		
		std::vector<double> featurevector; //!< Feature vector for classification
		int selected_classifier; //!< Index of the classifier corresponding to the feature vector
		
		static const std::string valid_args; //!< Valid command line arguments
};

inline std::string Fingerprint::getId() const{ return id; }

inline int Fingerprint::getWidth() const { return w; }

inline int Fingerprint::getHeight() const { return h; }

inline std::vector<Minutia> Fingerprint::getMinutiae() const { return minutiae; }

inline std::vector<double> Fingerprint::getFeatureVector() const { return featurevector; }

inline Minutia Fingerprint::getMinutia(unsigned int index) const { return minutiae[index];}

inline double Fingerprint::getFeature(unsigned int index) const { return featurevector[index];}

inline Matrix<int> Fingerprint::getNeighbourhood() const { return neighbourhood; }

inline Matrix<float> Fingerprint::getDistanceMatrix() const { return distanceMatrix; }

inline Matrix<int> Fingerprint::getRidgeCountMatrix() const { return ridgeCount; }

inline unsigned int Fingerprint::size() const{ return minutiae.size();}

inline void Fingerprint::setId(std::string val) { id = val;}

inline char Fingerprint::getClass() const { return fpclass;}

inline void Fingerprint::setClass(char fpc) { fpclass = fpc;}

inline bool Fingerprint::better(float score1, float score2) {return score1 > score2;}

inline bool Fingerprint::betterOrEqual(float score1, float score2) {return score1 >= score2;}

inline bool Fingerprint::better(const Score &score1, const Score &score2) {return better(score1.getScore(), score2.getScore());}

inline bool Fingerprint::betterOrEqual(const Score &score1, const Score &score2) {return betterOrEqual(score1.getScore(), score2.getScore());}

inline std::string Fingerprint::getValidArgs() { return valid_args; }

inline int Fingerprint::configureAlgorithm(int argc, char *argv[]) { return 0; }

inline int Fingerprint::readFile(const std::string & file, unsigned int quality)
{
	id = file;
	return readMinutiaeFile(file, quality);
}

inline void Fingerprint::setFeature(double val, int index) { featurevector[index] = val; }

inline int Fingerprint::getSelectedClassifier() const {return selected_classifier;}

#endif
