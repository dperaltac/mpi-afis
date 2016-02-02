/**
 * \file    Matrix.h
 * \author  Joaquin Derrac <jderrac@decsai.ugr.es>
 * \version 2.0
 *
 * \section DESCRIPTION
 *
 * A template for building generic matrixes
 */

#ifndef MATRIX_H
#define MATRIX_H

#include<limits>
#include<iostream>
#include<list>
#include<algorithm>
#include <sstream>
#include <fstream>

template<typename T> class Matrix;
template<typename T> std::ostream& operator<<(std::ostream& out, const Matrix<T>& M);

template <typename T>

/**
 * @class Matrix
 *
 * The Matrix class allows to define generic matrixes, able to store any of the primitive types.
 * It supports also some additional functions, such as filling or printing (<<).
 */
class Matrix{

public:

    /** Default Constructor */
    Matrix();

    /**
     * Constructor
     * \param rows Rows of the matrix
     * \param cols Columns of the matrix
     */
    Matrix(int rows, int cols);

    /** Default destructor */
    ~Matrix();

    /**
     * Copy constructor
     * \param o Object to copy from
     */
    Matrix(const Matrix& o);

    /**
     * Assignment operator
     * \param o Object to assign from
     * \return A reference to this object
     */
    Matrix & operator= (const Matrix &o);

    /**
     * Erase all the contents of the matrix
     */
    void clear();

    /**
     * Changes the size of the matrix. Data is kept if dimensions are not reduced.
     * \param rows New number of rows of the matrix
     * \param columns New number of columns of the matrix
     */
    void resize(int rows, int columns);

    /**
     * Access operator
     * \param i Row index
     * \param j Column index
     * \return A reference the value stored in (i,j)
     */
    const T& operator() (int i, int j) const;

    /**
     * Get vector operator
     * \param i Row index
     * \return A pointer to the i-th row
     */
    const T* operator[] (int i) const;

    /**
     * Access operator
     * \param i Row index
     * \param j Column index
     * \return A reference the value stored in (i,j)
     */
    T& operator() (int i, int j);

    /**
     * Get vector operator
     * \param i Row index
     * \return A pointer to the i-th row
     */
    T* operator[] (int i);


    /**
     * Get a value of the matrix
     * \param i Row index
     * \param j Column index
     * \return A copy of the value stored in (i,j)
     */
    T Get(int i, int j) const;

    /**
     * Get the number of rows of the matrix
     * \return Number of rows of the matrix
     */
    int rows() const;

    /**
     * Get the number of columns of the matrix
     * \return Number of columns of the matrix
     */
    int cols() const;

    /**
     * Get the number of elements of the matrix (rows*columns)
     * \return Number of elements of the matrix
     */
    int size() const;
    
    /**
     * Get if the matrix is empty or not
     * \return true if the matrix if empty, false otherwise
     */
    bool empty() const;

    /**
     * Get the minimum dimension of the matrix
     * \return the minimum dimension of the matrix
     */
    int minsize(void) {
        return ((nRows < nCols) ? nRows : nCols);
    }

    /**
     * Fills the matrix with a given value
     * \param value Value used for filling the matrix
     */
    void fill(T value);

    /**
     * Normalizes the matrix, keeping the maximum value for each row and column
     */
    void normalizeMaxMatrix();

    T getNormalizedSum() const;

    //friend operators
    friend std::ostream& operator<< <>(std::ostream& out, const Matrix<T>& M);

    /**
     * Return the pointer of the matrix
     */
    T** getPointer();
		
		/**
		 * Reads a matrix from a csv
		 * \param filename Name of the csv file
		 * \param removelast Determines if the last column of the file should be read or not
		 */
		void readFromCSV(const std::string &filename, bool removelast = false);
    
private:
    int nRows; ///< Number of rows
    int nCols; ///< Number of columns
    T ** contents; ///< Matrix contents
};

template<typename T>
Matrix<T>::Matrix() : nRows(0), nCols(0), contents(0) {}

template <typename T>
Matrix<T>::Matrix(int row, int col) : nRows(row), nCols(col), contents(0) {

	if (row == 0 || col == 0)
	{
		nRows = 0;
		nCols = 0;
		contents = 0;
	}
	else
	{
		contents = new T* [nRows];
		contents[0]= new T [nRows*nCols];

		for (int i=1;i<nRows;i++)
			contents[i] = contents[i-1]+nCols;
	}
}

template <typename T>
Matrix<T>::~Matrix(){

    clear();

}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& o) : nRows(o.nRows), nCols(o.nCols), contents(0) {
	
	if (o.contents == 0) {
		nRows = 0;
		nCols = 0;
		contents = 0;
		return;
	}

	contents = new T* [nRows];
	contents[0]= new T [nRows*nCols];

	for (int i=1;i<nRows;++i)
		contents[i] = contents[i-1]+nCols;

	for(int i=0; i<nRows*nCols; i++)
		contents[0][i]=o.contents[0][i];
}

template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &o) {

	if (this == &o)
		return *this;
	else if (o.contents == 0)
	{
		clear();
		return *this;
	}

	const int MxN = o.nRows * o.nCols;

	if (MxN != nRows * nCols)
	{
		clear();
		contents = new T* [o.nRows];
		contents[0]= new T [MxN];

		for (int i=1; i<o.nRows; i++)
			contents[i] = contents[i-1]+o.nCols;
	}

	nRows = o.nRows;
	nCols = o.nCols;
	
	for(int i=0; i<MxN; i++)
		contents[0][i]=o.contents[0][i];

	return *this;
}

template <typename T>
void Matrix<T>::clear() {

	if ( contents != 0 ){
		delete [] contents[0];
		delete [] contents;
		contents = 0;
	}    

	nRows = 0;
	nCols = 0;
}

template <typename T>
void Matrix<T>::resize(int rows, int columns) {

	T ** new_matrix;

	if(rows == 0 || columns == 0)
	{
		clear();
		return;
	}
	
	new_matrix = new T* [rows]; // rows
	new_matrix[0]= new T [rows*columns];

	for (int i=1;i<rows;i++){
		new_matrix[i] = new_matrix[i-1]+columns;
	}

	// copy data from saved pointer to new arrays
	int minrows = std::min<int>(rows, nRows);
	int mincols = std::min<int>(columns, nCols);

	for ( int x = 0 ; x < minrows ; x++ ){
		for ( int y = 0 ; y < mincols ; y++ ){
			new_matrix[x][y] = contents[x][y];
		}
	}
	
	clear();

	contents = new_matrix;

	nRows = rows;
	nCols = columns;

}

template <typename T>
inline T* Matrix<T>::operator[] (int i) {return contents[i];}

template <typename T>
inline T& Matrix<T>::operator() (int i, int j) {return contents[i][j];}

template <typename T>
inline const T* Matrix<T>::operator[] (int i) const {return contents[i];}

template <typename T>
inline const T& Matrix<T>::operator() (int i, int j) const {return contents[i][j];}

template <typename T>
inline T Matrix<T>::Get(int i, int j) const{ return contents[i][j];}

template <typename T>
inline int Matrix<T>::rows() const { return nRows; }

template <typename T>
inline int Matrix<T>::cols() const { return nCols; }

template <typename T>
inline int Matrix<T>::size() const { return nRows * nCols; }

template <typename T>
inline bool Matrix<T>::empty() const { return (nRows == 0 || nCols == 0); }

template <typename T>
void Matrix<T>::fill(T value)
{
	for(int i=0; i<nRows*nCols;i++)
		contents[0][i]=value;
}


template<typename T>
void Matrix<T>::normalizeMaxMatrix()
{
	T maxML;
	std::list<int> validrows;
	std::list<int> validcols;
	std::list<int>::iterator maxitJ, maxitK;

	for(int i = 0; i < nRows; ++i)
		validrows.push_back(i);

	for(int i = 0; i < nCols; ++i)
		validcols.push_back(i);

	do
	{
		maxML = 0;

		// Find the maximum value of the matrix that has not been explored yet
		for (std::list<int>::iterator j = validrows.begin(); j != validrows.end(); ++j)
			for (std::list<int>::iterator k = validcols.begin(); k != validcols.end(); ++k)
				if(contents[*j][*k] > maxML)
				{
					maxML = contents[*j][*k];
					maxitJ = j;
					maxitK = k;
				}

		// If a value was found, set its whole row and column to 0
		if(maxML > 0)
		{
			for (std::list<int>::iterator l = validrows.begin(); l != validrows.end(); ++l)
				contents[*l][*maxitK]=0;

			for (std::list<int>::iterator l = validcols.begin(); l != validcols.end(); ++l)
				contents[*maxitJ][*l]=0;

			contents[*maxitJ][*maxitK] = maxML;

			validrows.erase(maxitJ);
			validcols.erase(maxitK);
		}

	} while (maxML > 0 && !validrows.empty() && !validcols.empty());
}


template<typename T>
T Matrix<T>::getNormalizedSum() const
{
	T maxML;
	T res = 0;
	std::list<int> validrows;
	std::list<int> validcols;
	std::list<int>::iterator maxitJ, maxitK;

	for(int i = 0; i < nRows; ++i)
		validrows.push_back(i);

	for(int i = 0; i < nCols; ++i)
		validcols.push_back(i);

	do
	{
		maxML = 0;

		// Find the maximum value of the matrix that has not been explored yet
		for (std::list<int>::iterator j = validrows.begin(); j != validrows.end(); ++j)
			for (std::list<int>::iterator k = validcols.begin(); k != validcols.end(); ++k)
				if(maxML < contents[*j][*k])
				{
					maxML = contents[*j][*k];
					maxitJ = j;
					maxitK = k;
				}

		// If a value was found, set its whole row and column to 0
		if(maxML > 0)
		{
			res += maxML;
			validrows.erase(maxitJ);
			validcols.erase(maxitK);
		}

	} while (maxML > 0);

	return res;
}







/**
 * Output operator <<
 * \param out Output stream
 * \param M Matrix to print
 * \return The output stream
 */
template <class T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& M) {

    out << "Matrix size: (" << M.nRows << "," << M.nCols << ")" << std::endl;
    for(int i=0; i<M.nRows;i++){
        for(int j=0; j<M.nCols;j++){
            out << M.contents[i][j] << " ";
        }
        out << std::endl;
    }

    return out;
}

template <typename T>
inline T** Matrix<T>::getPointer() { return contents; }


template <typename T>
void Matrix<T>::readFromCSV(const std::string &filename, bool removelast)
{
	std::ifstream mfile(filename.c_str());
	std::string buffer, token;
	
	// Read first line
	mfile >> buffer;
	
	int numcols = std::count(buffer.begin(), buffer.end(), ',');
	
	if(removelast)
		numcols--;
	
	int numlines = std::count(std::istreambuf_iterator<char>(mfile),
														std::istreambuf_iterator<char>(), '\n');
	
	// Initialize the matrix
	this->resize(numlines, numcols);
	
	// Reset the flags and position of the file
	mfile.clear();
	mfile.seekg(0, std::ios::beg);
	
	std::stringstream ss;
	
	for(int i = 0; i < numlines; ++i)
	{
		mfile >> buffer;
		ss.str(buffer);
		
		for(int j = 0; j < numcols; ++j)
		{
			std::getline(ss, token, ',');
			
			if(token == "NA")
				contents[i][j] = std::numeric_limits<double>::min();
			else
				contents[i][j] = atof(token.c_str());
		}
	}
}



#endif
