/**
 * \file    Matrix.h
 * \author  Daniel Peralta <daniel.peralta@irc.vib-ugent.be>
 * \version 2.0
 *
 * \section DESCRIPTION_MATRIX
 *
 * A template for building generic matrixes
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <limits>
#include <iostream>
#include <list>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstring>

#include "Functions.h"

template<typename T> class Matrix;
template<typename T> std::ostream& operator<<(std::ostream& out, const Matrix<T>& M);

template <typename T>

/**
 * @class Matrix
 *
 * The Matrix class allows to define generic matrixes, able to store any of the primitive types.
 * It is meant to be more efficient than dynamically reserved vectors of vectors.
 * It supports also some additional functions, such as filling or printing (<<).
 * \todo Maybe turn this into a specialization of Tensor
 */
class Matrix {

public:
	
	
	// iterator
	typedef T* iterator;	///< Iterator
	typedef const T* const_iterator;	///< Constant iterator
	
	// iterator functions
	/**
	 * Iterator pointing to the first element of the matrix
	 * \return Iterator pointing to the first element of the matrix
	 */
	inline iterator begin() { return contents[0]; }
	
	/**
	 * Iterator pointing to the end of the matrix
	 * \return Iterator pointing to the end of the matrix
	 */
	inline iterator end()   { return contents[0] + size(); }
	
	/**
	 * Iterator pointing to the first element of the matrix
	 * \return Iterator pointing to the first element of the matrix
	 */
	inline const_iterator begin() const { return contents[0]; }
	
	/**
	 * Iterator pointing to the end of the matrix
	 * \return Iterator pointing to the end of the matrix
	 */
	inline const_iterator end() const   { return contents[0] + size(); }
	
	/**
	 * Iterator pointing to the beginining of row \p row
	 * \param row Selected row
	 * \return Iterator pointing to row \p row
	 */
	inline iterator begin(int row) { return contents[row]; }
	
	/**
	 * Iterator pointing to the end of row \p row
	 * \param row Selected row
	 * \return Iterator pointing to row \p row
	 */
	inline iterator end(int row)   { return contents[row+1]; }
	
	/**
	 * Iterator pointing to the beginining of row \p row
	 * \param row Selected row
	 * \return Iterator pointing to row \p row
	 */
	inline const_iterator begin(int row) const { return contents[row]; }
	
	/**
	 * Iterator pointing to the end of row \p row
	 * \param row Selected row
	 * \return Iterator pointing to row \p row
	 */
	inline const_iterator end(int row) const   { return contents[row+1]; }

    /** Default Constructor */
	Matrix();

	/**
	 * Constructor
	 * \param rows Rows of the matrix
	 * \param cols Columns of the matrix
	 */
	Matrix(int rows, int cols);

	/**
	 * Constructor
	 * \param rows Rows of the matrix
	 * \param cols Columns of the matrix
	 * \param value Initial value
	 */
	Matrix(int rows, int cols, const T &value);

	/**
	 * Constructor. Builds a matrix out from a vector, reusing the memory.
	 * \param v Vector
	 * \param rows Rows of the matrix
	 * \param cols Columns of the matrix
	 */
	Matrix(T *v, int rows, int cols);

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
     * Changes the size of the matrix.
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
    int minsize() {
        return ((nRows < nCols) ? nRows : nCols);
    }

    /**
     * Fills the matrix with a given value
     * \param value Value used for filling the matrix
     */
    void fill(const T &value);

    /**
		 * Output operator
		 */
		friend std::ostream& operator<< <>(std::ostream& out, const Matrix<T>& M);
		
		/**
		 * Return the pointer of the matrix
		 * \return Pointer to the start of the matrix
		 */
		T** getPointer();
		
		/**
		 * Return the pointer of the matrix, from where all the data can be read as a single vector.
		 * \return Pointer to the start of the matrix
		 */
		T* asVector();
		
		/**
		 * Efficiently swaps the contents of the matrix by that of \p o
		 * \param o Matrix to be swapped by \p this
		 */
		void swap(Matrix &o);
		
		
		/** Transposes the matrix, swapping rows and columns.
		 * Whenever the matrix is square, this is done efficiently by avoiding to re-allocate the memory.
		 * For non-square matrices a new matrix is reserved in memory and the old one is deleted.
		 */
		void transpose();
		
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
    
    int reserved;
};

template<typename T>
Matrix<T>::Matrix() : nRows(0), nCols(0), contents(0), reserved(0) {}

template <typename T>
Matrix<T>::Matrix(int row, int col) : nRows(row), nCols(col), contents(0), reserved(row*col) {

	if (row == 0 || col == 0)
	{
		nRows = 0;
		nCols = 0;
		contents = 0;
		reserved = 0;
	}
	else
	{
		contents = new T* [nRows];
		contents[0]= new T [reserved];

		for (int i=1;i<nRows;++i)
			contents[i] = contents[i-1]+nCols;
	}
 }

template <typename T>
Matrix<T>::Matrix(int row, int col,  const T &value) : nRows(row), nCols(col), contents(0), reserved(row*col) {

	if (row == 0 || col == 0)
	{
		nRows = 0;
		nCols = 0;
		contents = 0;
		reserved = 0;
	}
	else
	{
		contents = new T* [nRows];
		contents[0]= new T [reserved];

		for (int i=1;i<nRows;i++)
		contents[i] = contents[i-1]+nCols;
	}
	
	fill(value);
 }
 
 

template <typename T>
Matrix<T>::Matrix(T *v, int row, int col) : nRows(row), nCols(col), reserved(row*col) {

	if (row == 0 || col == 0)
	{
		nRows = 0;
		nCols = 0;
		contents = 0;
		reserved = 0;
	}
	else
	{
		contents = new T* [nRows];
		contents[0]= v;

		for (int i=1;i<nRows;i++)
			contents[i] = contents[i-1]+nCols;
	}
 }

template <typename T>
Matrix<T>::~Matrix(){

	clear();

}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& o) : nRows(o.nRows), nCols(o.nCols), contents(0), reserved(o.reserved) {

	if (o.contents == 0) {
		nRows = 0;
		nCols = 0;
		contents = 0;
		reserved = 0;
		return;
	}

	contents = new T* [nRows];
	contents[0]= new T [reserved];

	for (int i=1;i<reserved/nCols;++i)
		contents[i] = contents[i-1]+nCols;

	for(int i=0; i<nRows*nCols; i++)
		contents[0][i]=o.contents[0][i];
}



template <typename T>
void Matrix<T>::swap(Matrix<T>& o) {
	
	std::swap(nRows, o.nRows);
	std::swap(nCols, o.nCols);
	std::swap(reserved, o.reserved);
	std::swap(contents, o.contents);
}


template <typename T>
void Matrix<T>::transpose() {
	
	// If the matrix is square, the process is simple
	if(nRows == nCols)
	{
		for(int i=0; i < nRows; i++)
			for(int j = 0; j < i; j++)
				std::swap(contents[i][j], contents[j][i]);
	}
	
	// For a non-square matrix, it's much more complex
	else
	{
		// Create new row pointers
// 		T *newcontents = new T [nCols];
// 		
// 		newcontents[0] = contents[0];
// 		for (int i=1;i<nCols;i++)
// 			newcontents[i] = newcontents[i-1]+nRows;
// 		
// 		// Create binary array to check the cycles
// 		checked = std::bitset<>
// 		
// 		// Transpose data by cycles
// 		while(true)
// 		{
// 			
// 		}
		
		Matrix<T> new_matrix(nCols, nRows);
		
		for(int i=0; i < nRows; i++)
			for(int j = 0; j < nCols; j++)
				new_matrix.contents[j][i] = contents[i][j];
		
		swap(new_matrix);
	}
	
// 	for(int i = 0; i < nRows; ++i)
// 		for(int j = 0; j < nCols; ++j)
// 		{
// 			
// 			T tmp = newcontents[i][j];
// 			newcontents[i][j] = contents[j][i];
// 			contents[j][i] = tmp;
// 		}
	
// 	std::swap(nRows, nCols);
// 	
// 	if(nRows != nCols)
// 		delete [] contents;
// 	
// 	contents = newcontents;
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
	
	if (reserved != o.reserved)
	{
		clear();
		contents = new T* [o.nRows];
		contents[0]= new T [o.reserved];

		for (int i=1; i<o.reserved/o.nCols; i++)
			contents[i] = contents[i-1]+o.nCols;
	}

	nRows = o.nRows;
	nCols = o.nCols;

	for(int i=0; i<nRows*nCols; i++)
		contents[0][i]=o.contents[0][i];
	
	reserved = o.reserved;

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
	reserved = 0;
}

template <typename T>
void Matrix<T>::resize(int rows, int columns)
{
	if(rows == 0 || columns == 0)
	{
		clear();
	}
	else if(columns == nCols && columns*rows <= reserved)
	{
		nRows = rows;
	}
	else if(columns*rows < reserved)
	{
		nRows = rows;
		nCols = columns;
		
		for (int i=1; i < nRows; ++i)
			contents[i] = contents[i-1]+nCols;
	}
	else if(columns*rows > reserved)
	{
		T ** new_matrix;
		
		reserved = rows*columns;
		
		new_matrix   = new T* [rows]; // rows
		new_matrix[0]= new T  [reserved];
		
		for (int i=1; i < rows; ++i)
			new_matrix[i] = new_matrix[i-1]+columns;
		
		if(!empty())
		{
			// copy data from saved pointer to new arrays
			int minrows = std::min<int>(rows, nRows);
			int mincols = std::min<int>(columns, nCols);
			
			for ( int x = 0 ; x < minrows ; x++ )
				memcpy(new_matrix[x], contents[x], mincols * sizeof(T));
			
			clear();
		}
		
		contents = new_matrix;
		
		nRows = rows;
		nCols = columns;
	}
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
inline int Matrix<T>::rows() const { return nRows; }

template <typename T>
inline int Matrix<T>::cols() const { return nCols; }

template <typename T>
inline int Matrix<T>::size() const { return nRows * nCols; }

template <typename T>
inline bool Matrix<T>::empty() const { return (nRows == 0 || nCols == 0); }

template <typename T>
inline void Matrix<T>::fill(const T &value)
{
	memset(contents[0], value, nRows * nCols * sizeof(T));
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
inline T* Matrix<T>::asVector() { return (contents == 0) ? 0 : contents[0]; }

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
