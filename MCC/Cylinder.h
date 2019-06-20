/**
 * \file    Cylinder.h
 * \author  Salvador Garcia <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file for the cylinder class
 */

#ifndef CYLINDER_H
#define CYLINDER_H

#include <vector>
#include <iostream>
#include "Constants.h"
#include "Minutia.h"
#include "Functions.h"
#include "ConvexHull.h"
#include "GrahamScanConvexHull.h"

using namespace std;

/**
 * @class Cylinder
 *
 * The cylinder class represents the cylinder created from a minutia which is extended.
 */
class Cylinder : public Minutia
{

    public:
        /** Default constructor */
        Cylinder();

        /**
         * Primitive constructor
         * \param index Index of the cylinder
         * \param valX  X component
         * \param valY  Y component
         * \param valT  T component
         */
        Cylinder(int index, int valX, int valY, float valT);

        /** Default destructor */
        ~Cylinder();

        /**
         * Copy constructor
         * \param other Object to copy from
         */
        Cylinder(const Cylinder& other);

        /**
         * Transformation constructor from minutia
         * \param other Minutia to copy from
         */
        Cylinder(const Minutia& other);

        /**
         * Assignment operator
         * \param other Object to assign from
         * \return A reference to this object
         */
        Cylinder& operator=(const Cylinder& other);

        /**
         * Operator ==
         * \param o Cylinder to compare
         * \return True if the minutiae are equal. False, otherwise
         */
        bool operator== (const Cylinder &o);

        /**
         * Operator !=
         * \param o Cylinder to compare
         * \return True if the cylinder are not equal. False, otherwise
         */
        bool operator!= (const Cylinder &o);

        //friend operators
        friend std::ostream& operator<<(std::ostream& output, const Cylinder& C);


        /** Parameter configuration
         * \param ns Number of cells in each side of the cylinder base
         * \param pbit Flag that defines the use of bit operations
         * \param nhs Flag that defines the use of the NHS consolidation
				 * \return 0 if there is no error, a different number otherwise
         */
        static int configureAlgorithm(unsigned int ns, bool pbit, bool nhs);

        /** Parameter configuration for compatibility with the command line
         * One parameter is expected, as in the previous method:
         * \param ns Number of cells in each side of the cylinder base
         * \param pbit Flag that defines the use of bit operations
         * \param nhs Flag that defines the use of the NHS consolidation
				 * \return 0 if there is no error, a different number otherwise
         */
        static int configureAlgorithm(int argc, char *argv[]);

        /**
         * Copies the minutia location
         * \param other Minutia to copy from
         */
        void setMinutia(const Minutia& other);

        /**
         * Get the cmVector value at position i
         * \param i index of the cmVector
         * \return the cm[i] value
         */
        float getCM(int i) const;

        /**
         * Get the complete cmVector
         * \return the cm vector
         */
        vector<float> getCMVector() const;

        /**
         * Get the complete cb1Vector
         * \return the b1 vector
         */
        vector<bool> getB1Vector() const;

        /**
         * Get the complete bit mask Vector
         * \return the bit mask vector
         */
        vector<bool> getB2Vector() const;

        /**
         * Get the cmBit1 value at position i
         * \param i index of the cmBit1
         * \return the cm[i] value
         */
        bool getB1(int i) const;

        /**
         * Get the cmBit2 value at position i
         * \param i index of the cmBit2
         * \return the cm[i] value
         */
        bool getB2(int i) const;

        /**
         * Returns the number of cells of the cylinder
         * \return the number of cells of the cylinder
         */
        static unsigned int getNumCells();

        /**
         * Set the cmVector value at position i
         * \param i index of the cmVector
         * \param v is the value to set
         */
        void setCM(int i, float v);

        /**
         * Set the cmBit1 value at position i
         * \param i index of the cmBit1
         * \param v is the value to set
         */
        void setB1(int i, bool v);

        /**
         * Set the cmBit2 value at position i
         * \param i index of the cmBit2
         * \param v is the value to set
         */
        void setB2(int i, bool v);

        /**
         * Set the validity of the cylinder
         * \param v is the validity
         */
        void setValidity(bool v);

        /**
         * Get the validity of the cylinder
         * \return the validity
         */
        bool getValidity() const;

        /** Computes the angle associated to all cells at heigh k
         * \param k is the height of the cell
         * \return A float value estimating the angle associated
         */
        static float dFik(unsigned int k);

        /** Computes the two-dimensional point corresponding to the center of the cell (i,j)
         * \param i i-index
         * \param j j-index
         * \param i_out i-index computed
         * \param j_out j-index computed
         */
        void pij(int i, int j, float & i_out, float & j_out);

        /** Computes the euclidean distance between the cylinder center and a (i,j) point
         * \param i i-index
         * \param j j-index
         * \return the euclidean distance
         */
        float ds(float i, float j) const;

        /** Computes the euclidean distance between the cylinder center and a (i,j) point avoiding the sqrt calculation
         * \param i i-index
         * \param j j-index
         * \return the squared euclidean distance
         */
        float dss(float i, float j) const;

        /** Computes the difference between two angles
         * \param ang1 is the angle 1
         * \param ang2 is the angle 2
         * \return The difference between the two angles in the interval [-PI,PI)
         */
        static float dFi(float ang1, float ang2);

        /** Checks either a cell is valid or not
         * \param i i-index
         * \param j j-index
         * \param convex is the convex hull already computed
         * \return correctness of the cell
         */
        bool xiM(float i, float j, vector <point2d> const & convex);

        /** Checks either a cell is valid or not without convexhull
         * \param i i-index
         * \param j j-index
         * \return correctness of the cell
         */
        bool xiM(float i, float j);

        /** Computes the spatial contribution of the minutia with respect an absolute (i,j) position
         * \param p_i i-index
         * \param p_j j-index
         * \return spatial contribution value
         */
        float contributionS(float p_i, float p_j) const;

        /** Computes the directional contribution of the minutia with respect an angle and height angle
         * \param angle is the angle of comparison (in radians)
         * \param k_angle is the angle associated with k height (in radians)
         * \return directional contribution value
         */
        float contributionD(float angle, float k_angle) const ;
				
				/** Computes the CM vector and sets the validity of the cylinder
         * \param cylinders representes the rest of cylinders involved in the computation
         * \param convex is the convex hull of the minutiae set
				 */
				void computeCMVector(vector <Cylinder> const & cylinders, vector <point2d> const & convex);
				
				/** Computes the CM vector and sets the validity of the cylinder
         * \param cylinders represents the rest of cylinders involved in the computation
				 */
				void computeCMVector(vector <Cylinder> const & cylinders);
				
				/** Computes the CM vector and sets the validity of the cylinder, for NHS consolidation
         * \param cylinders represents the rest of cylinders involved in the computation
				 */
				void computeNHSBitVectors(vector <Cylinder> const & cylinders);

        /** Computes cm value for a cell given by its coordinates (i, j, k)
         * \param i_abs i relative coordinate
         * \param j_abs j relative coordinate
         * \param k k-index
         * \param cylinders represents the rest of cylinders involved in the computation
         * \param convex is the convex hull of the minutiae set
         * \return value of C_m(i,j,k) function
         */
        float cm(float i_abs, float j_abs, unsigned int k, vector <Cylinder> const & cylinders, vector <point2d> const & convex);

        /** Computes cm value for a cell given by its coordinates (i, j, k) without convex hull
         * \param i_abs i relative coordinate
         * \param j_abs j relative coordinate
         * \param k k-index
         * \param cylinders representes the rest of cylinders involved in the computation
         * \return value of C_m(i,j,k) function
         */
        float cm(float i_abs, float j_abs, unsigned int k, vector <Cylinder> const & cylinders);

        /** Sigmoid function that limits the contribution of dense minutiae clusters
         * \param v is the value to be transformed
         * \return Sigmoid transformation of the input
         */
        static float psi(float v, float par1, float par2);

        /** Step function that limits the contribution of dense minutiae clusters (bit implementation)
         * \param v is the value to be transformed
         * \return Step transformation of the input
         */
        static bool psiBit(float v);

        /** Computes the similarity with other cylinder, using the Hamming distance
         * \param c cylinder to compare
         * \return similarity degree
         */
        float nhs(const Cylinder & c) const;

        /** Computes the similarity with other cylinder
         * \param c cylinder to compare
         * \return similarity degree
         */
        float similarity(const Cylinder & c) const;

    protected:
    private:

        // Parameters
        static unsigned int Ns; //!< Number of cells along the cylinder diameter
        static bool bit; //!< Flag that defines the use of bit-based operations

        // Constants
    static const unsigned int R = 70; //!< Cylinder radius (in pixel)
    static const unsigned int Nd = 6; //!< Number of cylinder sections
    static constexpr float SIGMAS = 9.3333333; //!< Standard deviation in the spacial contribution
    static const int OMEGA = 50; //!< Offset applied to enlarge the convex hull (in pixels)
    static constexpr float MINME = 0.6; //!< Minimum number of matching elements in two matchable cylinders
    static constexpr float TAUPSI = 400; //!< Sigmoid parameter 2 for function Psi
    static constexpr float DELTAZETA = 1.57079633; //!< Maximum global rotation allowed between two templates
    static constexpr float DELTAXY = 256; //!< Maximum distance between two templates
    
		static const unsigned int MINM = 4; //!< Minimum number of minutiae for a cylinder to be valid
		static constexpr float RNEIGHBORHOODRADIUS = 98; //R+NEIGHBORHOODRADIUS;
		static constexpr float MINVC = 0.75; //!< Minimum number of valid cells for a cylinder to be valid


    // Constants for performance optimization
    static constexpr float ISIGMAD = 1.43239448783; //1.0/SIGMAD;
    static constexpr float DELTAD = 1.047197551; //(2.0*PI)/Nd;
    static constexpr float GAUSSIANDEN = 0.042743816; //1.0 / (SIGMAS * sqrt(2*PI));
    static constexpr float ISIGMASQ2 = 0.005739796; //1.0 / (2.0*SIGMAS*SIGMAS);
    static constexpr float NEIGHBORHOODRADIUS = 28; //3*SIGMAS;
    
    // Constants depending on parameters
    static unsigned int NUMCELLS;
    static float DELTAS;
    static float RELSHIFT;
    static unsigned int MINCELLS;
    static float MUPSI; //!< Sigmoid parameter 1 for function Psi
    
    static unsigned int MAXNVCELLS;

        bool valid; //!< validity of the cylinder
        vector <float> cmVector; //!< Cylinder codification
        vector <bool> cmBit1; //!< Bit Cylinder codification
        vector <bool> cmBit2; //!< Bit-mask

        /** Checks either a point falls inside of a convex hull
         * \param testx i-index
         * \param testy j-index
         * \param convex is the convex hull already computed
         * \return true if the point is inside the convex polygon
         */
        bool pnpoly(vector <point2d> const & convex, float testx, float testy);
};

inline bool Cylinder::operator== (const Cylinder &o){

    return (X == o.X && Y == o.Y && T == o.T && cmVector == o.cmVector && valid == o.valid && cmBit1 == o.cmBit1 && cmBit2 == o.cmBit2);

}

inline bool Cylinder::operator!= (const Cylinder &o){

    return !(*this == o);
}

inline float Cylinder::getCM(int i) const
{
    return cmVector[i];
}

inline vector<float> Cylinder::getCMVector() const
{
    return cmVector;
}


inline void Cylinder::setCM(int i, float v)
{
    cmVector[i] = v;
}

inline bool Cylinder::getB1(int i) const
{
    return cmBit1[i];
}

inline unsigned int Cylinder::getNumCells()
{
	return NUMCELLS;
}


inline vector<bool> Cylinder::getB1Vector() const
{
    return cmBit1;
}


inline vector<bool> Cylinder::getB2Vector() const
{
    return cmBit2;
}



inline void Cylinder::setB1(int i, bool v)
{
    cmBit1[i] = v;
}

inline bool Cylinder::getB2(int i) const
{
    return cmBit2[i];
}

inline void Cylinder::setB2(int i, bool v)
{
    cmBit2[i] = v;
}

inline void Cylinder::setValidity(bool v)
{
    valid = v;
}

inline bool Cylinder::getValidity() const
{
    return valid;
}

inline float Cylinder::psi(float v, float par1, float par2)
{
	return 1.0 / (1.0 + exp((par2*(par1-v))));
}

inline bool Cylinder::psiBit(float v)
{
		return (v>=MUPSI);
}

inline float Cylinder::dFik(unsigned int k)
{
    return (k - 0.5)*DELTAD - PI;
}

inline float Cylinder::ds(float i, float j) const
{
	return sqrt(dss(i,j));
}

inline float Cylinder::dss(float i, float j) const
{
	return square(this->X - i) + square(this->Y - j);
}

inline float Cylinder::contributionS(float p_i, float p_j) const
{
    return GAUSSIANDEN / exp(dss(p_i,p_j)*ISIGMASQ2);
}

#endif
