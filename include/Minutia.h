/**
 * \file    Minutia.h
 * \author  Joaquin Derrac <jderrac@decsai.ugr.es>
 * \author  Salvador Garcia <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file for the minutia class
 */

#ifndef MINUTIA_H
#define MINUTIA_H

#include <iostream>
#include "Constants.h"

/*! \enum typeMin
 * Enumeration of the type of minutia.
 */
enum typeMin {RIG=1, BIF=2, OTH=3};

/**
 * @class Minutia
 *
 * The minutia class represents a minutia as an index, X Y T components (position and orientation),
 * a quality value, and a type (RIG or BIF).
 */
class Minutia{

    public:
        /** Default constructor */
        Minutia();

        /**
         * Primitive constructor
         * \param index Index of the minutia
         * \param valX  X component
         * \param valY  Y component
         * \param valT  T component
         * \param valQ  Quality
         * \param valM  Minutia type
         */
        Minutia(int index, int valX, int valY, float valT, unsigned int valQ, typeMin valM);

        /** Default destructor */
        virtual ~Minutia();

        /**
         * Copy constructor
         * \param other Object to copy from
         */
        Minutia(const Minutia& other);

        /**
         * Assignment operator
         * \param other Object to assign from
         * \return A reference to this object
         */
        Minutia& operator=(const Minutia& other);

        /**
         * Access index
         * \return Index of the minutia
         */
        int getIndex() const;

        /**
         * Access X
         * \return X value of the minutia
         */
        int getX() const;

        /**
         * Access Y
         * \return Y value of the minutia
         */
        int getY() const;

        /**
         * Access T
         * \return T value of the minutia
         */
        float getT() const;

        /**
         * Access quality
         * \return The quality of the minutia
         */
        unsigned int getQuality() const;

        /**
         * Access type
         * \return The type of the minutia
         */
        typeMin getType() const;

        /**
         * Set index
         * \param val Index to set
         */
        void setIndex(int val);

        /**
         * Set X
         * \param val X value to set
         */
        void setX(int val);

        /**
         * Set Y
         * \param val Y value to set
         */
        void setY(int val);

        /**
         * Set T
         * \param val T value to set
         */
        void setT(float val);

        /**
         * Set quality
         * \param val Quality to set
         */
        void setQuality(unsigned int val);

        /**
         * Set type
         * \param val Type to set
         */
        void setType(typeMin val);

        /**
         * Get Theta (T) measured clockwise
         * \return T value clockwise
         */
        float getcT() const;

        /**
         * Get Theta (T) expressed in radians
         * \return T value in radians
         */
        float getrT() const;

        /**
         * Get Theta (T) normalized in the interval [-180,180)
         * \return T value normalized
         */
        float getnT() const;

        /**
         * Get Theta (T) normalized in the interval [-180,180), and expressed in radians
         * \return T value normalized and expressed in radians
         */
        float getrnT() const;

        /**
         * Get Theta (T) expressed in radians and measured clockwise
         * \return T value in radians and clockwise
         */
        float getcrT() const;

        /**
         * Get Theta (T) normalized in the interval [-180,180) and measured clockwise
         * \return T value normalized and clockwise
         */
        float getcnT() const;

        /**
         * Get Theta (T) normalized in the interval [-180,180), expressed in radians and measured clockwise
         * \return T value normalized, expressed in radians and clockwise
         */
        float getcrnT() const;

        /**
         * Operator ==
         * \param o Minutia to compare
         * \return True if the minutiae are equal. False, otherwise
         */
        bool operator== (const Minutia &o) const;

        /**
         * Operator !=
         * \param o Minutia to compare
         * \return True if the minutiae are not equal. False, otherwise
         */
        bool operator!= (const Minutia &o) const;

        //friend operators
        friend std::ostream& operator<<(std::ostream& output, const Minutia& M);

    protected:
        unsigned int index; //!< Index of the minutia
        int X; //!< X position
        int Y; //!< Y position
        float T; //!< T orientation
    private:
        unsigned int quality; //!< Quality of the minutia
        typeMin type; //!< Type  of the minutia

};

inline int Minutia::getIndex() const { return index; }

inline int Minutia::getX() const { return X; }

inline int Minutia::getY() const { return Y; }

inline float Minutia::getT() const { return T; }

inline unsigned int Minutia::getQuality() const { return quality; }

inline typeMin Minutia::getType() const { return type; }

inline void Minutia::setIndex(int val) { index = val; }

inline void Minutia::setX(int val) { X = val; }

inline void Minutia::setY(int val) { Y = val; }

inline void Minutia::setT(float val) { T = val; }

inline void Minutia::setQuality(unsigned int val) { quality = val; }

inline void Minutia::setType(typeMin val) { type = val; }

inline float Minutia::getrT() const { return (T * REGTORAD); }

inline float Minutia::getnT() const { return T-180; }

inline float Minutia::getrnT() const { return (getnT() * REGTORAD); }

inline float Minutia::getcT() const { return T==0 ? 0 : 360-T; }

inline float Minutia::getcrT() const { return (getcT() * REGTORAD); }

inline float Minutia::getcnT() const { return getcT() - 180; }

inline float Minutia::getcrnT() const { return (getcnT() * REGTORAD); }

inline bool Minutia::operator== (const Minutia &o) const {

    return (X == o.X && Y == o.Y && T == o.T &&
            quality == o.quality && type == o.type);

}

inline bool Minutia::operator!= (const Minutia &o) const {

    return !(*this == o);

}

#endif
