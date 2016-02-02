/**
 * \file    Parameters.h
 * \author  Salvador García <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Parameters for the Jiang matching algorithm.
 */

#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

#include "Constants.h"

const int NN = 2; //!< Number of neighboring minutiae considered
const float W1 = 1; //!< Weight factor associated to the distance of minutiae
const float W2 = (0.3*180/PI); //!< Weight factor associated to radial angle
const float W3 = (0.3*180/PI); //!< Weight factor associated to minutiae direction
const float W4 = 3; //!< Weight factor associated to ridge count
const float W5 = 3; //!< Weight factor associated to minutiae types
const int BEST = 5; //!< Number of iterations for the consolidation step
const float BG1 = 8; //!< Tolerance box for distance of minutiae
const float BG2 = PI/6.0; //!< Tolerance box for radial angle
const float BG3 = PI/6.0; //!< Tolerance box for minutiae direction

// Constants for improving performance
const int FVSIZE = 3*NN; //!< Size of the feature vector
const int NEWFVSIZE = 2*NN; //!< Size of the feature vector
const int NN2 = 2*NN;


const float BL = 6*FVSIZE; //!< Threshold for similarity values

#endif
