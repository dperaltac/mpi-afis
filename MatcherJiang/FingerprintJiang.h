/**
 * \file    FingerprintJiang.h
 * \author  Salvador Garcia <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file for the FingerprintJiang class
 */

#ifndef FINGERPRINTJIANG_H
#define FINGERPRINTJIANG_H

#include "Fingerprint.h"

/**
 * @class FingerprintJiang
 *
 * The Jiang et al. (2000) fingerprint matching algorithm. It outputs the matching score of two
 * fingerprints in the domain [0, 100]. 100 is the perfect matching and 0 is the lowest matching.
 */
class FingerprintJiang : public Fingerprint
{
    public:

        /** Matches the fingerprint with a second one.
         * \param f A second fingerprint
         * \return A float value estimating the score of the match
         */
        float match(const Fingerprint & f) const;

				/** Calls the following methods, in this order:
				 * - computeDistances()
				 * - computeNeighbourhood()
				 * - computeFeatureVectors()
				 */
				void initialize();
				

    protected:
    private:

        Matrix <float> Flk; //!< Matrix storing the values for feature vectors for all minutiae

        /** Computes the similarity value of two feature vectors
         * \param i is the minutiae in the *this fingerprint
         * \param p is the second fingerprint
				 * \param j is the minutiae in \p p
         * \return 
         */
        float similarity (int i, const FingerprintJiang &p, int j) const;

        /** Computes the matching local value between two transformed feature vectors
         * \param Fgi is the first transformed feature vector
         * \param Fgj is the second transformed feature vector
         * \param slij is the similarity value of the two original feature vectors
         * \return the matching local value
         */
        float ML (const float *Fgi, const float *Fgj, float slij) const;

				/** Computes the feature vectors for all minutiae of the fingerprint
				 */

				void computeFeatureVectors();
};

#endif // FINGERPRINTJIANG_H
