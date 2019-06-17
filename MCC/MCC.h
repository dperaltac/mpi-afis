/**
 * \file    MCC.h
 * \author  Salvador Garcia <sglopez@ujaen.es>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Header file for the MCC class
 */

#ifndef MCC_H
#define MCC_H

#include <cmath>

#include "Constants.h"
#include "Fingerprint.h"
#include "ConvexHull.h"
#include "GrahamScanConvexHull.h"
#include "Cylinder.h"
#include "Functions.h"

/**
 * @class MCC
 *
 * The MCC fingerprint matching algorithm. It outputs the matching score of two
 * fingerprints in the domain [0, 100]. 100 is the perfect matching and 0 is the lowest matching.
 */
class MCC : public Fingerprint
{
	public:

		/** Default constructor
			*/
		MCC();

		/** Copy constructor
			* \param fp Other fingerprint
			*/
		MCC(const MCC &fp);

		/** Parameter configuration
			* \param ns Number of cells in each side of the cylinder base
			* \param cons Global consolidation type
			* \param ch Determines whether the convex hull is extracted or not
			* \param fusion Determines the type of fusion
			* \param pbit Flag that determines the use of bit operations
			* \return 0 if there is no error, a different number otherwise
			*/
		static int configureAlgorithm(unsigned int ns, int cons, bool ch, int fus, bool pbit);

		/** Parameter configuration for compatibility with the command line
			* The same parameters as in the previous method are expected:
			* \param ns Number of cells in each side of the cylinder base
			* \param cons Global consolidation type
			* \param ch Determines whether the convex hull is extracted or not
			* \param fusion Determines the type of fusion
			* \return 0 if there is no error, a different number otherwise
			*/
		static int configureAlgorithm(int argc, char *argv[]);

		/** Matches the fingerprint with a second one
			* \param f A second fingerprint
			* \return A float value estimating the score of the match
			*/
		float match(const Fingerprint & f) const;

		/** Adds the structure of the fingerprints to the BD. Must be implemented in the subclass
			*/
		void addStructureBD();

		/** Computes the convex hull from the minutiae set
			* \param convex the set of points of the convex hull (output parameter)
			*/
		void computeConvexHull(vector<point2d> &convex);

		/** Generates and computes the set of cylinders from the set of minutiae
			*/
		void computeCylinders();

		/** Generates and computes the set of cylinders from the set of minutiae
			*/
		void computeCylindersNHS();

		/** After the cylinders are computed, this method removes the non-valid cylinders
			*/
		void compactCylinders();

		/** Loads the DB of cylinders from a binary file
			* \param name Path of the file
			*/
		void loadCylinder(string name);

		/** Loads the DB of cylinders from a binary file
			* \param ifstream input
			*/
		void loadCylinder(ifstream &input);

		/** Stores the DB of cylinders into a binary file
			* \param name Path of the file
			*/
		void writeCylinder(string name);
		
		/**
		* Returns the valid arguments for getopt
		* \returns String of valid arguments for getopt
		*/
		static string getValidArgs();

		/** Calls the following methods, in this order:
			* - computeCylinders()
			* - compactCylinders()
			*/
		virtual void initialize();
		
		//friend operators
		friend std::ostream& operator<<(std::ostream& output, const MCC& F);
		
		
	protected:
		
		enum typeConsolidation {LSS=1, LSA=2, LSSR=3, LSAR=4, LGS=5, LGSR=6, NHS=7};

		// Constants
		static const unsigned int Nd = 6; //!< Number of cylinder sections
		static const float SIGMAD = 0.698131700798; //!< Standard deviation in the directional contribution
		static const float MUP = 20; //!< Sigmoid parameter 1 in the computation of n_P
		static const float TAUP = 0.4; //!< Sigmoid parameter 2 in the computation of n_P
		static const unsigned int MINNP = 4; //!< Minimun number of minutiae in the computation of n_P
		static const unsigned int MAXNP = 12; //!< Maximum number of minutiae in the computation of n_P
		static const float WR = 0.5; //!< Weight parameter in Relaxation computation
		static const float MUP1 = 5; //!< Sigmoid parameter 1 in the computation of d_1 relaxation
		static const float TAUP1 = -1.6; //!< Sigmoid parameter 2 in the computation of d_1 relaxation
		static const float MUP2 = 0.261799387799; //!< Sigmoid parameter 1 in the computation of d_2 relaxation
		static const float TAUP2 = -30; //!< Sigmoid parameter 2 in the computation of d_2 relaxation
		static const float MUP3 = 0.261799387799; //!< Sigmoid parameter 1 in the computation of d_3 relaxation
		static const float TAUP3 = -30; //!< Sigmoid parameter 2 in the computation of d_3 relaxation
		static const unsigned int NREL = 5; //!<Number of relaxation iterations for LSSR and LSAR
		
		static const string valid_args;

		// Constants for performance optimization
		static const float RNEIGHBORHOODRADIUS = 98; //R+NEIGHBORHOODRADIUS;

		// Parameters
		static unsigned int Ns; //!< Number of cells along the cylinder diameter
		static typeConsolidation consolidation; //!< Type of consolidation used in MCC
		static bool convexhull; //!< Determines whether the convex hull is extracted or not
		static bool bit; //!< Flag that defines the use of bit-based operations
		
		vector <Cylinder> cylinders; //!< Set of Cylinders
		
		/**
		* Returns the type of consolidation coded in \a cinput
		* \param cinput String that contains a consolidation type
		* \returns Type of consolidation
		*/
		static typeConsolidation getConsolidationType(const string &cons);

		/** Computes the global matching score from an assignment of cylinders pairs
			* \param elements is the vector of assignment of cylinders
			* \param gamma is the similarity matrix
			* \return the global matching score
			*/
		float computeGlobalScore(const vector <pair <int,int> > & elements, Matrix <float> & gamma) const;

		/** Computes the n_P value for the size of the top similarities
			* \param n_A number of cylinders of the first fingerprint
			* \param n_B number of cylinders of the second fingerprint
			* \return the n_P value
			*/
		static int computeNP(int n_A, int n_B);

		/** Computes the rho value for the relaxation consolidation
			* \param t_a is the first part of the pair of base cylinder
			* \param t_b is the second part of the pair of base cylinder
			* \param k_a is the fisrt part of the pair of the reference cylinder
			* \param k_b is the second part of the pair of reference cylinder
			* \return the rho value for one interation
			*/
		static float rho(Cylinder const & t_a, Cylinder const & t_b, Cylinder const & k_a, Cylinder const & k_b);

/** Computes the rho value for the relaxation consolidation tuned for fvc_ongoing
			* \param t_a is the first part of the pair of base cylinder
			* \param t_b is the second part of the pair of base cylinder
			* \param k_a is the fisrt part of the pair of the reference cylinder
			* \param k_b is the second part of the pair of reference cylinder
			* \return the rho value for one interation
			*/
		static float rho_tun(Cylinder const & t_a, Cylinder const & t_b, Cylinder const & k_a, Cylinder const & k_b);
		
		/** Computes the radial angle of two cylinders
			* \param c1 is the first cylinder
			* \param c2 is the second cylinder
			* \return the radial angle
			*/
		static float dR(Cylinder const & c1, Cylinder const & c2);

		/** Computes the assignment of cylinders by LSS
			* \param gamma is the similarity matrix
			* \param nP is the number of pairs to be found
			* \param nA is the number of cylinders of the first fingerprint
			* \param nB is the number of cylinders of the second fingerprint
			* \return the vector of assignment of cylinders
			*/
		void consolidationLSS(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, unsigned int nA, unsigned int nB) const;

		/** Computes the assignment of cylinders by LSA
			* \param gamma is the similarity matrix
			* \param nP is the number of pairs to be found
			* \param nA is the number of cylindes of the first fingerprint
			* \param nB is the number of cylindes of the second fingerprint
			* \return the vector of assignment of cylinders
			*/
		void consolidationLSA(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, unsigned int nA, unsigned int nB) const;
		
		/** Computes the assignment of cylinders by LGS
			* \param gamma is the similarity matrix
			* \param nP is the number of pairs to be found
			* \param nA is the number of cylinders of the first fingerprint
			* \param nB is the number of cylinders of the second fingerprint
			* \return the vector of assignment of cylinders
			*/
		void consolidationLGS(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, unsigned int nA, unsigned int nB) const;

		/** Computes the assignment of cylinders by LSSR
			* \param gamma is the similarity matrix
			* \param nP is the number of pairs to be found
			* \param cA are the cylinders of the first fingerprint
			* \param cB are the cylinders of the first fingerprint
			* \return the vector of assignment of cylinders
			*/
		void consolidationLSSR(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, vector <Cylinder> const & cA, vector <Cylinder> const & cB) const;

		/** Computes the assignment of cylinders by LSAR
			* \param gamma is the similarity matrix
			* \param nP is the number of pairs to be found
			* \param cA are the cylinders of the first fingerprint
			* \param cB are the cylinders of the first fingerprint
			* \return the vector of assignment of cylinders
			*/
		void consolidationLSAR(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, vector <Cylinder> const & cA, vector <Cylinder> const & cB) const;
		
		/** Computes the assignment of cylinders by LGSR
			* \param gamma is the similarity matrix
			* \param nP is the number of pairs to be found
			* \param cA are the cylinders of the first fingerprint
			* \param cB are the cylinders of the first fingerprint
			* \return the vector of assignment of cylinders
			*/
		void consolidationLGSR(vector <pair <int,int> > &best, Matrix <float> & gamma, unsigned int nP, vector <Cylinder> const & cA, vector <Cylinder> const & cB) const;
		
};


inline int MCC::computeNP(int n_A, int n_B)
{
  return MINNP + roundInt(Cylinder::psi((float)min(n_A,n_B),MUP,TAUP*(MAXNP-MINNP)));
}


inline float MCC::dR(Cylinder const & c1, Cylinder const & c2)
{
    return c1.dFi(c1.getrT(), atan2(static_cast<float>(c1.getY()-c2.getY()), static_cast<float>(c2.getX()-c1.getX())));
}

inline string MCC::getValidArgs() { return valid_args; }

#endif // MCC_H
